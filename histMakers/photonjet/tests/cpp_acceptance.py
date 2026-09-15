#!/usr/bin/env python3
"""Plan and run bounded direct-maker acceptance against preserved numeric evidence."""
import argparse
from argparse import Namespace
from copy import deepcopy
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time
from ruamel.yaml import YAML

PJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PJ / 'scripts'))
import inputs
import workflow as flow
from cpp_validation import CRITERIA, compare, root_snapshot

PARTS = {'pp_data': list(range(200)), 'auau_data': list(range(11650, 11850)),
         'pp_signal': [0, 715, 1429], 'pp_inclusive': [0, 715, 2858],
         'auau_signal': [0, 718, 840], 'auau_inclusive': [0, 718, 2150]}
AUDITS = {'pp_signal': (471, '059c0acbaed2bdfaf6a9ced622aa2d065565b3439098064f4c47547639e18136'),
          'auau_signal': (1193, '8d316378a678e21ab32d0e28d365d6af9291ea008dc5554b80d857685bf82168')}
ORIGINAL_FIELDS = ('original_part_index', 'original_part_id', 'source_path', 'source_uuid',
    'source_md5', 'source_bytes', 'events', 'photons', 'truth_photons', 'source_link_rows',
    'truth_jet_rows', 'links_available', 'truth_jets_available')
TIMING_LIMITS = (
    'Bounded inputs and uncontrolled filesystem caches; no full-release throughput or maximum-memory claim. '
    'Each new direct run has an isolated process. module_load_seconds includes PyROOT and C++ initialization; '
    'maker_call_seconds includes input checks, reading, filling and checked publication; analysis_seconds ends '
    'before final input rechecks/output writing. process_wall_seconds additionally includes imports, verification '
    'and snapshot export. Old direct times are separate ROOT processes; old accepted one-worker times are '
    'in-process calls on already converted inputs and exclude conversion. Their reported RSS is a cumulative '
    'process maximum. These different execution boundaries do not establish a production speedup.')


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    inputs.save(path, value)


def read_json(path):
    return json.loads(Path(path).read_text())


def scratch_base(path):
    path = Path(path).resolve()
    if path.parent != Path('/tmp') or not path.name.startswith('ppg12-2-cpp-acceptance-'):
        raise ValueError('acceptance output must be a new /tmp/ppg12-2-cpp-acceptance-* directory')
    return path


def checkout():
    status = subprocess.check_output(['git', 'status', '--porcelain', '--untracked-files=all'], cwd=PJ, text=True)
    if status.strip():
        raise RuntimeError('commit the completed source changes before acceptance planning/running')
    return subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PJ, text=True).strip()


def check_archived_hash(path, archived):
    path = Path(path)
    expected = archived.get(str(path))
    if expected is None:
        aliases = [value for name, value in archived.items() if Path(name).resolve() == path.resolve()]
        if len(set(aliases)) == 1:
            expected = aliases[0]
    if expected is None or digest(path) != expected:
        raise ValueError('missing or changed archived dependency: ' + str(path))


def golden_generation(old, tag):
    """Read an already sealed historical output; no historical backend is loaded."""
    import uproot
    cfg = inputs.config(old['config'])
    paths = flow.output_names(cfg, old['product'], tag)
    marker_path = Path(str(paths['main']) + '.complete.yaml')
    marker = YAML(typ='safe').load(marker_path.read_text())
    if marker['format_version'] != 1 or marker['backend'] != 'photonjet-rdf-v1' or marker['product'] != old['product']:
        raise ValueError('unexpected accepted golden generation')
    membership = {entry['kind']: entry for entry in marker['files']}
    if set(membership) != set(paths) or len(membership) != len(marker['files']):
        raise ValueError('incomplete historical output pair')
    common, provenance = None, None
    for kind, path in paths.items():
        entry = membership[kind]
        if Path(entry['path']).resolve() != path.resolve() or path.stat().st_size != entry['bytes'] or inputs.md5(path) != entry['md5']:
            raise ValueError('historical output differs from its completion seal: ' + str(path))
        with uproot.open(path) as file:
            if str(file['config']) != Path(old['config']).read_text():
                raise ValueError('historical output/config mismatch')
            text = str(file['photonjet_provenance'])
        if common is not None and common != text:
            raise ValueError('mixed historical main/response generation')
        common = text
        provenance = YAML(typ='safe').load(text)
        if (provenance['run_id'] != marker['run_id'] or provenance['workers'] != 1 or provenance['graph_runs'] != 1 or
                provenance['format_version'] != 2 or provenance['backend'] != marker['backend'] or
                provenance['schema_version'] != 1 or provenance['rng_version'] != 'fnv1a64-splitmix64-mixmax17-gaus-v1'):
            raise ValueError('historical baseline is not the accepted schema/RNG and one worker/pass')
        ids = [record['original_part_id'] for record in provenance['inputs']]
        if len(ids) != len(set(ids)) or sorted(ids) != sorted(marker['original_part_ids']):
            raise ValueError('historical original coverage differs from its completion seal')
    return ({kind: str(path) for kind, path in paths.items()}, marker_path, provenance)


def copy_config(source, directory, variation, no_smear=False):
    cfg = inputs.config(source)
    cfg['output']['var_type'] = variation
    for key, stem in (('data_outfile', 'data'), ('eff_outfile', 'eff'),
                      ('response_outfile', 'response'), ('final_outfile', 'final')):
        cfg['output'][key] = str(directory / stem)
    cfg['photonjet']['external_mbd_eff_file'] = ''
    if no_smear:
        for kind in ('data', 'mc'):
            for index in range(3):
                cfg['analysis'][f'cluster_eres_{kind}_p{index}'] = 0.
    directory.mkdir(parents=True, exist_ok=False)
    path = directory / 'config.yaml'
    with path.open('w') as stream:
        YAML().dump(cfg, stream)
    return path, cfg


def selected_originals(cfg, product, info, known):
    """Seal the actual original files, matching historical identities and counts."""
    master, paths = inputs.master_parts(cfg, product)
    selected = []
    for old in info['inputs']:
        index = old['original_part_index']
        path = Path(paths[index])
        key = str(path.resolve())
        if key not in known:
            current = inputs.inspect_source(path, index, cfg, product, info['release_id'])
            known[key] = dict(record=current, sha256=digest(path))
        current = known[key]['record']
        for field in ORIGINAL_FIELDS:
            if field == 'source_path':
                if Path(current[field]).resolve() != Path(old[field]).resolve():
                    raise ValueError('original path differs from historical coverage')
            elif current[field] != old[field]:
                raise ValueError('original identity/count differs: ' + field + ' ' + str(path))
        if current['original_part_index'] != index or current['original_part_id'] != old['original_part_id']:
            raise ValueError('one original path has conflicting canonical identities')
        selected.append(deepcopy(current))
    return master, selected


def write_list(path, records):
    path.write_text(''.join(f"{record['original_part_index']} {record['source_path']}\n" for record in records))


def plan(base, golden, reference):
    head = checkout()
    base = scratch_base(base)
    base.mkdir(parents=True, exist_ok=False)
    old_plan = read_json(golden / 'plan.json')
    old_results = read_json(golden / 'results.json')
    first_plan, first_results = read_json(reference / 'plan.json'), read_json(reference / 'results.json')
    if not old_results.get('passed') or not first_results.get('passed'):
        raise ValueError('historical acceptance did not pass')
    historical_commands = read_json(reference / 'commands.json')
    files = {golden / 'plan.json', golden / 'results.json', reference / 'plan.json',
             reference / 'results.json', reference / 'commands.json'}
    cases, known = [], {}
    first_cases = {item['name']: item for item in first_plan['cases']}
    first_outcomes = {item['name']: item for item in first_results['cases']}
    for category in ('reference', 'matrix', 'scaling'):
        for old in old_plan[category]:
            check_archived_hash(old['config'], old_plan['sha256'])
            tag = 'reference' if category == 'reference' else 'workers1' if category == 'matrix' else 'first1'
            golden_paths, marker, info = golden_generation(old, tag)
            files.update((Path(old['config']), marker, *(Path(path) for path in golden_paths.values())))
            mode = 'no_smear' if category == 'reference' or old['name'].endswith('_no_smear') else 'full'
            if category == 'matrix':
                stem = old['name'].removesuffix('_' + old['product'] + '_' + mode)
                system = info['system']
                source = PJ / 'configs' / system / (stem + '.yaml')
                check_archived_hash(source, old_plan['sha256'])
            elif category == 'scaling':
                source = PJ / 'configs' / info['system'] / ('config_pj_pp_nom.yaml' if info['system'] == 'pp' else 'config_pj_auau_c00_20.yaml')
                check_archived_hash(source, old_plan['sha256'])
                if old['indices'] != PARTS[old['name']]:
                    raise ValueError('historical scaling selection changed')
            else:
                source = Path(old['config'])
            directory = base / category / old['name']
            cfg_path, cfg = copy_config(source, directory, 'PPG12_2_cpp_' + category + '_' + old['name'], mode == 'no_smear')
            master, records = selected_originals(cfg, old['product'], info, known)
            listing = directory / 'originals.list'
            write_list(listing, records)
            item = dict(category=category, name=old['name'], product=old['product'], mode=mode,
                config=str(cfg_path), source_config=str(source), input=str(listing), indices=[p['original_part_index'] for p in records],
                originals=records, release_id=info['release_id'], golden=golden_paths)
            files.update((source, cfg_path, master, listing))
            if category == 'reference':
                snapshot = Path(old['snapshot'])
                check_archived_hash(snapshot, old_plan['sha256'])
                item['reference_snapshot'] = str(snapshot)
                item['legacy'] = dict(zip(golden_paths, first_cases[old['name']]['reference']['outputs']))
                for path in item['legacy'].values():
                    if digest(path) != first_outcomes[old['name']]['outputs_sha256'][path]:
                        raise ValueError('frozen reference ROOT changed: ' + path)
                files.update((snapshot, *(Path(path) for path in item['legacy'].values())))
            elif category == 'matrix' and mode == 'no_smear':
                old_cfg = inputs.config(old['config'])
                item['legacy'] = {kind: str(path) for kind, path in flow.output_names(old_cfg, old['product'], 'legacy').items()}
                files.update(Path(path) for path in item['legacy'].values())
                successful = [r for r in old_results['comparisons'] if r['name'] == old['name'] + ' frozen maker']
                if len(successful) != 1 or successful[0]['differences']:
                    raise ValueError('missing accepted frozen-maker comparison')
            elif category == 'scaling':
                if item['indices'] != PARTS[item['name']]:
                    raise ValueError('golden scaling ROOT coverage differs from the planned indices')
                reverse = directory / 'reversed.list'
                write_list(reverse, list(reversed(records)))
                item['reversed_input'] = str(reverse)
                files.add(reverse)
                item['partitions'] = []
                for size in ((50, 73) if item['product'] == 'data' else (1, 2)):
                    partition_cfg, _ = copy_config(cfg_path, directory / f'partition_{size}',
                        'PPG12_2_cpp_' + item['name'] + '_partition_' + str(size))
                    item['partitions'].append(dict(config=str(partition_cfg), size=size))
                    files.add(partition_cfg)
                if item['product'] == 'sim_signal':
                    # This JSON is preserved evidence only; no converted input is consumed.
                    audit = PJ / '.analysis-context/PPG12-2-rdf-completion-20260915' / ('acceptance-a5-' + item['name'] + '-response-audit.json')
                    rows = read_json(audit)
                    count, checksum = AUDITS[item['name']]
                    if len(rows) != count or hashlib.sha256('\n'.join(rows).encode()).hexdigest() != checksum:
                        raise ValueError('accepted response audit changed')
                    item['golden_audit'] = str(audit)
                    files.add(audit)
            files.update(Path(path) for path in inputs.dependencies(cfg, item['product']))
            cases.append(item)
    scalar_cases = {(Path(item['source_config']).resolve(), item['product']): item for item in cases
                    if item['category'] == 'matrix' and item['mode'] == 'no_smear'}
    for item in cases:
        source = item if item['category'] == 'reference' else scalar_cases[
            (Path(item['source_config']).resolve(), item['product'])]
        # Full-mode RDF values remain the numeric baseline. Titles retain the
        # frozen scalar-maker contract, including TEfficiency's inner histograms.
        item['title_reference'] = dict(category=source['category'], name=source['name'],
            files=deepcopy(source['legacy']),
            sha256={path: digest(path) for path in source['legacy'].values()})
    matrix = [item for item in cases if item['category'] == 'matrix']
    expected = {(path.stem, product, mode) for system in ('pp', 'auau')
                for path in (PJ / 'configs' / system).glob('*.yaml')
                for product in inputs.PRODUCTS for mode in ('no_smear', 'full')}
    observed = {(item['name'].removesuffix('_' + item['product'] + '_' + item['mode']), item['product'], item['mode']) for item in matrix}
    if len(matrix) != 36 or observed != expected:
        raise ValueError('acceptance must cover all six shipped configs, three products and both smearing modes')
    scaling_records = [record for item in cases if item['category'] == 'scaling' for record in item['originals']]
    if len(scaling_records) != 412:
        raise ValueError('bounded scaling coverage must contain 412 original parts')
    totals = tuple(sum(record[key] for record in scaling_records) for key in ('events', 'photons', 'truth_photons'))
    if totals != (164254, 60135, 43561):
        raise ValueError('bounded original counts differ from accepted 412-part coverage')
    code = [*(PJ / name for name in inputs.CODE_SOURCES), PJ / 'scripts/inputs.py', PJ / 'scripts/workflow.py',
            PJ / 'scripts/make_chunks.py', PJ / 'scripts/CheckChunks.C', PJ / 'histmakers/InjectMbdEff.C',
            PJ / 'tests/cpp_validation.py', Path(__file__)]
    files.update(code)
    timings = dict(original_direct=[row for row in historical_commands if row['label'] == 'reference'],
                   accepted_one_worker=[row for row in old_results['graphs'] if row['workers'] == 1],
                   limits=TIMING_LIMITS)
    booking_config = inputs.config(PJ / 'configs/auau/config_pj_auau_c00_20.yaml')['analysis']
    save(base / 'plan.json', dict(format_version=1, git_head=head, criteria=CRITERIA, cases=cases,
        booking_axes={key: booking_config[key] for key in ('pT_bins', 'pT_bins_truth')},
        originals=known, sha256={str(path): digest(path) for path in sorted(files)}, historical_timings=timings,
        scaling_coverage={key: sum(record[key] for record in scaling_records)
                          for key in ('events', 'photons', 'truth_photons', 'source_bytes')},
        scope='Six reference cases, 36 config/product/mode cases and 412-part bounded order/chunk acceptance'))
    print(base / 'plan.json', flush=True)


def snapshot_files(paths):
    return {kind: root_snapshot(path, titles=True)['objects'] for kind, path in paths.items()}


def snapshot_save(path, value):
    with gzip.open(path, 'wt') as stream:
        json.dump(value, stream, allow_nan=False)


def snapshot_read(path):
    with gzip.open(path, 'rt') as stream:
        return json.load(stream)


def without_titles(value):
    if isinstance(value, dict):
        return {key: without_titles(item) for key, item in value.items() if key != 'title'}
    if isinstance(value, list):
        return [without_titles(item) for item in value]
    return value


def titles(value, prefix=''):
    result = {}
    if isinstance(value, dict):
        for key, item in value.items():
            path = prefix + '/' + key
            if key == 'title':
                result[path] = item
            else:
                result.update(titles(item, path))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            result.update(titles(item, prefix + '/' + str(index)))
    return result


def equal(expected, actual, results, name, compare_titles=True, expected_titles=None):
    differences = compare(without_titles(expected), without_titles(actual))
    title_differences = compare(titles(expected) if expected_titles is None else expected_titles,
                                titles(actual)) if compare_titles else []
    result = dict(name=name, numeric_differences=len(differences), title_differences=len(title_differences),
                  examples=(differences + title_differences)[:10])
    results['comparisons'].append(result)
    if differences or title_differences:
        raise AssertionError((name, result['examples']))



BOOKING_PROBE = r"""
#include <fstream>
#include <set>
#include <sys/resource.h>
#include <unistd.h>
std::string PJCppBookingProbe(unsigned cells, const std::vector<double> &pt, const std::vector<double> &truth_pt) {
    if(cells!=1 && cells!=6) throw std::runtime_error("booking probe requires one or six cells");
    PJ::Cuts cuts; cuts.system="auau"; cuts.pT_bins=pt; cuts.pT_bins_truth=truth_pt;
    cuts.centrality_bins = cells==1 ? std::vector<double>{0,80} : std::vector<double>{0,20,50,80};
    cuts.eta_bins = cells==1 ? std::vector<double>{-.7,.7} : std::vector<double>{-.7,0,.7};
    const PJ::BinLayout bins(cuts);
    auto rss=[]() -> ULong64_t {
        std::ifstream stream("/proc/self/statm"); ULong64_t total=0,resident=0;
        if(!(stream>>total>>resident)) throw std::runtime_error("cannot read resident memory");
        return resident*static_cast<ULong64_t>(sysconf(_SC_PAGESIZE))/1024;
    };
    auto peak=[]() -> ULong64_t {
        struct rusage usage{};
        if(getrusage(RUSAGE_SELF,&usage)) throw std::runtime_error("cannot read peak memory");
        return usage.ru_maxrss;
    };
    YAML::Node result; result["layout"]=bins.Metadata(); result["cells"]=cells;
    result["resident_before_kib"]=rss(); result["peak_before_kib"]=peak();
    auto qa=PJ::BookJobQA(cuts);
    auto histograms=PJ::BookHistograms(cuts,bins,true);
    result["resident_after_kib"]=rss(); result["peak_after_kib"]=peak();
    std::set<const TH1 *> seen;
    auto count=[&](YAML::Node group,const TH1 *hist) {
        if(!hist || !seen.insert(hist).second) throw std::runtime_error("missing or shared histogram buffer");
        if(!dynamic_cast<const TH1D *>(hist) && !dynamic_cast<const TH2D *>(hist))
            throw std::runtime_error("booking probe only accounts double histogram storage");
        group["histograms"]=group["histograms"].as<ULong64_t>()+1;
        group["content_bytes"]=group["content_bytes"].as<ULong64_t>()+sizeof(double)*hist->GetNcells();
        group["sumw2_bytes"]=group["sumw2_bytes"].as<ULong64_t>()+sizeof(double)*hist->GetSumw2N();
    };
    for(const char *name:{"global_buffers","cell_buffers"})
        for(const char *key:{"histograms","content_bytes","sumw2_bytes"}) result[name][key]=ULong64_t(0);
    for(auto *h:{qa.cutflow.get(),qa.flags.get(),qa.vertex.get(),qa.weight.get(),qa.centrality.get(),qa.run.get(),qa.recipe.get()})
        count(result["global_buffers"],h);
    for(const auto &row:histograms) for(const auto &h:row) {
        auto group=result["cell_buffers"];
        for(auto *hist:{h.reco.all.get(),h.reco.common.get(),h.reco.tight.get(),h.reco.signal_all.get(),h.reco.signal_tight.get(),
                       h.truth.spectrum.get(),h.truth.novtx.get(),h.truth.vertexcut.get(),h.truth.mbd.get(),h.truth.north.get(),
                       h.truth.south.get(),h.truth.only_north.get(),h.truth.only_south.get(),h.truth.neither.get(),
                       h.response.reco.get(),h.response.truth.get()}) count(group,hist);
        for(int i=0;i<4;++i)
            for(auto *hist:{h.reco.abcd[i].get(),h.reco.signal[i].get(),h.reco.unmatched[i].get()}) count(group,hist);
        for(auto *hist:{h.qa.isolation.get(),h.qa.score.get(),h.response.matrix.get()}) count(group,hist);
        for(auto *eff:{h.truth.reco.get(),h.truth.iso.get(),h.truth.id.get(),h.truth.all.get(),h.truth.converts.get()}) {
            count(group,eff->GetPassedHistogram()); count(group,eff->GetTotalHistogram());
        }
        for(const auto *hist:std::array<const TH1 *,4>{h.response.unfold->Hresponse(),h.response.unfold->Hmeasured(),
                       h.response.unfold->Htruth(),h.response.unfold->Hfakes()}) count(group,hist);
    }
    result["booked_main_objects"]=7+33*cells; result["booked_response_objects"]=4*cells;
    return YAML::Dump(result);
}
"""


def child_booking(job_path):
    """Book empty signal histograms; read no original trees and write no ROOT file."""
    job = read_json(job_path)
    if Path(job['result']).exists():
        raise ValueError('refusing to replace booking evidence')
    ROOT = flow.root_module('histmakers/PhotonJetHistograms.h')
    if not ROOT.gInterpreter.Declare(BOOKING_PROBE):
        raise RuntimeError('cannot load booking memory probe')
    record = YAML(typ='safe').load(str(ROOT.PJCppBookingProbe(job['cells'], job['pT_bins'], job['pT_bins_truth'])))
    record.update(ROOT_version=str(ROOT.gROOT.GetVersion()), axes={key: job[key] for key in ('pT_bins', 'pT_bins_truth')},
        limits='Empty booking only. RSS includes allocator/ROOT effects; buffers count double contents and Sumw2, '
               'including TEff passed/total and Roo internals, excluding axes, labels, object headers and lazy fit caches. '
               'Global QA is counted once; input metadata and event rows are not allocated.')
    record['resident_delta_kib'] = record['resident_after_kib'] - record['resident_before_kib']
    record['peak_delta_kib'] = record['peak_after_kib'] - record['peak_before_kib']
    save(job['result'], record)


def booking_resources(base, plan_data, results):
    for cells in (1, 6):
        directory = base / 'booking' / str(cells)
        directory.mkdir(parents=True, exist_ok=False)
        job = dict(plan_data['booking_axes'], cells=cells, result=str(directory / 'result.json'))
        path = directory / 'job.json'
        save(path, job)
        timing, log = directory / 'time.json', directory / 'run.log'
        command = ['/usr/bin/time', '-f', '{"process_wall_seconds":%e,"peak_rss_kib":%M}', '-o', str(timing),
                   sys.executable, str(Path(__file__).resolve()), 'booking', '--job', str(path)]
        with log.open('w') as stream:
            process = subprocess.run(command, cwd=PJ, stdout=stream, stderr=subprocess.STDOUT)
        results['commands'].append(dict(command=command, cwd=str(PJ), log=str(log), returncode=process.returncode))
        if process.returncode:
            raise RuntimeError('booking resource probe failed; see ' + str(log))
        record = read_json(job['result'])
        record.update(read_json(timing), result_path=job['result'], result_sha256=digest(job['result']))
        results['booking'].append(record)
        print('[cpp acceptance] booking', cells, record['resident_delta_kib'], 'KiB resident delta', record['cell_buffers'], flush=True)
    one, six = results['booking']
    if one['global_buffers'] != six['global_buffers'] or any(six['cell_buffers'][key] != 6 * value for key, value in one['cell_buffers'].items()):
        raise AssertionError('histogram buffer storage does not scale with the configured cell count')

def child_case(job_path):
    """Internal isolated worker, invoked by run after its complete preflight."""
    job = read_json(job_path)
    os.environ['PJ_RELEASE_ID'] = job['release_id']
    start = time.perf_counter()
    ROOT = flow.root_module('histmakers/PhotonJetHistMaker.C')
    module_seconds = time.perf_counter() - start
    paths = flow.output_names(inputs.config(job['config']), job['product'], job['tag'])
    if any(path.exists() for path in [*paths.values(), Path(str(paths['main']) + '.complete.yaml')]):
        raise RuntimeError('refusing to replace existing acceptance outputs')
    started = time.perf_counter()
    ROOT.PJ.Run(job['config'], job['product'], job['input'], job['tag'], job.get('audit', ''))
    maker_seconds = time.perf_counter() - started
    tick = time.perf_counter()
    info = flow.checked_generation(ROOT, paths['main'], job['config'], job['product'])
    verification_seconds = time.perf_counter() - tick
    if info['workers'] != 1 or info['event_loops'] != 1 or info['events_processed'] != sum(p['events'] for p in info['inputs']):
        raise AssertionError('direct serial accounting differs from actual originals')
    tick = time.perf_counter()
    snapshot = snapshot_files(paths)
    snapshot_save(job['snapshot'], snapshot)
    record = dict(category=job['category'], name=job['name'], tag=job['tag'], module_load_seconds=module_seconds,
        maker_call_seconds=maker_seconds, verification_seconds=verification_seconds,
        snapshot_seconds=time.perf_counter() - tick, process_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        ROOT_version=str(ROOT.gROOT.GetVersion()), snapshot=job['snapshot'], snapshot_sha256=digest(job['snapshot']),
        output_paths={kind: str(path) for kind, path in paths.items()},
        output_sha256={str(path): digest(path) for path in paths.values()},
        provenance={key: info[key] for key in ('preflight_seconds', 'loading_seconds', 'processing_seconds',
            'analysis_seconds', 'events_processed', 'parts_processed', 'event_loops', 'workers', 'selected_events_by_centrality')})
    save(job['result'], record)


def isolated(item, tag, results, input_path=None, audit=False):
    directory = Path(item['config']).parent / 'runs' / tag
    directory.mkdir(parents=True, exist_ok=False)
    job = {key: item[key] for key in ('category', 'name', 'config', 'product', 'release_id')}
    job.update(input=str(input_path or item['input']), tag=tag, snapshot=str(directory / 'snapshot.json.gz'),
               result=str(directory / 'result.json'))
    if audit:
        job['audit'] = str(directory / 'response-audit.txt')
    path = directory / 'job.json'
    save(path, job)
    timing_path, log = directory / 'time.json', directory / 'run.log'
    command = ['/usr/bin/time', '-f', '{"process_wall_seconds":%e,"peak_rss_kib":%M}', '-o', str(timing_path),
               sys.executable, str(Path(__file__).resolve()), 'case', '--job', str(path)]
    started = time.perf_counter()
    with log.open('w') as stream:
        process = subprocess.run(command, cwd=PJ, stdout=stream, stderr=subprocess.STDOUT)
    call = dict(command=command, cwd=str(PJ), log=str(log), returncode=process.returncode,
                seconds=time.perf_counter() - started)
    results['commands'].append(call)
    if process.returncode:
        raise RuntimeError('direct acceptance worker failed; see ' + str(log))
    timing = json.loads(timing_path.read_text().splitlines()[-1])
    record = read_json(job['result'])
    record.update(process_wall_seconds=timing['process_wall_seconds'], isolated_peak_rss_kib=timing['peak_rss_kib'])
    if audit:
        record['audit'] = job['audit']
    results['runs'].append(record)
    print('[cpp acceptance]', item['category'], item['name'], tag, record['provenance'], flush=True)
    return snapshot_read(job['snapshot']), record


def verify_plan(plan_data):
    if checkout() != plan_data['git_head']:
        raise ValueError('checkout changed after acceptance planning')
    for path, checksum in plan_data['sha256'].items():
        if digest(path) != checksum:
            raise ValueError('planned dependency changed: ' + path)
    for item in plan_data['originals'].values():
        path = item['record']['source_path']
        if digest(path) != item['sha256']:
            raise ValueError('planned original input changed: ' + path)
    cases = {(item['category'], item['name']): item for item in plan_data['cases']}
    for item in plan_data['cases']:
        reference = item['title_reference']
        source = cases[(reference['category'], reference['name'])]
        if (source['mode'] != 'no_smear' or source['product'] != item['product'] or
                Path(source['source_config']).resolve() != Path(item['source_config']).resolve() or
                reference['files'] != source['legacy'] or set(reference['files']) != set(item['golden']) or
                reference['sha256'] != {path: plan_data['sha256'][path] for path in source['legacy'].values()}):
            raise ValueError('title reference differs from its sealed scalar source case')


def audit_equal(expected, paths, results, name):
    actual = sorted(line for path in paths for line in Path(path).read_text().splitlines())
    if actual != expected or len(actual) != len(set(actual)):
        raise AssertionError('exact candidate identity/seed/Gaussian/reco/truth/weight audit differs: ' + name)
    results['audits'].append(dict(name=name, matched_candidates=len(actual),
        sha256=hashlib.sha256('\n'.join(actual).encode()).hexdigest(), paths=[str(path) for path in paths]))


def partitions(base, item, baseline, results):
    from make_chunks import make
    ROOT = flow.root_module('histmakers/PhotonJetHistMaker.C', 'scripts/CheckChunks.C')
    for partition in item['partitions']:
        cfg = inputs.config(partition['config'])
        chunk_dir = base / 'jobs'
        arguments = Namespace(config=Path(partition['config']), products=item['product'],
            chunk_size=item['product'] + '=' + str(partition['size']), outdir=chunk_dir,
            release_id=item['release_id'], indices=','.join(map(str, item['indices'])))
        tick = time.perf_counter()
        job_plan = make(arguments)
        results['commands'].append(dict(kind='python_function', function='make_chunks.make',
            arguments={key: str(value) if isinstance(value, Path) else value for key, value in vars(arguments).items()},
            cwd=str(PJ), completed=True))
        planning_seconds = time.perf_counter() - tick
        args = Namespace(config=arguments.config, chunk_dir=chunk_dir, products=item['product'], threads=1)
        tick = time.perf_counter()
        flow.local(args)
        results['commands'].append(dict(kind='python_function', function='workflow.local',
            arguments={key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
            cwd=str(PJ), completed=True))
        processing_seconds = time.perf_counter() - tick
        tick = time.perf_counter()
        flow.merge(args)
        results['commands'].append(dict(kind='python_function', function='workflow.merge',
            arguments={key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
            cwd=str(PJ), completed=True))
        merge_seconds = time.perf_counter() - tick
        paths = flow.output_names(cfg, item['product'])
        equal(baseline, snapshot_files(paths), results, item['name'] + ' checked merge size ' + str(partition['size']))
        info = flow.checked_generation(ROOT, paths['main'], partition['config'], item['product'])
        if info['parts_processed'] != len(item['indices']) or info['events_processed'] != sum(p['events'] for p in item['originals']):
            raise AssertionError('repartitioned input accounting differs')
        plan_path = chunk_dir / cfg['output']['var_type'] / 'jobs.json'
        results['partitions'].append(dict(name=item['name'], chunk_size=partition['size'],
            jobs=len(job_plan['products'][item['product']]['jobs']), plan=str(plan_path), plan_sha256=digest(plan_path),
            planning_seconds=planning_seconds, local_seconds=processing_seconds, merge_seconds=merge_seconds,
            output_sha256={str(path): digest(path) for path in paths.values()}))
        if item['product'] == 'sim_signal' and partition['size'] == 2:
            # Audit actual response decisions with these same changed boundaries.
            # The audit copies use new tags/lists and leave sealed worker outputs intact.
            audit_paths = []
            for index, job in enumerate(job_plan['products'][item['product']]['jobs']):
                listing = Path(item['config']).parent / f'audit_partition_{index}.list'
                listing.write_text(Path(job['input']).read_text())
                _, record = isolated(item, 'audit_partition_' + str(index), results, listing, audit=True)
                audit_paths.append(record['audit'])
            audit_equal(read_json(item['golden_audit']), audit_paths, results, item['name'] + ' regrouped candidates')


def run(base):
    base = scratch_base(base)
    plan_data = read_json(base / 'plan.json')
    if (base / 'results.json').exists():
        raise ValueError('preserve this acceptance attempt; use a new output directory for a retry')
    preflight = time.perf_counter()
    verify_plan(plan_data)
    results = dict(passed=False, criteria=CRITERIA, preflight_seconds=time.perf_counter() - preflight,
        environment=dict(hostname=platform.node(), platform=platform.platform(), python=sys.version,
                         cpu_affinity=sorted(os.sched_getaffinity(0)), git_head=plan_data['git_head']),
        runs=[], comparisons=[], partitions=[], audits=[], booking=[], commands=[], historical_timings=plan_data['historical_timings'],
        scaling_coverage=plan_data['scaling_coverage'], timing_limits=TIMING_LIMITS)
    started = time.perf_counter()
    try:
        booking_resources(base, plan_data, results)
        cases = {(item['category'], item['name']): item for item in plan_data['cases']}
        title_baselines = {}
        for item in plan_data['cases']:
            golden = snapshot_files(item['golden'])
            if 'legacy' in item:
                legacy = snapshot_files(item['legacy'])
                # A saved old-maker file has no completion seal. Cross-check all
                # values with the sealed accepted baseline before using it.
                equal(legacy, golden, results, item['name'] + ' historical baseline cross-check', compare_titles=False)
                expected = legacy
            else:
                expected = golden
            reference = item['title_reference']
            source_key = (reference['category'], reference['name'])
            if source_key not in title_baselines:
                source = cases[source_key]
                if source is item:
                    scalar = legacy
                else:
                    scalar = snapshot_files(reference['files'])
                    # Check the title source against its own sealed no-smear
                    # generation, even for full-mode or larger-input cases.
                    equal(scalar, snapshot_files(source['golden']), results,
                          item['name'] + ' scalar title-source cross-check', compare_titles=False)
                title_baselines[source_key] = titles(scalar)
            actual, record = isolated(item, 'direct', results, audit=item['category'] == 'scaling' and item['product'] == 'sim_signal')
            equal(expected, actual, results, item['category'] + ' ' + item['name'] + ' golden',
                  expected_titles=title_baselines[source_key])
            if 'reference_snapshot' in item:
                frozen = read_json(item['reference_snapshot'])
                snapshots = {kind: value['objects'] for kind, value in zip(actual, frozen)}
                equal(snapshots, actual, results, item['name'] + ' frozen snapshot', compare_titles=False)
            if item['category'] == 'scaling':
                warm, _ = isolated(item, 'warm', results)
                equal(actual, warm, results, item['name'] + ' repeated direct run')
                reverse, reversed_record = isolated(item, 'reversed', results, item['reversed_input'], audit=item['product'] == 'sim_signal')
                equal(actual, reverse, results, item['name'] + ' reversed originals')
                if item['product'] == 'sim_signal':
                    audit = read_json(item['golden_audit'])
                    audit_equal(audit, [record['audit']], results, item['name'] + ' accepted candidates')
                    audit_equal(audit, [reversed_record['audit']], results, item['name'] + ' reversed candidates')
                partitions(base, item, actual, results)
            save(base / 'results.json', results)
        verify_plan(plan_data)
        results['passed'] = True
    except BaseException as error:
        results['error'] = repr(error)
        raise
    finally:
        results['total_seconds'] = time.perf_counter() - started
        results['parent_peak_rss_kib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        save(base / 'results.json', results)
        save(base / 'commands.json', results['commands'])
        save(base / 'timings.json', dict(limits=TIMING_LIMITS, current=[{key: value for key, value in row.items()
            if key not in ('output_sha256', 'output_paths', 'snapshot', 'snapshot_sha256')} for row in results['runs']],
            historical=plan_data['historical_timings'], partitions=results['partitions'], booking=results['booking']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('plan', 'run', 'case', 'booking'))
    parser.add_argument('--outdir', type=Path)
    parser.add_argument('--golden', type=Path, default=Path('/tmp/ppg12-2-rdf-acceptance-20260915-a5'))
    parser.add_argument('--reference', type=Path, default=Path('/tmp/ppg12-2-rdf-phase1-20260914-a2'))
    parser.add_argument('--job', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.stage in ('case', 'booking'):
        if args.job is None:
            parser.error('internal stage requires --job')
        (child_case if args.stage == 'case' else child_booking)(args.job)
    else:
        if args.outdir is None:
            parser.error('--outdir is required')
        if args.stage == 'plan':
            plan(args.outdir, args.golden.resolve(), args.reference.resolve())
        else:
            run(args.outdir)


if __name__ == '__main__':
    main()
