"""Direct-input chunk/merge safety, populated pp yield and scalable grid workflows."""
from argparse import Namespace
from copy import deepcopy
import builtins
from contextlib import redirect_stdout
import io
import json
import math
import os
import shutil
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from ruamel.yaml import YAML

PJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PJ / 'scripts'))
import inputs
import workflow as flow
from make_chunks import make
from cpp_fixtures import KEY, GRID_EXPECTED, GRID_SELECTED_EVENTS, grid_rows, write_part
from cpp_validation import root_snapshot, compare
from fixtures import ev, ph, tr, correction_file, MBD_NAMES

ROOT = flow.root_module('histmakers/PhotonJetHistMaker.C', 'scripts/CheckChunks.C', 'histmakers/InjectMbdEff.C')


def write_config(path, cfg):
    with Path(path).open('w') as stream:
        YAML().dump(cfg, stream)


def populated(product):
    """Populate every truth/reco bin with known ABCD signal and background."""
    rows = {name: [] for name in ('events', 'photons', 'truthPhotons', 'recoTruthLinks', 'truthJets')}
    for pt in (9., 11., 13., 15., 17., 19., 21., 23., 25., 27., 30., 34., 40.):
        if product == 'data' and pt < 16:
            continue
        for region, signal, background in (('A', 80, 20), ('B', 8, 40), ('C', 8, 40), ('D', 4, 80)):
            count = signal if product == 'sim_signal' else signal + background
            for item in range(count):
                index = len(rows['events'])
                matched = product != 'data' and item < signal
                event = dict(ev, event_id_lo=index + 1, source_entry=index,
                             event_sequence=index, physical_event_sequence=index, event_weight=2.)
                photon = dict(ph, **event, candidate_id_hi=11, candidate_id_lo=index + 1,
                              photon_et=pt, photon_encounter_ordinal=0)
                tight, isolated = region in ('A', 'B'), region in ('A', 'C')
                photon.update(bdt_score=.9 if tight else .6, bdt_is_tight=int(tight),
                    bdt_is_nontight=int(not tight), bdt_is_not_tight=int(not tight),
                    bdt_tight_threshold=.8, bdt_nontight_low_threshold=.2, bdt_nontight_high_threshold=.7,
                    iso_r03=0. if isolated else 5., iso_r04=0. if isolated else 5.,
                    iso_r03_threshold=1., iso_r04_threshold=1.,
                    iso_r03_nonisolated_threshold=2., iso_r04_nonisolated_threshold=2.,
                    iso_r03_pass=int(isolated), iso_r04_pass=int(isolated),
                    truth_matched=int(matched), truth_barcode=42 if matched else -1)
                rows['events'].append(event)
                rows['photons'].append(photon)
                if matched:
                    rows['truthPhotons'].append(dict(tr, **{k: event[k] for k in KEY}, truth_photon_pt=pt,
                        truth_photon_id_hi=12, truth_photon_id_lo=index + 1))
                    rows['recoTruthLinks'].append(dict({k: event[k] for k in KEY},
                        reco_type=1, truth_type=1, link_class=0, reco_index=0, truth_index=0,
                        reco_id_hi=11, reco_id_lo=index + 1, truth_id_hi=12, truth_id_lo=index + 1))
    return rows


def plan_jobs(config, chunk_dir, products):
    return make(Namespace(config=config, products=','.join(products),
        chunk_size=','.join(product + '=1' for product in products), outdir=chunk_dir,
        release_id='fixture-v1', indices=None))


def generation_paths(cfg, product):
    names = flow.output_names(cfg, product)
    return [*names.values(), Path(str(names['main']) + '.complete.yaml')]


class WorkflowTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = Path(os.environ.get('PJ_TEST_OUTDIR') or tempfile.mkdtemp(prefix='ppg12-2-cpp-workflow-'))
        cls.base.mkdir(parents=True, exist_ok=True)
        cls.output_alias = cls.base.parent / (cls.base.name + '-output-alias')
        cls.output_alias.symlink_to(cls.base, target_is_directory=True)
        cls.cfg = inputs.config(PJ / 'configs/pp/config_pj_pp_nom.yaml')
        analysis = cls.cfg['analysis']
        analysis.update(apply_trigger_eff_correction=0, mc_iso_scale=1., mc_iso_shift=0.,
                        reco_iso_max_b=1., reco_iso_max_s=0., reco_noniso_min_shift=1.)
        analysis['tight'].update(bdt_min_intercept=.8, bdt_min_slope=0.)
        analysis['non_tight'].update(bdt_min_intercept=.2, bdt_min_slope=0., bdt_max_intercept=.7, bdt_max_slope=0.)
        analysis['unfold']['reweight'] = 0
        for kind in ('data', 'mc'):
            for index in range(3):
                analysis[f'cluster_eres_{kind}_p{index}'] = 0.
        cls.cfg['photonjet'].update(apply_shower_shape_windows=0, flag_check='strict',
                                   match_source='links', release='fixture-v1')
        cls.cfg['output']['var_type'] = 'PPG12_2_cpp_workflow'
        for key, stem in (('data_outfile', 'data'), ('eff_outfile', 'eff'),
                          ('response_outfile', 'response'), ('final_outfile', 'final')):
            cls.cfg['output'][key] = str(cls.output_alias / stem)
        cls.source = cls.base / 'mbd.root'
        correction_file(cls.source, edges=analysis['pT_bins_truth'], corrected=True)
        cls.cfg['photonjet']['external_mbd_eff_file'] = str(cls.source)
        # Archive the actual correction selection, so compatibility is checked.
        with ROOT.TFile.Open(str(cls.source), 'UPDATE') as file:
            ROOT.TObjString(json.dumps(cls.cfg)).Write('config')
        cls.originals = {}
        for product, key in inputs.PRODUCTS.items():
            folder = cls.base / product
            folder.mkdir()
            rows = populated(product)
            cls.originals[product] = [folder / f'original_{index}.root' for index in range(2)]
            for path in cls.originals[product]:
                write_part(path, rows)
            master = folder / 'files.txt'
            master.write_text(''.join(str(path) + '\n' for path in cls.originals[product]))
            cls.cfg['photonjet'][key] = str(master)
        cls.config = cls.base / 'config.yaml'
        write_config(cls.config, cls.cfg)
        cls.args = Namespace(config=cls.config, chunk_dir=cls.base / 'chunks', products=None,
                             threads=1, log_dir=cls.base / 'logs', allow_subset=False)
        cls.plan = plan_jobs(cls.config, cls.args.chunk_dir, inputs.PRODUCTS)
        cls.plan_path = cls.args.chunk_dir / cls.cfg['output']['var_type'] / 'jobs.json'
        flow.local(cls.args)
        flow.merge(cls.args)
        cls.paths = [path for product in inputs.PRODUCTS for path in generation_paths(cls.cfg, product)]
        cls.baseline = {path: path.read_bytes() for path in cls.paths}

    def tearDown(self):
        self.restore_outputs()

    def restore_outputs(self):
        for path, content in self.baseline.items():
            path.write_bytes(content)

    def assert_outputs_unchanged(self):
        self.assertEqual({path: path.read_bytes() for path in self.paths}, self.baseline)

    def test_planning_is_stat_only_and_workers_validate_original_payloads(self):
        folder = self.base / 'stat-only-plan'
        folder.mkdir()
        original = folder / 'invalid.root'
        original.write_bytes(b'This has file metadata, but is not a ROOT file.\n')
        master = folder / 'files.txt'
        master.write_text(str(original) + '\n')
        cfg = deepcopy(self.cfg)
        cfg['photonjet']['data_files_list'] = str(master)
        cfg['output']['var_type'] = 'PPG12_2_stat_only'
        cfg['output']['data_outfile'] = str(folder / 'output')
        config = folder / 'config.yaml'
        write_config(config, cfg)
        # A fresh interpreter also catches eager imports hidden by the suite's
        # already-loaded ROOT module. Importing the planner needs neither reader.
        import_check = subprocess.run([sys.executable, '-c',
            "import sys; sys.modules.update(ROOT=None, uproot=None); sys.path.insert(0, sys.argv[1]); import make_chunks",
            str(PJ / 'scripts')], capture_output=True, text=True, timeout=20)
        self.assertEqual(import_check.returncode, 0, import_check.stdout + import_check.stderr)
        real_path_open, real_open, real_import = Path.open, builtins.open, builtins.__import__

        def guarded_path_open(path, *args, **kwargs):
            self.assertNotEqual(Path(path).resolve(), original.resolve(), 'planner opened an original payload')
            return real_path_open(path, *args, **kwargs)

        def guarded_open(path, *args, **kwargs):
            if isinstance(path, (str, bytes, os.PathLike)):
                self.assertNotEqual(Path(os.fsdecode(path)).resolve(), original.resolve(), 'planner opened an original payload')
            return real_open(path, *args, **kwargs)

        def guarded_import(name, *args, **kwargs):
            self.assertNotIn(name.split('.')[0], ('ROOT', 'uproot'), 'planning imported a ROOT reader')
            return real_import(name, *args, **kwargs)

        with patch.object(Path, 'open', guarded_path_open), patch.object(builtins, 'open', guarded_open), \
             patch.object(builtins, '__import__', guarded_import):
            plan = plan_jobs(config, folder / 'chunks', ['data'])
        part = plan['products']['data']['jobs'][0]['parts'][0]
        self.assertEqual(part, dict(original_part_index=0,
            original_part_id='fixture-v1/data/pp/part/0', source_path=str(original),
            source_bytes=original.stat().st_size, source_mtime_ns=original.stat().st_mtime_ns))
        self.assertEqual(plan['format_version'], 4)
        self.assertEqual(plan['input_verification'], 'worker-md5-stat-v1')
        job = plan['products']['data']['jobs'][0]
        with self.assertRaises(Exception):
            flow.worker(Namespace(config=config, product='data', input=job['input'], tag=job['tag'], threads=1))
        self.assertFalse(flow.output_names(cfg, 'data', job['tag'])['main'].exists())
        self.assert_outputs_unchanged()

    @unittest.skipUnless(shutil.which('strace'), 'strace is required to observe native original-file access')
    def test_normal_status_and_merge_do_not_open_originals(self):
        trace = self.base / 'status-merge-open.trace'
        program = ("import sys; from pathlib import Path; from argparse import Namespace; "
            "sys.path.insert(0, sys.argv[1]); import workflow; "
            "args=Namespace(config=Path(sys.argv[2]),chunk_dir=Path(sys.argv[3]),"
            "log_dir=Path(sys.argv[4]),products='data',threads=1,allow_subset=False); "
            "workflow.status(args); workflow.merge(args)")
        result = subprocess.run(['strace', '-f', '-qq', '-s', '4096', '-e', 'trace=open,openat', '-o', str(trace),
            sys.executable, '-c', program, str(PJ / 'scripts'), str(self.config),
            str(self.args.chunk_dir), str(self.args.log_dir)], cwd=PJ,
            capture_output=True, text=True, timeout=180)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for product in inputs.PRODUCTS:
            self.assertIn(f'{product}: 2/2 checked complete generations', result.stdout)
        opened = trace.read_text()
        self.assertIn(json.dumps(str(self.config.resolve())), opened, 'trace did not capture expected config opens')
        data_job = self.plan['products']['data']['jobs'][0]
        self.assertIn(json.dumps(str(flow.output_names(self.cfg, 'data', data_job['tag'])['main'])), opened,
                      'trace did not capture the merged chunk being read')
        for paths in self.originals.values():
            for path in paths:
                self.assertNotIn(json.dumps(str(path)), opened, 'normal verification opened an original payload')
                self.assertNotIn(json.dumps(str(path.resolve())), opened, 'normal verification opened a canonical original payload')
        print('[workflow] native status + merge opened no originals; trace:', trace, flush=True)

    def test_worker_post_read_hash_rejects_same_stat_payload_mutation(self):
        # Invoke the actual worker pre/post verification boundaries around a
        # deterministic mutation, without timing races or production hooks.
        self.assertTrue(ROOT.gInterpreter.Declare(r'''
            namespace PJVerificationTest {
            struct Prepared {
                PJ::Inputs inputs;
                YAML::Node provenance;
                Prepared(const std::string &path, const std::string &list) {
                    const auto config = PJ::ReadConfig(path);
                    const auto cuts = PJ::LoadCuts(config.yaml);
                    inputs = PJ::ReadInputs(list, config, PJ::Product::Data, cuts);
                    const auto model = PJ::Configure(config, PJ::Product::Data, inputs);
                    provenance = PJ::Provenance(config, PJ::Product::Data, inputs, model);
                }
                void Recheck() const { PJ::RecheckInputs(inputs, provenance); }
            };
            }
        '''))
        job = self.plan['products']['data']['jobs'][0]
        prepared = ROOT.PJVerificationTest.Prepared(str(self.config), job['input'])
        original = self.originals['data'][0]
        saved, saved_stat = original.read_bytes(), original.stat()
        try:
            modified = bytearray(saved)
            modified[-1] ^= 1
            original.write_bytes(modified)
            os.utime(original, ns=(saved_stat.st_atime_ns, saved_stat.st_mtime_ns))
            self.assertEqual(original.stat().st_size, saved_stat.st_size)
            self.assertEqual(original.stat().st_mtime_ns, saved_stat.st_mtime_ns)
            with self.assertRaisesRegex(Exception, 'fingerprint|hash|digest|changed'):
                prepared.Recheck()
        finally:
            original.write_bytes(saved)
            os.utime(original, ns=(saved_stat.st_atime_ns, saved_stat.st_mtime_ns))
        prepared.Recheck()
        self.assert_outputs_unchanged()

    def test_duplicate_uuid_across_jobs_is_rejected_at_merge(self):
        folder = self.base / 'duplicate-uuid'
        folder.mkdir()
        first, second = folder / 'first.root', folder / 'second.root'
        shutil.copy2(self.originals['data'][0], first)
        shutil.copy2(first, second)
        master = folder / 'files.txt'
        master.write_text(f'{first}\n{second}\n')
        cfg = deepcopy(self.cfg)
        cfg['photonjet']['data_files_list'] = str(master)
        cfg['output']['var_type'] = 'PPG12_2_duplicate_uuid'
        cfg['output']['data_outfile'] = str(folder / 'output')
        config = folder / 'config.yaml'
        write_config(config, cfg)
        args = Namespace(config=config, chunk_dir=folder / 'chunks', products='data',
                         threads=1, log_dir=folder / 'logs', allow_subset=False)
        plan = plan_jobs(config, args.chunk_dir, ['data'])
        self.assertEqual(len(plan['products']['data']['jobs']), 2)
        flow.local(args)
        records = [flow.checked_generation(ROOT, flow.output_names(cfg, 'data', job['tag'])['main'], config, 'data')
                   for job in plan['products']['data']['jobs']]
        self.assertEqual(records[0]['inputs'][0]['source_uuid'], records[1]['inputs'][0]['source_uuid'])
        with self.assertRaisesRegex(Exception, 'duplicate|coverage|UUID'):
            flow.merge(args)
        self.assertFalse(flow.output_names(cfg, 'data')['main'].exists())
        self.assert_outputs_unchanged()

    def test_source_identity_uses_directory_uuid_not_file_header_uuid(self):
        import uuid
        import uproot
        path = self.base / 'distinct-root-uuids.root'
        write_part(path, {name: [] for name in ('events', 'photons', 'truthPhotons', 'recoTruthLinks', 'truthJets')})
        file = ROOT.TFile.Open(str(path))
        expected = str(file.GetUUID().AsString())
        file.Close()
        with uproot.open(path) as source:
            header_uuid = source.file.fUUID
        content = bytearray(path.read_bytes())
        offset = content[:100].find(header_uuid)
        self.assertGreaterEqual(offset, 0)
        replacement = uuid.UUID('11111111-2222-3333-8444-555555555555')
        content[offset:offset+16] = replacement.bytes
        path.write_bytes(content)
        with uproot.open(path) as source:
            self.assertEqual(str(source.file.uuid), str(replacement))
        self.assertNotEqual(expected, str(replacement))
        record = inputs.inspect_source(path, 0, self.cfg, 'data', 'fixture-v1')
        self.assertEqual(record['source_uuid'], expected)
        file = ROOT.TFile.Open(str(path))
        self.assertEqual(record['source_uuid'], str(file.GetUUID().AsString()))
        file.Close()
        self.assert_outputs_unchanged()

    def test_aliased_destinations_are_checked_before_publication(self):
        cfg = deepcopy(self.cfg)
        cfg['output']['response_outfile'] = str(self.base / 'eff')
        path = self.base / 'same-destination.yaml'
        write_config(path, cfg)
        job = self.plan['products']['sim_signal']['jobs'][0]
        with self.assertRaisesRegex(Exception, 'destinations must differ|paths must differ'):
            ROOT.PJ.Run(str(path), 'sim_signal', '', '')
        self.assert_outputs_unchanged()
        cfg = deepcopy(self.cfg)
        cfg['photonjet']['external_mbd_eff_file'] = str(flow.output_names(cfg, 'sim_signal')['main'])
        with patch.object(flow, 'load_plan', return_value=(cfg, self.plan_path, self.plan)):
            with self.assertRaisesRegex(ValueError, 'source and final target'):
                flow.merge(self.args)
        self.assert_outputs_unchanged()

    def test_local_wrapper_propagates_worker_failure(self):
        package = self.base / 'wrapper'
        (package / 'scripts').mkdir(parents=True, exist_ok=True)
        shutil.copy2(PJ / 'oneforall_photonjet.sh', package / 'oneforall_photonjet.sh')
        (package / 'scripts/workflow.py').write_text("import sys; print('simulated worker failure'); sys.exit(17)\n")
        result = subprocess.run(['bash', str(package / 'oneforall_photonjet.sh'), str(self.config), 'local'],
                                capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 17, result.stdout + result.stderr)

    def test_complete_merge_matches_single_run(self):
        self.assertEqual(self.plan['format_version'], 4)
        self.assertEqual(self.plan['input_verification'], 'worker-md5-stat-v1')
        self.assertEqual(self.plan['backend'], inputs.BACKEND)
        for product in inputs.PRODUCTS:
            ROOT.PJ.Run(str(self.config), product, '', 'together')
            for kind, path in flow.output_names(self.cfg, product).items():
                extra = ('mbd_eff_source', 'mbd_eff_provenance') if product == 'sim_signal' and kind == 'main' else ()
                merged = root_snapshot(path, extra)['objects']
                together = root_snapshot(flow.output_names(self.cfg, product, 'together')[kind])['objects']
                if extra:
                    for name in MBD_NAMES:
                        merged.pop(name)
                        together.pop(name)
                self.assertEqual(compare(together, merged), [])
            record = flow.checked_generation(ROOT, flow.output_names(self.cfg, product)['main'], self.config, product)
            self.assertEqual(record['format_version'], 4)
            self.assertEqual(record['input_verification'], 'worker-md5-stat-v1')
            self.assertEqual(record['job_plan'], dict(path=str(self.plan_path.resolve()), md5=inputs.md5(self.plan_path),
                assignments=[dict(input=str(Path(job['input']).resolve()), tag=job['tag'])
                             for job in self.plan['products'][product]['jobs']]))
            marker = json.loads(Path(str(flow.output_names(self.cfg, product)['main']) + '.complete.yaml').read_text())
            self.assertEqual(marker['format_version'], 3)
            self.assertEqual(marker['input_verification'], 'worker-md5-stat-v1')
            for part in record['inputs']:
                source = Path(part['source_path'])
                self.assertEqual(part['source_md5'], inputs.md5(source))
                self.assertEqual(part['source_bytes'], source.stat().st_size)
                self.assertEqual(part['source_mtime_ns'], source.stat().st_mtime_ns)
                self.assertGreater(part['events'], 0)
                with ROOT.TFile.Open(str(source)) as original:
                    self.assertEqual(part['source_uuid'], str(original.GetUUID().AsString()))
                planned = next(entry for job in self.plan['products'][product]['jobs'] for entry in job['parts']
                               if entry['original_part_index'] == part['original_part_index'])
                self.assertEqual({key: part[key] for key in planned}, planned)
            self.assertEqual(record['workers'], 1)
            self.assertEqual(record['parts_processed'], 2)
            self.assertEqual(record['events_processed'], sum(part['events'] for part in record['inputs']))
        with ROOT.TFile.Open(str(flow.output_names(self.cfg, 'sim_signal')['main'])) as file:
            self.assertEqual(file.Get(MBD_NAMES[1]).GetBinContent(1) / file.Get(MBD_NAMES[0]).GetBinContent(1), .4)
            self.assertTrue(YAML(typ='safe').load(str(file.Get('mbd_eff_provenance').GetString()))['source_selection_checked'])

    def test_completion_config_and_version_rejections(self):
        main = flow.output_names(self.cfg, 'sim_signal')['main']
        marker = Path(str(main) + '.complete.yaml')
        original = marker.read_bytes()
        marker.unlink()
        with self.assertRaises(Exception):
            flow.checked_generation(ROOT, main, self.config, 'sim_signal')
        marker.write_bytes(original)
        with main.open('ab') as stream:
            stream.write(b'changed')
        with self.assertRaises(Exception):
            flow.checked_generation(ROOT, main, self.config, 'sim_signal')
        main.write_bytes(self.baseline[main])
        alternate = self.base / 'changed.yaml'
        alternate.write_text(self.config.read_text() + '\n# stale config\n')
        with self.assertRaises(Exception):
            flow.checked_generation(ROOT, main, alternate, 'sim_signal')
        for field in ('rng_version', 'schema_version', 'backend', 'code', 'dependencies', 'format_version', 'input_verification', 'job_plan'):
            seal = json.loads(original)
            for entry in seal['files']:
                with ROOT.TFile.Open(entry['path'], 'UPDATE') as file:
                    info = YAML(typ='safe').load(str(file.Get('photonjet_provenance').GetString()))
                    if field == 'code':
                        info[field][next(iter(info[field]))] = '0' * 32
                    elif field == 'dependencies':
                        dependency = self.base / 'changed-dependency.txt'
                        dependency.write_text('before')
                        info[field][str(dependency)] = inputs.md5(dependency)
                        dependency.write_text('after')
                    elif field == 'job_plan':
                        info[field]['assignments'][0]['tag'] = 'unassigned'
                    else:
                        info[field] = 3 if field == 'format_version' else 'wrong'
                    file.Delete('photonjet_provenance;*')
                    ROOT.TObjString(json.dumps(info)).Write('photonjet_provenance')
                entry.update(md5=inputs.md5(entry['path']), bytes=Path(entry['path']).stat().st_size)
            inputs.save(marker, seal)
            with self.assertRaises(Exception, msg=field):
                flow.checked_generation(ROOT, main, self.config, 'sim_signal')
            self.restore_outputs()
        for field, value in (('format_version', 2), ('backend', 'wrong'), ('input_verification', 'wrong')):
            seal = json.loads(original)
            seal[field] = value
            inputs.save(marker, seal)
            with self.assertRaises(Exception, msg=field):
                flow.checked_generation(ROOT, main, self.config, 'sim_signal')
            self.restore_outputs()

    def test_duplicate_missing_and_extra_original_coverage(self):
        jobs = self.plan['products']['sim_signal']['jobs']
        mains = [flow.output_names(self.cfg, 'sim_signal', job['tag'])['main'] for job in jobs]
        listing = self.base / 'bad-coverage.list'
        for paths in ([mains[0], mains[0]], [mains[0]]):
            listing.write_text(''.join(str(path) + '\n' for path in paths))
            with self.assertRaises(Exception):
                ROOT.PJ.ValidateChunks(str(listing), str(self.config), 'sim_signal', str(self.plan_path))
        hidden = mains[0].with_suffix('.missing')
        mains[0].rename(hidden)
        try:
            with self.assertRaisesRegex(ValueError, 'missing|unplanned'):
                flow.merge(self.args)
        finally:
            hidden.rename(mains[0])
        extra = mains[0].parent / (flow.output_names(self.cfg, 'sim_signal')['main'].stem + '.sim_signal_unplanned.root')
        extra.write_bytes(mains[0].read_bytes())
        try:
            with self.assertRaisesRegex(ValueError, 'missing|unplanned'):
                flow.merge(self.args)
        finally:
            extra.unlink()
        # Exact indices and the expected chunk filename do not establish that
        # an otherwise valid direct run belongs to this sealed job plan.
        names = flow.output_names(self.cfg, 'sim_signal', jobs[0]['tag'])
        marker = Path(str(names['main']) + '.complete.yaml')
        saved = {path: path.read_bytes() for path in [*names.values(), marker]}
        unsealed = self.base / 'unbound-original.list'
        unsealed.write_text(f"0 {self.originals['sim_signal'][0]}\n")
        try:
            ROOT.PJ.Run(str(self.config), 'sim_signal', str(unsealed), jobs[0]['tag'])
            record = flow.checked_generation(ROOT, names['main'], self.config, 'sim_signal')
            self.assertEqual([part['original_part_index'] for part in record['inputs']], [0])
            self.assertFalse(record.get('job_plan'))
            with self.assertRaisesRegex(Exception, 'plan|assignment|bound'):
                flow.merge(self.args)
        finally:
            for path, content in saved.items():
                path.write_bytes(content)
        self.assert_outputs_unchanged()

    def test_stale_plans_lists_and_original_dependencies(self):
        saved_plan = self.plan_path.read_bytes()
        job = self.plan['products']['sim_signal']['jobs'][0]
        for field, value in (('format_version', 3), ('backend', 'wrong'), ('schema_version', 'wrong'),
                             ('input_verification', 'wrong'),
                             ('config_md5', '0' * 32), ('code', {})):
            changed = deepcopy(self.plan)
            changed[field] = value
            inputs.save(self.plan_path, changed)
            try:
                with self.assertRaises(Exception, msg=field):
                    flow.load_plan(self.config, self.args.chunk_dir)
                with self.assertRaises(Exception, msg=field):
                    ROOT.PJ.Run(str(self.config), 'sim_signal', job['input'], job['tag'])
                with self.assertRaises(Exception, msg=field):
                    flow.checked_generation(ROOT, flow.output_names(self.cfg, 'sim_signal')['main'], self.config, 'sim_signal')
            finally:
                self.plan_path.write_bytes(saved_plan)
        job = self.plan['products']['sim_signal']['jobs'][0]
        worker = Namespace(config=self.config, product='sim_signal', input=job['input'], tag=job['tag'], threads=1)
        paths = [Path(job['input']), Path(self.plan['products']['sim_signal']['master']),
                 self.originals['sim_signal'][0], self.source]
        for path in paths:
            original = path.read_bytes()
            original_stat = path.stat()
            try:
                with path.open('ab') as stream:
                    stream.write(b'\n# changed input\n')
                with self.assertRaises(Exception, msg=str(path)):
                    flow.worker(worker)
                with self.assertRaises(Exception, msg=str(path)):
                    flow.checked_generation(ROOT, flow.output_names(self.cfg, 'sim_signal')['main'], self.config, 'sim_signal')
            finally:
                path.write_bytes(original)
                os.utime(path, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
        # A same-size source with a new mtime is stale even when its bytes match.
        source = self.originals['sim_signal'][0]
        saved_stat = source.stat()
        try:
            os.utime(source, ns=(saved_stat.st_atime_ns, saved_stat.st_mtime_ns + 1000000000))
            with self.assertRaises(Exception):
                flow.worker(worker)
            with self.assertRaises(Exception):
                flow.checked_generation(ROOT, flow.output_names(self.cfg, 'sim_signal')['main'], self.config, 'sim_signal')
            with self.assertRaises(Exception):
                flow.merge(self.args)
            output = io.StringIO()
            with redirect_stdout(output):
                flow.status(self.args)
            self.assertIn('sim_signal: 1/2 checked complete generations', output.getvalue())
        finally:
            os.utime(source, ns=(saved_stat.st_atime_ns, saved_stat.st_mtime_ns))
        self.assert_outputs_unchanged()

    def test_mbd_failures_and_empty_source_preserve_output(self):
        target = flow.output_names(self.cfg, 'sim_signal')['main']
        source_bytes = self.source.read_bytes()
        source_stat = self.source.stat()
        try:
            for mode in ('missing', 'axis', 'numerator', 'selection'):
                self.source.write_bytes(source_bytes)
                with ROOT.TFile.Open(str(self.source), 'UPDATE') as file:
                    if mode == 'missing':
                        file.Delete(MBD_NAMES[-1] + ';*')
                    elif mode == 'axis':
                        file.Delete(MBD_NAMES[0] + ';*')
                        ROOT.TH1D(MBD_NAMES[0], 'bad axis', 1, 8, 45).Write()
                    elif mode == 'numerator':
                        histogram = file.Get(MBD_NAMES[1])
                        histogram.SetBinContent(1, 101)
                        histogram.Write(MBD_NAMES[1], ROOT.TObject.kOverwrite)
                    else:
                        cfg = deepcopy(self.cfg)
                        cfg['analysis']['vertex_cut'] = 20.
                        ROOT.TObjString(json.dumps(cfg)).Write('config', ROOT.TObject.kOverwrite)
                with self.assertRaises(Exception, msg=mode):
                    ROOT.PJ.InjectCorrection(str(self.config), str(target))
                self.assertEqual(target.read_bytes(), self.baseline[target])
            cfg = deepcopy(self.cfg)
            cfg['photonjet']['external_mbd_eff_file'] = ''
            path = self.base / 'empty-mbd.yaml'
            write_config(path, cfg)
            ROOT.PJ.InjectCorrection(str(path), str(target))
            self.assertEqual(target.read_bytes(), self.baseline[target])
        finally:
            self.source.write_bytes(source_bytes)
            os.utime(self.source, ns=(source_stat.st_atime_ns, source_stat.st_mtime_ns))

    def test_processing_and_mid_publication_failures(self):
        original = self.originals['sim_signal'][0]
        original_bytes = original.read_bytes()
        original_stat = original.stat()
        try:
            rows = populated('sim_signal')
            rows['photons'][0]['vertex_z'] += 1.
            write_part(original, rows)
            with self.assertRaisesRegex(Exception, 'photon/event disagreement'):
                ROOT.PJ.Run(str(self.config), 'sim_signal', '', '')
            self.assert_outputs_unchanged()
        finally:
            original.write_bytes(original_bytes)
            os.utime(original, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
        with patch.object(flow, 'merge_root', side_effect=RuntimeError('injected merge failure')):
            with self.assertRaisesRegex(RuntimeError, 'injected'):
                flow.merge(self.args)
        self.assert_outputs_unchanged()
        replace, count = os.replace, 0
        def fail_second_root(source, target):
            nonlocal count
            if str(target).endswith('.root') and Path(target) in self.paths:
                count += 1
                if count == 2:
                    raise OSError('injected second publication failure')
            return replace(source, target)
        with patch.object(flow.os, 'replace', side_effect=fail_second_root):
            with self.assertRaisesRegex(OSError, 'second publication'):
                flow.merge(self.args)
        self.assertEqual(count, 2)
        for product in inputs.PRODUCTS:
            main = flow.output_names(self.cfg, product)['main']
            self.assertFalse(Path(str(main) + '.complete.yaml').exists())
            with self.assertRaises(Exception):
                flow.checked_generation(ROOT, main, self.config, product)
        with self.assertRaises(Exception):
            flow.yield_outputs(self.args)

    def test_worker_cpu_allocation(self):
        job = self.plan['products']['data']['jobs'][0]
        args = Namespace(config=self.config, product='data', input=job['input'], tag=job['tag'], threads=2)
        with self.assertRaisesRegex(ValueError, 'one|1|serial'):
            flow.worker(args)
        ad = self.base / 'job.ad'
        ad.write_text('RequestCpus = 4\n')
        args.threads = 1
        with patch.dict(os.environ, {'_CONDOR_JOB_AD': str(ad)}):
            with self.assertRaisesRegex(ValueError, 'RequestCpus'):
                flow.worker(args)
        self.assert_outputs_unchanged()

    def test_yield_and_mc_closure(self):
        flow.yield_outputs(self.args)
        marker = Path(self.cfg['output']['final_outfile'] + '_' + self.cfg['output']['var_type'] + '.root.complete.yaml')
        result = json.loads(marker.read_text())
        self.assertEqual({entry['kind'] for entry in result['files']}, {'data', 'mc'})
        for entry in result['files']:
            self.assertEqual(entry['md5'], inputs.md5(entry['path']))
            with ROOT.TFile.Open(entry['path']) as file:
                histogram = file.Get('h_unfold_sub_result')
                self.assertGreater(histogram.Integral(), 0.)
                for index in range(histogram.GetNcells()):
                    self.assertTrue(math.isfinite(histogram.GetBinContent(index)))
                    self.assertTrue(math.isfinite(histogram.GetBinError(index)))
                self.assertEqual(str(file.Get('config').GetString()), self.config.read_text())
        final_paths = [Path(entry['path']) for entry in result['files']] + [marker]
        saved = {path: path.read_bytes() for path in final_paths}
        with patch.object(flow.subprocess, 'run', side_effect=RuntimeError('injected yield processing failure')):
            with self.assertRaisesRegex(RuntimeError, 'injected yield'):
                flow.yield_outputs(self.args)
        self.assertEqual({path: path.read_bytes() for path in final_paths}, saved)
        print('[workflow] complete populated pp yield and MC closure', result['files'], flush=True)

    def test_multibin_chunk_merge_and_consumer_rejections(self):
        folder = self.base / 'grid'
        folder.mkdir()
        cfg = deepcopy(self.cfg)
        cfg['analysis']['eta_bins'] = [-.7, 0., .7]
        cfg['photonjet'].update(system='auau', centrality_bins=[0., 20., 50., 80.],
            weight_mode='sample_map', threshold_source='tree', external_mbd_eff_file='')
        sample_map = folder / 'samples.yaml'
        write_config(sample_map, {'ranges': [{'first': 0, 'last': 1, 'sample': 'photon10'}]})
        cfg['photonjet']['sample_map_signal'] = str(sample_map)
        sample_weight = float(ROOT.PPG12.GetSampleConfig('photon10').weight)
        cfg['output']['var_type'] = 'PPG12_2_cpp_workflow_grid'
        for key in ('data_outfile', 'eff_outfile', 'response_outfile', 'final_outfile'):
            cfg['output'][key] = str(folder / key)
        paths = [folder / f'original_{index}.root' for index in range(2)]
        rows = grid_rows()
        for path in paths:
            write_part(path, rows, auau=True)
        master = folder / 'files.txt'
        master.write_text(''.join(str(path) + '\n' for path in paths))
        cfg['photonjet']['sim_signal_files_list'] = str(master)
        config = folder / 'config.yaml'
        write_config(config, cfg)
        args = Namespace(config=config, chunk_dir=folder / 'chunks', products='sim_signal',
                         threads=1, log_dir=folder / 'logs', allow_subset=False)
        plan = plan_jobs(config, args.chunk_dir, ['sim_signal'])
        flow.local(args)
        flow.merge(args)
        ROOT.PJ.Run(str(config), 'sim_signal', '', 'together')
        names = flow.output_names(cfg, 'sim_signal')
        for kind, path in names.items():
            merged = root_snapshot(path)['objects']
            together = root_snapshot(flow.output_names(cfg, 'sim_signal', 'together')[kind])['objects']
            self.assertEqual(compare(together, merged), [])
            self.assertEqual(len(merged), 206 if kind == 'main' else 25)
        record = flow.checked_generation(ROOT, names['main'], config, 'sim_signal')
        self.assertEqual(record['layout'], inputs.bin_layout(cfg))
        self.assertEqual(record['events_processed'], 2 * len(rows['events']))
        self.assertEqual(record['parts_processed'], 2)
        self.assertEqual(record['selected_events_by_centrality'], [2 * count for count in GRID_SELECTED_EVENTS])
        with ROOT.TFile.Open(str(names['main'])) as file:
            self.assertEqual(file.Get('h_pj_cutflow').GetBinContent(1), 2 * len(rows['events']))
            self.assertEqual(file.Get('h_pj_vertex_z').GetEntries(), 2 * sum(GRID_SELECTED_EVENTS))
            for (centrality, eta), expected in GRID_EXPECTED.items():
                suffix = f'_cent{centrality}_eta{eta}'
                self.assertEqual(file.Get('h_all_cluster' + suffix).Integral(), 2 * sample_weight * expected['reco'])
                self.assertEqual(file.Get('h_truth_pT' + suffix).Integral(), 2 * sample_weight * expected['truth'])
        saved = {path: path.read_bytes() for path in generation_paths(cfg, 'sim_signal')}
        chunk = flow.output_names(cfg, 'sim_signal', plan['products']['sim_signal']['jobs'][0]['tag'])
        chunk_marker = Path(str(chunk['main']) + '.complete.yaml')
        chunk_saved = {path: path.read_bytes() for path in [*chunk.values(), chunk_marker]}
        try:
            seal = YAML(typ='safe').load(chunk_marker.read_text())
            for entry in seal['files']:
                with ROOT.TFile.Open(entry['path'], 'UPDATE') as file:
                    info = YAML(typ='safe').load(str(file.Get('photonjet_provenance').GetString()))
                    info['layout']['eta_edges'] = [-.7, .1, .7]
                    file.Delete('photonjet_provenance;*')
                    ROOT.TObjString(json.dumps(info)).Write('photonjet_provenance')
                entry.update(md5=inputs.md5(entry['path']), bytes=Path(entry['path']).stat().st_size)
            inputs.save(chunk_marker, seal)
            with self.assertRaisesRegex(Exception, 'layout'):
                flow.merge(args)
            self.assertEqual({path: path.read_bytes() for path in saved}, saved)
        finally:
            for path, content in chunk_saved.items():
                path.write_bytes(content)
        mbd_cfg = deepcopy(cfg)
        mbd_cfg['photonjet']['external_mbd_eff_file'] = str(self.source)
        mbd_config = folder / 'unsupported-mbd.yaml'
        write_config(mbd_config, mbd_cfg)
        with self.assertRaisesRegex(Exception, 'multiple|multi-cell|multi-bin'):
            ROOT.PJ.InjectCorrection(str(mbd_config), str(names['main']))
        self.assertEqual({path: path.read_bytes() for path in saved}, saved)
        # A pp two-eta plan must reject before looking for yield inputs. This
        # isolates layout support from the separately unsupported AuAu yield.
        pp_cfg = deepcopy(cfg)
        pp_cfg['photonjet'].update(system='pp', centrality_bins=[])
        pp_cfg['output']['var_type'] = 'PPG12_2_cpp_workflow_two_eta'
        pp_config = folder / 'two-eta.yaml'
        write_config(pp_config, pp_cfg)
        plan_jobs(pp_config, args.chunk_dir, ['sim_signal'])
        pp_args = Namespace(**vars(args))
        pp_args.config = pp_config
        with self.assertRaisesRegex(ValueError, 'multiple|multi-cell|multi-bin'):
            flow.yield_outputs(pp_args)
        self.assertFalse(Path(pp_cfg['output']['final_outfile'] + '_' + pp_cfg['output']['var_type'] + '.root').exists())
        self.assert_outputs_unchanged()


if __name__ == '__main__':
    unittest.main(verbosity=2)
