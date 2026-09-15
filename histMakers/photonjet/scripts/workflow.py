#!/usr/bin/env python3
"""Checked serial original-input workers, staged merge/MBD and the pp yield chain."""
import argparse
from contextlib import ExitStack
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import tempfile
from inputs import (PJ, PRODUCTS, BACKEND, SCHEMA, config, label, md5, save,
                    bin_layout, code_digests, require_single_cell, PLAN_FORMAT_VERSION,
                    PROVENANCE_FORMAT_VERSION, COMPLETION_FORMAT_VERSION, INPUT_VERIFICATION)


_loaded_sources=set()


def root_module(*sources):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load('libyaml-cpp', '', True)
    ROOT.gSystem.Load('libRooUnfold', '', True)
    for source in sources:
        if source in _loaded_sources: continue
        if not ROOT.gInterpreter.Declare('#include ' + json.dumps(str(PJ/source))):
            raise RuntimeError('cannot load ' + source)
        _loaded_sources.add(source)
    return ROOT


def output_path(name):
    """Resolve directory aliases while preserving the published filename."""
    path=Path(name)
    return path.parent.resolve()/path.name


def output_names(cfg, product, tag=''):
    out=cfg['output']; variation=label(str(out['var_type']))
    stem=out['data_outfile'] if product=='data' else out['eff_outfile']+('_jet' if product=='sim_inclusive' else '')
    names={'main':output_path(stem+'_'+variation+'.root')}
    if product=='sim_signal': names['response']=output_path(out['response_outfile']+'_'+variation+'.root')
    if len(set(names.values()))!=len(names): raise ValueError('main and response output paths must differ')
    if tag:
        label(tag)
        names={kind:output_path(path.parent/'chunks'/f'{path.stem}.{tag}.root') for kind,path in names.items()}
    return names


def load_plan(config_path, chunk_dir):
    cfg=config(config_path)
    path=Path(chunk_dir).resolve()/label(str(cfg['output']['var_type']))/'jobs.json'
    plan=json.loads(path.read_text())
    if (plan.get('format_version')!=PLAN_FORMAT_VERSION or plan.get('backend')!=BACKEND or
            plan.get('schema_version')!=SCHEMA or plan.get('input_verification')!=INPUT_VERIFICATION):
        raise ValueError('incompatible job plan or input verification policy; regenerate chunks')
    if Path(plan['config']).resolve()!=Path(config_path).resolve() or plan['config_md5']!=md5(config_path) or plan['layout']!=bin_layout(cfg):
        raise ValueError('stale config/layout in job plan')
    if plan['code']!=code_digests(): raise ValueError('maker code changed since chunk planning')
    for source,digest in plan['list_md5'].items():
        if md5(source)!=digest: raise ValueError('job list changed: '+source)
    for product,record in plan['products'].items():
        if product not in PRODUCTS or md5(record['master'])!=record['master_md5']:
            raise ValueError('master list changed since job planning')
        for source,digest in record['dependencies'].items():
            if md5(source)!=digest: raise ValueError('dependency changed since job planning: '+source)
    # Check all final/chunk destinations together, including directory aliases.
    destinations=[]
    for product,record in plan['products'].items():
        destinations.extend(output_names(cfg,product).values())
        for job in record['jobs']: destinations.extend(output_names(cfg,product,job['tag']).values())
    if len(destinations)!=len(set(destinations)): raise ValueError('output destinations collide across products or jobs')
    protected={Path(config_path).resolve(),path.resolve(),*(Path(name).resolve() for name in plan['list_md5'])}
    for record in plan['products'].values():
        protected.add(Path(record['master']).resolve())
        protected.update(Path(name).resolve() for name in record['dependencies'])
        protected.update(Path(part['source_path']).resolve() for job in record['jobs'] for part in job['parts'])
    if any(destination.resolve() in protected for destination in destinations):
        raise ValueError('output destination collides with an input or sealed plan')
    return cfg,path,plan


def checked_generation(ROOT, main, config_path, product):
    from ruamel.yaml import YAML
    info=YAML(typ='safe').load(str(ROOT.PJ.CheckCompletionText(str(main),str(config_path),product)))
    if (info.get('format_version')!=PROVENANCE_FORMAT_VERSION or info.get('backend')!=BACKEND or
            info.get('schema_version')!=SCHEMA or info.get('input_verification')!=INPUT_VERIFICATION):
        raise ValueError('incompatible histogram provenance or input verification policy')
    return info


def job_plan_binding(plan_path, jobs):
    path=Path(plan_path).resolve()
    return dict(path=str(path),md5=md5(path),
                assignments=[dict(input=str(Path(job['input']).resolve()),tag=job['tag']) for job in jobs])


def require_job_plan(info, plan_path, jobs):
    if info.get('job_plan')!=job_plan_binding(plan_path,jobs):
        raise ValueError('histogram output is not bound to the current job plan and assignments')


def require_serial(threads):
    if threads != 1:
        raise ValueError('the direct maker runs serially; request one worker/CPU and use Condor jobs for parallelism')


def worker(args):
    require_serial(args.threads)
    ad=os.environ.get('_CONDOR_JOB_AD')
    if ad:
        allocation=None
        for line in Path(ad).read_text().splitlines():
            key,separator,value=line.partition('=')
            if separator and key.strip().lower()=='requestcpus': allocation=int(value.strip())
        if allocation!=1: raise ValueError('the direct worker requires Condor RequestCpus = 1')
    cfg,plan_path,plan=load_plan(args.config,Path(args.input).resolve().parent.parent)
    matches=[j for j in plan['products'][args.product]['jobs'] if j['tag']==args.tag and Path(j['input']).resolve()==Path(args.input).resolve()]
    if len(matches)!=1: raise ValueError('worker input/tag is not in the sealed job plan')
    ROOT=root_module('histmakers/PhotonJetHistMaker.C')
    ROOT.PJ.Run(str(args.config),args.product,str(Path(args.input).resolve()),args.tag)
    info=checked_generation(ROOT,output_names(cfg,args.product,args.tag)['main'],args.config,args.product)
    require_job_plan(info,plan_path,matches)


def merge_root(ROOT, paths, output):
    merger=ROOT.TFileMerger()
    if not merger.OutputFile(str(output),'RECREATE'): raise RuntimeError('cannot open staged merge output')
    for path in paths:
        if not merger.AddFile(str(path)): raise RuntimeError('cannot add merge input: '+str(path))
    if not merger.Merge(): raise RuntimeError('ROOT merge failed: '+str(output))
    # Destroy the merger before reopening its output for finalization.
    del merger


def write_metadata(ROOT, path, original_config, info):
    file=ROOT.TFile.Open(str(path),'UPDATE')
    if not file or file.IsZombie(): raise RuntimeError('cannot update merged metadata')
    try:
        # ROOT combines TEfficiency's optional global sample weight harmonically.
        # These shards partition one sample: merge the counts and retain the
        # neutral global weight of a single complete run.
        for name in [key.GetName() for key in file.GetListOfKeys()]:
            obj=file.Get(name)
            if obj.InheritsFrom('TEfficiency'):
                obj.SetWeight(1.)
                if obj.Write(name,ROOT.TObject.kOverwrite)<=0: raise RuntimeError('efficiency finalization failed')
        for name,text in {'config':original_config, 'photonjet_provenance':json.dumps(info,allow_nan=False),
                          'input_manifest':''.join(f"{p['original_part_index']} {p['source_path']}\n" for p in info['inputs'])}.items():
            file.Delete(name+';*')
            if ROOT.TObjString(text).Write(name,ROOT.TObject.kOverwrite)<=0: raise RuntimeError('metadata write failed')
        file.Delete('h_pj_input_parts;*')
        parts=ROOT.TH1D('h_pj_input_parts','Input completeness',3,0,3)
        parts.SetDirectory(0)
        for i,name in enumerate(('expected','processed','failed'),1): parts.GetXaxis().SetBinLabel(i,name)
        parts.SetBinContent(1,len(info['inputs']));parts.SetBinContent(2,len(info['inputs']))
        if parts.Write('h_pj_input_parts',ROOT.TObject.kOverwrite)<=0: raise RuntimeError('accounting write failed')
        if file.TestBit(ROOT.TFile.kWriteError): raise RuntimeError('ROOT write error')
    finally:
        file.Close()


def check_root(ROOT,path,kind):
    file=ROOT.TFile.Open(str(path))
    if not file or file.IsZombie() or file.TestBit(ROOT.TFile.kRecovered): raise RuntimeError('invalid ROOT output: '+str(path))
    try:
        from ruamel.yaml import YAML
        provenance=file.Get('photonjet_provenance')
        if not provenance: raise RuntimeError('missing merged provenance')
        info=YAML(typ='safe').load(str(provenance.GetString().Data()))
        if (info.get('format_version')!=PROVENANCE_FORMAT_VERSION or info.get('backend')!=BACKEND or
                info.get('schema_version')!=SCHEMA or info.get('input_verification')!=INPUT_VERIFICATION):
            raise RuntimeError('incompatible merged provenance or input verification policy')
        layout=info['layout']
        required=['config','photonjet_provenance','input_manifest','h_pj_input_parts']
        one=layout['n_centrality']==1 and layout['n_eta']==1
        for icent in range(layout['n_centrality']):
            for ieta in range(layout['n_eta']):
                suffix='_0' if one else f'_cent{icent}_eta{ieta}'
                bases=('response_matrix_full','h_response_full','h_pT_truth_response','h_pT_reco_response') if kind=='response' else ('eff_all_eta' if one else 'eff_all','h_truth_pT','h_tight_iso_cluster')
                required.extend(base+suffix for base in bases)
        if kind=='main': required.append('h_pj_cutflow')
        for name in required:
            if not file.Get(name): raise RuntimeError('missing merged object: '+name)
        import math
        keys=[key.GetName() for key in file.GetListOfKeys()]
        if len(keys)!=len(set(keys)): raise RuntimeError('duplicate ROOT object cycles after merge')
        cells=layout['n_centrality']*layout['n_eta']
        expected=ROOT.PJ.ExpectedRootObjectCount(cells,kind=='response',bool(file.Get('mbd_eff_provenance')))
        if len(keys)!=expected: raise RuntimeError('incomplete merged histogram inventory')
        for key in file.GetListOfKeys():
            obj=file.Get(key.GetName())
            histograms=[obj] if obj.InheritsFrom('TH1') else [obj.GetPassedHistogram(),obj.GetTotalHistogram()] if obj.InheritsFrom('TEfficiency') else []
            if obj.InheritsFrom('RooUnfoldResponse'):
                histograms=[getattr(obj,method)() for method in ('Hresponse','Htruth','Hmeasured','Hfakes')]
                if obj.HasFakes(): raise RuntimeError('matched-only response acquired fake entries')
            for hist in histograms:
                if not math.isfinite(hist.GetEntries()): raise RuntimeError('nonfinite histogram entries')
                for i in range(hist.GetNcells()):
                    if not math.isfinite(hist.GetBinContent(i)) or not math.isfinite(hist.GetBinError(i)):
                        raise RuntimeError('nonfinite merged histogram: '+key.GetName())
        return str(file.GetUUID().AsString())
    finally:
        file.Close()


MERGE_SOURCES=('scripts/workflow.py','scripts/inputs.py','scripts/CheckChunks.C','histmakers/InjectMbdEff.C')


def merge(args):
    cfg,plan_path,plan=load_plan(args.config,args.chunk_dir)
    products=args.products.split(',') if args.products else list(plan['products'])
    if len(set(products))!=len(products): raise ValueError('duplicate merge product')
    ROOT=root_module('scripts/CheckChunks.C','histmakers/PhotonJetHistograms.h','histmakers/InjectMbdEff.C')
    source_hashes={str(PJ/p):md5(PJ/p) for p in MERGE_SOURCES}
    staged=[]
    with ExitStack() as cleanup:
        for product in products:
            jobs=plan['products'][product]['jobs']
            chunks=[output_names(cfg,product,j['tag']) for j in jobs]
            names=output_names(cfg,product)
            source=cfg['photonjet'].get('external_mbd_eff_file','') if product=='sim_signal' else ''
            if source: require_single_cell(plan['layout'],'MBD injection')
            if source and (output_path(source)==names['main'] or
                           (Path(source).exists() and names['main'].exists() and Path(source).samefile(names['main']))):
                raise ValueError('MBD source and final target must differ')
            for kind,final in names.items():
                expected={entry[kind] for entry in chunks}
                observed={output_path(p) for p in (final.parent/'chunks').glob(f'{final.stem}.{product}_*.root')}
                if expected!=observed: raise ValueError(f'missing or unplanned {product}/{kind} chunk outputs')
            with tempfile.NamedTemporaryFile('w',suffix='.list') as stream:
                stream.write(''.join(str(chunk['main'])+'\n' for chunk in chunks));stream.flush()
                ROOT.PJ.ValidateChunks(stream.name,str(args.config),product,str(plan_path))
            records=[checked_generation(ROOT,chunk['main'],args.config,product) for chunk in chunks]
            for record,job in zip(records,jobs): require_job_plan(record,plan_path,[job])
            info=deepcopy(records[0]);info['run_id']=str(ROOT.TUUID().AsString());info['stage']='merged'
            info['job_plan']['assignments']=[assignment for record in records for assignment in record['job_plan']['assignments']]
            info['inputs']=sorted([p for record in records for p in record['inputs']],key=lambda p:p['original_part_index'])
            info['dependencies']={key:value for record in records for key,value in record['dependencies'].items()}
            info['dependencies'][str(plan_path)]=md5(plan_path)
            info['dependencies'].update(source_hashes)
            info['merge_inputs']=[dict(main=str(chunk['main']),md5=md5(chunk['main']),run_id=record['run_id']) for chunk,record in zip(chunks,records)]
            info['workers']=1
            for name in ('events_processed','parts_processed','event_loops','analysis_seconds'):
                info[name]=sum(r[name] for r in records)
            info['selected_events_by_centrality']=[sum(r['selected_events_by_centrality'][i] for r in records)
                                                  for i in range(info['layout']['n_centrality'])]
            info['chunk_execution']=[{key:r[key] for key in ('run_id','workers','events_processed','parts_processed','analysis_seconds')} for r in records]
            temporary={}
            for kind,final in names.items():
                final.parent.mkdir(parents=True,exist_ok=True)
                directory=Path(cleanup.enter_context(tempfile.TemporaryDirectory(prefix=final.stem+'.merge.',dir=final.parent)))
                temporary[kind]=directory/final.name
                merge_root(ROOT,[chunk[kind] for chunk in chunks],temporary[kind])
            if product=='sim_signal':
                file=ROOT.TFile.Open(str(temporary['response']),'UPDATE')
                try:
                    ROOT.PJ.FinalizeMergedResponses(file,ROOT.PJ.LoadCuts(ROOT.PJ.ReadConfig(str(args.config)).yaml))
                finally:
                    file.Close()
                source=cfg['photonjet'].get('external_mbd_eff_file','')
                if source and info['dependencies'].get(str(Path(source).resolve()))!=md5(source):
                    raise ValueError('MBD source changed since job planning')
            for path in temporary.values(): write_metadata(ROOT,path,Path(args.config).read_text(),info)
            if product=='sim_signal' and cfg['photonjet'].get('external_mbd_eff_file',''):
                ROOT.PJ.InjectCorrection(str(args.config),str(temporary['main']))
            marker=dict(format_version=COMPLETION_FORMAT_VERSION,backend=BACKEND,input_verification=INPUT_VERIFICATION,
                        layout=info['layout'],stage='merged',product=product,run_id=info['run_id'],
                        original_part_ids=[p['original_part_id'] for p in info['inputs']],files=[])
            for kind,path in temporary.items():
                uuid=check_root(ROOT,path,kind)
                marker['files'].append(dict(kind=kind,path=str(names[kind]),md5=md5(path),bytes=path.stat().st_size,ROOT_UUID=uuid))
            for path,digest in info['dependencies'].items():
                if md5(path)!=digest: raise ValueError('dependency changed during merge: '+path)
            marker_temp=temporary['main'].parent/'completion.yaml';save(marker_temp,marker)
            staged.append((product,names,temporary,marker_temp))
        # Every merge and MBD import has succeeded before invalidating any prior
        # generation. A partial rename leaves no acceptable mixed generation.
        for _,names,_,_ in staged: Path(str(names['main'])+'.complete.yaml').unlink(missing_ok=True)
        for _,names,temporary,_ in staged:
            for kind,path in temporary.items(): os.replace(path,names[kind])
        for _,names,_,marker in staged: os.replace(marker,Path(str(names['main'])+'.complete.yaml'))
        for product,names,_,_ in staged:
            checked_generation(ROOT,names['main'],args.config,product)
            print(f"[merge] published complete {product}: {names['main']}",flush=True)


def checked_merged(ROOT, cfg, plan, config_path, allow_subset=False, plan_path=None):
    records={}
    for product in PRODUCTS:
        record=checked_generation(ROOT,output_names(cfg,product)['main'],config_path,product)
        assignment=plan['products'][product]
        source_plan=plan_path or Path(assignment['jobs'][0]['input']).parent/'jobs.json'
        require_job_plan(record,source_plan,assignment['jobs'])
        if record.get('stage')!='merged' or sorted(p['original_part_index'] for p in record['inputs'])!=sorted(assignment['indices']):
            raise ValueError('merged output does not cover the current plan: '+product)
        if assignment['scope']!='full' and not allow_subset:
            raise ValueError('yield requires full original coverage; --allow-subset is for explicitly bounded validation')
        records[product]=record
    return records


def check_final(ROOT,path):
    import math
    file=ROOT.TFile.Open(str(path))
    if not file or file.IsZombie() or file.TestBit(ROOT.TFile.kRecovered): raise RuntimeError('invalid yield output: '+str(path))
    try:
        for name in ('h_unfold_sub_result','h_unfold_sub_result_woeff','gpurity','gpurity_leak','g_mbd_eff'):
            if not file.Get(name): raise RuntimeError('missing yield output: '+name)
        for name in ('h_unfold_sub_result','h_unfold_sub_result_woeff'):
            hist=file.Get(name)
            for i in range(hist.GetNcells()):
                if not math.isfinite(hist.GetBinContent(i)) or not math.isfinite(hist.GetBinError(i)):
                    raise RuntimeError('nonfinite final yield: '+name)
        return str(file.GetUUID().AsString())
    finally:
        file.Close()


def yield_outputs(args):
    from ruamel.yaml import YAML
    cfg,plan_path,plan=load_plan(args.config,args.chunk_dir)
    require_single_cell(plan['layout'],'CalculatePhotonYield')
    if cfg['photonjet']['system']!='pp': raise ValueError('CalculatePhotonYield is the pp yield/closure chain')
    ROOT=root_module('histmakers/PhotonJetIO.h')
    records=checked_merged(ROOT,cfg,plan,args.config,args.allow_subset,plan_path)
    variation=label(str(cfg['output']['var_type']))
    prefix=output_path(cfg['output']['final_outfile'])
    finals={kind:Path(str(prefix)+'_'+variation+suffix+'.root') for kind,suffix in (('data',''),('mc','_mc'))}
    marker=Path(str(finals['data'])+'.complete.yaml')
    input_hashes={str(path):md5(path) for product in PRODUCTS for path in output_names(cfg,product).values()}
    input_hashes.update({str(output_names(cfg,p)['main'])+'.complete.yaml':md5(str(output_names(cfg,p)['main'])+'.complete.yaml') for p in PRODUCTS})
    dependencies={str(path):md5(path) for path in (Path(args.config),plan_path,PJ/'support/CalculatePhotonYield.C',PJ/'support/CrossSectionWeights.h',Path(__file__))}
    purity=cfg['analysis'].get('mc_purity_correction_file','')
    if cfg['analysis'].get('mc_purity_correction',0) and purity:
        dependencies[str(Path(purity).resolve())]=md5(purity)
    logdir=Path(args.log_dir);logdir.mkdir(parents=True,exist_ok=True)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=prefix.name+'.yield.',dir=prefix.parent) as directory:
        stage=Path(directory); runtime=deepcopy(cfg)
        for key in ('data_outfile','eff_outfile','response_outfile'):
            runtime['output'][key]=str(output_path(runtime['output'][key]))
        runtime['output']['final_outfile']=str(stage/'final')
        runtime_path=stage/'runtime.yaml'
        with runtime_path.open('w') as stream: YAML().dump(runtime,stream)
        result=dict(format_version=1,stage='yield',run_id=str(ROOT.TUUID().AsString()),product='pp_data_and_mc_closure',
                    config_md5=md5(args.config),inputs=input_hashes,dependencies=dependencies,
                    histogram_generations={p:r['run_id'] for p,r in records.items()},files=[])
        outputs={}
        for kind,suffix in (('data',''),('mc','_mc')):
            call=str(PJ/'support/CalculatePhotonYield.C')+'('+json.dumps(str(runtime_path))+','+('true' if kind=='mc' else 'false')+')'
            log=logdir/f'photonjet_yield_{variation}{suffix}.log'
            with log.open('w') as stream:
                subprocess.run(['root','-l','-b','-q',call],cwd=PJ,stdout=stream,stderr=subprocess.STDOUT,check=True)
            path=stage/f'final_{variation}{suffix}.root';check_final(ROOT,path)
            file=ROOT.TFile.Open(str(path),'UPDATE')
            try:
                for name,value in dict(config=Path(args.config).read_text(),runtime_config=runtime_path.read_text(),
                                       photonjet_yield_provenance=json.dumps({k:v for k,v in result.items() if k!="files"},allow_nan=False)).items():
                    file.Delete(name+';*')
                    if ROOT.TObjString(value).Write(name,ROOT.TObject.kOverwrite)<=0: raise RuntimeError('yield provenance write failed')
                if file.TestBit(ROOT.TFile.kWriteError): raise RuntimeError('yield ROOT write error')
            finally: file.Close()
            outputs[kind]=path
            result['files'].append(dict(kind=kind,path=str(finals[kind]),md5=md5(path),bytes=path.stat().st_size,ROOT_UUID=check_final(ROOT,path)))
        checked_merged(ROOT,cfg,plan,args.config,args.allow_subset,plan_path)
        for path,digest in {**input_hashes,**dependencies}.items():
            if md5(path)!=digest: raise ValueError('yield dependency changed: '+path)
        temporary=stage/'completion.yaml';save(temporary,result)
        marker.unlink(missing_ok=True)
        for kind,path in outputs.items(): os.replace(path,finals[kind])
        os.replace(temporary,marker)
    print('[yield] published data yield and MC closure: '+str(marker),flush=True)


def local(args):
    require_serial(args.threads)
    _,_,plan=load_plan(args.config,args.chunk_dir)
    products=args.products.split(',') if args.products else list(plan['products'])
    for product in products:
        for job in plan['products'][product]['jobs']:
            worker(argparse.Namespace(config=args.config,product=product,input=job['input'],tag=job['tag'],threads=args.threads))


def status(args):
    cfg,plan_path,plan=load_plan(args.config,args.chunk_dir)
    ROOT=root_module('histmakers/PhotonJetIO.h')
    for product,assignment in plan['products'].items():
        complete=0
        for job in assignment['jobs']:
            main=output_names(cfg,product,job['tag'])['main']
            try:
                info=checked_generation(ROOT,main,args.config,product)
                require_job_plan(info,plan_path,[job]);complete+=1
            except Exception as error: print(f'[status] {product}/{job["tag"]}: incomplete ({error})')
        print(f'[status] {product}: {complete}/{len(assignment["jobs"])} checked complete generations',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step',choices=('worker','local','merge','yield','status','check'))
    parser.add_argument('config',type=Path)
    parser.add_argument('--chunk-dir',default=os.environ.get('PJ_CHUNK_DIR',str(PJ/'chunks')))
    parser.add_argument('--log-dir',default=os.environ.get('PJ_LOG_DIR',str(PJ/'logs')))
    parser.add_argument('--products')
    parser.add_argument('--product',choices=PRODUCTS)
    parser.add_argument('--input',type=Path)
    parser.add_argument('--tag')
    parser.add_argument('--threads',type=int,default=int(os.environ.get('PJ_THREADS','1')))
    parser.add_argument('--allow-subset',action='store_true',help='allow explicitly planned subsets for bounded yield validation')
    args=parser.parse_args()
    require_serial(args.threads)
    args.config=args.config.resolve();args.chunk_dir=Path(args.chunk_dir).resolve();args.log_dir=Path(args.log_dir).resolve()
    if args.input: args.input=args.input.resolve()
    os.chdir(PJ)
    if args.step=='worker' and not all((args.input,args.tag,args.product)):
        parser.error('worker needs --product, --input and --tag')
    if args.step=='check':
        cfg,plan_path,plan=load_plan(args.config,args.chunk_dir)
        checked_merged(root_module('histmakers/PhotonJetIO.h'),cfg,plan,args.config,args.allow_subset,plan_path)
        print('[check] all histogram generations and planned coverage are valid')
    else:
        {'worker':worker,'local':local,'merge':merge,'yield':yield_outputs,'status':status}[args.step](args)


if __name__=='__main__':
    main()
