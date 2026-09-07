#!/usr/bin/env python3
"""Generate explicit staged inputs. This program never submits or runs simulations."""
from __future__ import annotations
import argparse, copy, hashlib, json, math, os, re, shutil, subprocess, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
CONTRACT=json.loads((ROOT/'campaign_contract.json').read_text())
BASELINE=CONTRACT['model_baseline_sha']
EXEC_PLACEHOLDER=CONTRACT['execution_source_sha_policy']['planning_placeholder']
IMAGE_PLACEHOLDER=CONTRACT['digest_policy']['planning_placeholder']
SEEDS=CONTRACT['seeds']
SHA40=re.compile(r'[0-9a-f]{40}')
IMAGE_RE=re.compile(CONTRACT['execution']['image_digest_pattern'])
CORE_FILES=['campaign_contract.json','campaign_decision_record.json','CURSOR_HANDOFF.md','prepare.py','preflight.py','analyze.py','aws_commands.py','README.md','AWS_HANDOFF.md','COMPLETION_NOTE.md','prepare_and_preflight.sh','tests/test_package.py']

def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    data=json.dumps(obj,indent=2,sort_keys=True)+'\n'
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f: f.write(data); f.flush(); os.fsync(f.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def git(*args): return subprocess.check_output(['git','-C',str(REPO),*args],text=True,stderr=subprocess.STDOUT).strip()
def validate_planning_baseline():
    manifest=(ROOT/CONTRACT['model_baseline_snapshot_manifest']).resolve()
    if manifest.exists():
        try: observed=json.loads(manifest.read_text()).get('commit_sha')
        except Exception as e: raise SystemExit(f'REFUSED: model-baseline retrieval manifest is unreadable: {e}')
        if observed!=BASELINE: raise SystemExit(f'REFUSED: retrieval manifest SHA {observed} does not match model_baseline_sha {BASELINE}')
        return
    try: head=git('rev-parse','HEAD'); rc=subprocess.call(['git','-C',str(REPO),'merge-base','--is-ancestor',BASELINE,head],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    except Exception as e: raise SystemExit(f'REFUSED: planning requires a matching baseline retrieval manifest or Git ancestry: {e}')
    if rc!=0: raise SystemExit(f'REFUSED: model_baseline_sha {BASELINE} is not an ancestor of planning checkout HEAD {head}')
def validate_deployment_identity(execution_sha,image_digest):
    if not SHA40.fullmatch(execution_sha or ''): raise SystemExit('REFUSED: --deployment requires --execution-source-sha as a full lowercase 40-hex commit')
    if not IMAGE_RE.fullmatch(image_digest or ''): raise SystemExit('REFUSED: --deployment requires --image-digest sha256:<64 lowercase hex>')
    try: inside=git('rev-parse','--is-inside-work-tree'); head=git('rev-parse','HEAD')
    except Exception as e: raise SystemExit(f'REFUSED: deployment generation requires a normal git checkout: {e}')
    if inside!='true' or head!=execution_sha: raise SystemExit(f'REFUSED: --execution-source-sha {execution_sha} does not match git HEAD {head}')
    try: subprocess.check_call(['git','-C',str(REPO),'merge-base','--is-ancestor',BASELINE,head],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    except Exception: raise SystemExit(f'REFUSED: audited model baseline {BASELINE} is not an ancestor of execution source {head}')
    rel=ROOT.relative_to(REPO)
    for name in CORE_FILES:
        try: git('ls-files','--error-unmatch','--',str(rel/name))
        except Exception: raise SystemExit(f'REFUSED: {rel/name} is not committed at execution_source_sha; commit the campaign package first')
    dirty=git('status','--porcelain','--',*[str(rel/name) for name in CORE_FILES])
    if dirty: raise SystemExit('REFUSED: provenance-bearing campaign files differ from HEAD; commit them before deployment generation:\n'+dirty)

def common(execution_sha):
    return {
      '_comment':['GPU-only receptor-selection campaign v1.',f'Audited design baseline: {BASELINE}.',f'Execution source: {execution_sha}.','Runtime /run_provenance/git_sha must match execution_source_sha.','AWS image digest is enforced by handoff/preflight, outside parser config.'],
      'total_time':21600.0,'bio_dt':60.0,'output_interval':300.0,'seed':0,
      'domain_x':1e-4,'domain_y':1e-4,'domain_z':1e-4,'grid_dx':2e-6,'mucus_thickness':1e-4,
      'radial_turnover':1e9,'distal_transit':1e9,'peristaltic_enabled':False,'crypts_enabled':False,'motility.enabled':False,
      'carbon_z_gradient':False,'carbon.boundary_conc':0.05,'metabolism.uptake_limit':'delivery','oxygen.k_ROS':0.0,
      'dysbiosis_threshold':1e10,'gpu_enabled':True,'gpu_device_id':0,'chemistry.toxin_evaluation':'grid','chemistry.toxin_lumping':'per_receptor',
      'initial_population.placement':'z_slab','initial_population.z_min':0.0,'initial_population.z_max':1e-4,
      'fixes':['metabolism','bacteriocin','receptor','mechanics'],'hdf5_file':'output.h5',
      'hdf5':{'enabled':True,'compression':'gzip','compression_level':4,'schedule':{'summary':1,'agents':10,'grid':0,'lineage':0,'genome':0,'provenance':10,'grid_species':[]}}
    }
def strain(t,plasmids=None,receptor=None):
    x={'type':t,'count':60,'mu_max':5.5e-4,'plasmids':plasmids or [],'conjugative':False}
    if receptor is not None: x['receptor_expression']={'BtuB':receptor}
    return x
def annotate(c,stage,run_id,arm,seed,axes,gate_locked,execution_sha):
    c['_campaign']={'campaign_id':CONTRACT['campaign_id'],'stage':stage,'run_id':run_id,'arm':arm,'paired_seed':seed,'axes':axes,'model_baseline_sha':BASELINE,'execution_source_sha':execution_sha,'gate_locked':gate_locked}
    return c
def make(stage,run_id,arm,seed,axes,strains,execution_sha,updates=None,grid=False,gate_locked=True):
    c=common(execution_sha); c['seed']=seed; c['initial_strains']=strains
    for k,v in (updates or {}).items(): c[k]=v
    if grid: c['hdf5']['schedule']['grid']=60; c['hdf5']['schedule']['grid_species']=['bacteriocin_BtuB']
    return annotate(c,stage,run_id,arm,seed,axes,gate_locked,execution_sha)
def planned_runs(prom,execution_sha):
    out={s:[] for s in 'QABCD'}
    base={'kd_corrinoid_btuB':1e-6,'kd_colicinE_btuB':5e-7,'b12_initial_conc':1e-3,'bacteriocin.mucin_charge.amplitude':60}
    for rep in (1,2):
      rid=f'Q_repeat{rep}_s{SEEDS[0]}'; out['Q'].append((rid,make('Q',rid,'qualification',SEEDS[0],{'repeat':rep},[strain(1,['ColE1']),strain(2)],execution_sha,base,gate_locked=False)))
    for kd in CONTRACT['axes']['A']['kd_corrinoid_btuB_mol_m3']:
      for arm in ('producer','null'):
       for seed in SEEDS:
        rid=f'A_kd{kd:.0e}_{arm}_s{seed}'; ss=[strain(1,['ColE1'] if arm=='producer' else []),strain(2)]
        out['A'].append((rid,make('A',rid,arm,seed,{'kd_corrinoid_btuB_mol_m3':kd},ss,execution_sha,{**base,'kd_corrinoid_btuB':kd})))
    kd=prom['selected_kd_corrinoid_btuB_mol_m3']
    for b12 in CONTRACT['axes']['B']['b12_initial_conc_mol_m3']:
     for seed in SEEDS:
      rid=f'B_b12{b12:.0e}_competition_s{seed}'; ss=[strain(1,[],1.0),strain(2,[],0.0)]
      out['B'].append((rid,make('B',rid,'btuB_normal_vs_null',seed,{'selected_kd_corrinoid_btuB_mol_m3':kd,'b12_initial_conc_mol_m3':b12},ss,execution_sha,{**base,'kd_corrinoid_btuB':kd,'b12_initial_conc':b12})))
    b12=prom['selected_b12_initial_conc_mol_m3']; amp0=prom['selected_mucin_charge_amplitude']; fixed={**base,'kd_corrinoid_btuB':kd,'b12_initial_conc':b12}
    for amp in CONTRACT['axes']['C']['bacteriocin.mucin_charge.amplitude']:
     for seed in SEEDS:
      rid=f'C_amp{amp}_producer_s{seed}'; out['C'].append((rid,make('C',rid,'producer',seed,{'amplitude':amp,'selected_kd':kd,'selected_b12':b12},[strain(1,['ColE1']),strain(2)],execution_sha,{**fixed,'bacteriocin.mucin_charge.amplitude':amp},grid=True)))
    for seed in SEEDS:
      rid=f'C_shared_null_amp{amp0}_s{seed}'; out['C'].append((rid,make('C',rid,'shared_null',seed,{'amplitude':amp0,'selected_kd':kd,'selected_b12':b12},[strain(1),strain(2)],execution_sha,{**fixed,'bacteriocin.mucin_charge.amplitude':amp0},grid=True)))
    # Twelve ColE1-carrier arms; P=0 remains a carrier and is not the null.
    for target,p in zip(CONTRACT['axes']['D']['realized_lysis_target_per_generation'],CONTRACT['axes']['D']['sos_lysis_prob']):
     for seed in SEEDS:
      rid=f'D_target{target:.3f}_producer_s{seed}'; updates={**fixed,'bacteriocin.mucin_charge.amplitude':amp0,'sos_basal_rate':0.0,'sos_lysis_prob':p}
      out['D'].append((rid,make('D',rid,'producer',seed,{'nominal_target_per_generation':target,'sos_lysis_prob':p,'selected_amplitude':amp0},[strain(1,['ColE1']),strain(2)],execution_sha,updates)))
    # Three formal same-revision controls, generated once per seed (not crossed with lysis targets).
    for seed in SEEDS:
      rid=f'D_plasmid_free_null_s{seed}'; updates={**fixed,'bacteriocin.mucin_charge.amplitude':amp0,'sos_basal_rate':0.0,'sos_lysis_prob':0.0}
      axes={'control':'plasmid_free_null','selected_amplitude':amp0,'selected_kd':kd,'selected_b12':b12}
      out['D'].append((rid,make('D',rid,'plasmid_free_null',seed,axes,[strain(1),strain(2)],execution_sha,updates)))
    return out
def main():
 p=argparse.ArgumentParser(); p.add_argument('--deployment',action='store_true'); p.add_argument('--execution-source-sha'); p.add_argument('--image-digest',default=IMAGE_PLACEHOLDER); p.add_argument('--selected-kd',type=float); p.add_argument('--selected-b12',type=float); p.add_argument('--selected-amplitude',type=float); p.add_argument('--clean',action='store_true'); a=p.parse_args()
 execution_sha=a.execution_source_sha or EXEC_PLACEHOLDER
 if a.deployment: validate_deployment_identity(execution_sha,a.image_digest)
 else:
  validate_planning_baseline()
  if a.execution_source_sha or a.image_digest!=IMAGE_PLACEHOLDER: raise SystemExit('REFUSED: real execution SHA/image digest require --deployment; planning generation uses explicit placeholders')
 prom=copy.deepcopy(CONTRACT['planning_promotions'])
 if a.selected_kd is not None: prom['selected_kd_corrinoid_btuB_mol_m3']=a.selected_kd; prom['status']='USER_RECORDED_PROMOTION'
 if a.selected_b12 is not None: prom['selected_b12_initial_conc_mol_m3']=a.selected_b12; prom['status']='USER_RECORDED_PROMOTION'
 if a.selected_amplitude is not None: prom['selected_mucin_charge_amplitude']=a.selected_amplitude; prom['status']='USER_RECORDED_PROMOTION'
 gen=ROOT/'generated'
 if a.clean and gen.exists(): shutil.rmtree(gen)
 allruns=planned_runs(prom,execution_sha); campaign=[]
 provenance={'model_baseline_sha':BASELINE,'execution_source_sha':execution_sha,'container_image_digest':a.image_digest}
 for stage,runs in allruns.items():
  sd=gen/f'stage_{stage}'; entries=[]
  for idx,(rid,cfg) in enumerate(runs):
   ip=sd/'jobs'/str(idx)/'input.json'; dump(ip,cfg)
   e={'array_index':idx,'run_id':rid,'stage':stage,'input_relpath':str(ip.relative_to(ROOT)),'input_sha256':digest(ip),'output_relpath':f'generated/stage_{stage}/results/{idx}/output.h5.gz'}; entries.append(e); campaign.append(e)
  manifest={'schema_version':2,'campaign_id':CONTRACT['campaign_id'],'stage':stage,**provenance,'mpi_ranks':1,'gpu_required':True,'attempt_timeout_s':7200,'gate_requires':CONTRACT['stages'][stage]['requires'],'gate_locked':CONTRACT['stages'][stage]['requires'] is not None,'promotion_inputs':prom,'runs':entries}
  dump(sd/'manifest.json',manifest)
 dump(gen/'campaign_manifest.json',{'schema_version':2,'campaign_id':CONTRACT['campaign_id'],**provenance,'promotion_inputs':prom,'runs':campaign})
 print('generated', {s:len(v) for s,v in allruns.items()}, 'total',len(campaign),'model_baseline_sha',BASELINE,'execution_source_sha',execution_sha,'digest',a.image_digest)
if __name__=='__main__': main()
