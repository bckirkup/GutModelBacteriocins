#!/usr/bin/env python3
"""Fail-closed structural preflight for generated campaign inputs/manifests."""
from __future__ import annotations
import argparse, hashlib, json, math, re, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent; REPO=ROOT.parents[1]
MUCIN_AMPLITUDE_KEY='bacteriocin.mucin_charge.amplitude'
OXYGEN_K_ROS_KEY='oxygen.k_ROS'
C=json.loads((ROOT/'campaign_contract.json').read_text()); G=ROOT/'generated'; ERR=[]; WARN=[]
BASELINE=C['model_baseline_sha']; EXEC_PLACEHOLDER=C['execution_source_sha_policy']['planning_placeholder']
SHA40=re.compile(r'[0-9a-f]{40}')
def fail(x): ERR.append(x)
def approx(a,b,tol=1e-12): return math.isclose(float(a),float(b),rel_tol=tol,abs_tol=tol*max(1,abs(float(b))))
def strip_meta(x):
 if isinstance(x,dict): return {k:strip_meta(v) for k,v in x.items() if not k.startswith('_')}
 if isinstance(x,list): return [strip_meta(v) for v in x]
 return x
def git(*args): return subprocess.check_output(['git','-C',str(REPO),*args],text=True,stderr=subprocess.STDOUT).strip()
def git_identity():
 try:
  if git('rev-parse','--is-inside-work-tree')!='true': return None
  return git('rev-parse','HEAD')
 except Exception: return None
def verify_baseline(head,deployment=False):
 if deployment:
  if not head: fail('deployment cannot verify model baseline without git HEAD'); return
  try:
   subprocess.check_call(['git','-C',str(REPO),'cat-file','-e',BASELINE+'^{commit}'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
   rc=subprocess.call(['git','-C',str(REPO),'merge-base','--is-ancestor',BASELINE,head],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
   if rc!=0: fail(f'audited model baseline {BASELINE} is not an ancestor of execution checkout HEAD {head}')
  except Exception as e: fail(f'cannot verify model baseline ancestry in deployment checkout: {e}')
  return
 manifest_path=(ROOT/C['model_baseline_snapshot_manifest']).resolve()
 if manifest_path.exists():
  try: observed=json.loads(manifest_path.read_text()).get('commit_sha')
  except Exception as e: fail(f'model-baseline retrieval manifest unreadable: {e}'); return
  if observed!=BASELINE: fail(f'model baseline mismatch in retrieval manifest: {observed} != {BASELINE}')
  return
 if head:
  try:
   subprocess.check_call(['git','-C',str(REPO),'cat-file','-e',BASELINE+'^{commit}'],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
   rc=subprocess.call(['git','-C',str(REPO),'merge-base','--is-ancestor',BASELINE,head],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
   if rc!=0: fail(f'audited model baseline {BASELINE} is not an ancestor of execution checkout HEAD {head}')
  except Exception as e: fail(f'cannot verify model baseline in git checkout: {e}')
 else: fail('cannot verify audited model baseline: no matching retrieval manifest and no git checkout')
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--deployment',action='store_true',help='require committed execution source at git HEAD and immutable image digest'); ap.add_argument('--execution-source-sha',help='required in deployment mode; must equal git HEAD and generated provenance'); a=ap.parse_args()
 head=git_identity(); verify_baseline(head,a.deployment)
 if a.deployment:
  if not a.execution_source_sha: fail('--deployment requires --execution-source-sha')
  elif not SHA40.fullmatch(a.execution_source_sha): fail('--execution-source-sha must be full lowercase 40-hex')
  if head is None: fail('deployment preflight requires a normal git checkout; retrieval-manifest-only/tarball provenance is planning-only')
  elif a.execution_source_sha!=head: fail(f'execution source SHA does not match git HEAD: {a.execution_source_sha} != {head}')
  else:
   core=['campaign_contract.json','campaign_decision_record.json','CURSOR_HANDOFF.md','prepare.py','preflight.py','analyze.py','aws_commands.py','README.md','AWS_HANDOFF.md','COMPLETION_NOTE.md','prepare_and_preflight.sh','tests/test_package.py']; rel=ROOT.relative_to(REPO)
   for name in core:
    try: git('ls-files','--error-unmatch','--',str(rel/name))
    except Exception: fail(f'deployment package file is not committed at execution_source_sha: {rel/name}')
   try: dirty=git('status','--porcelain','--',*[str(rel/name) for name in core])
   except Exception as e: fail(f'cannot inspect deployment package working-tree state: {e}'); dirty=''
   if dirty: fail('provenance-bearing campaign files differ from execution_source_sha:\n'+dirty)
 elif a.execution_source_sha: fail('--execution-source-sha is only accepted with --deployment')
 try:
  decision=json.loads((ROOT/C['authoritative_decision_record']).read_text()); by_gate={x['gate']:x for x in decision['decisions']}
  if by_gate['A_pass'].get('execution_source_sha')!='b884cc54b1c0684c079a90195c031e388fe70534' or not approx(by_gate['A_pass'].get('selected_kd_corrinoid_btuB_mol_m3'),1e-4): fail('authoritative Stage A decision drift')
  if by_gate['B_pass'].get('execution_source_sha')!='b884cc54b1c0684c079a90195c031e388fe70534' or not approx(by_gate['B_pass'].get('selected_b12_initial_conc_mol_m3'),1e-3): fail('authoritative Stage B decision drift')
  if by_gate['C_population_gate'].get('execution_source_sha')!='3f176b26c0d18a22a61db218e56106b1355b781e' or by_gate['C_population_gate'].get('amplitudes_tested')!=[0,15,60]: fail('authoritative ecological Stage C decision drift')
  if by_gate['PR416_intrinsic_transport_gate'].get('status')!='INVALID_FAILED_SUPERSEDED': fail('PR416 gate must remain superseded')
  if by_gate['single_source_transport_assay_v1'].get('status') not in (
      'PENDING', 'PENDING_RESUBMIT', 'PENDING_ADAPTIVE_RUNS'):
    fail('revised transport assay must remain pending until a recorded pass')
  if by_gate['D_release'].get('status')!='BLOCKED' or by_gate['D_release'].get('conditional_selected_mucin_charge_amplitude')!=15: fail('Stage D handoff must remain blocked/conditional')
 except Exception as e: fail(f'authoritative decision record missing/unreadable: {e}')
 expected={'Q':2,'A':24,'B':12,'C':12,'D':15}; configs={}; ids=[]; signatures=[]; image_digests=set(); execution_shas=set(); baseline_shas=set()
 try: campaign_manifest=json.loads((G/'campaign_manifest.json').read_text())
 except Exception as e: fail(f'campaign manifest missing/unreadable: {e}'); campaign_manifest={}
 # Source-backed audit of every deliberate scalar key. Nested HDF5/strain objects are parsed by config_json.cpp.
 try: parser_text=(REPO/'src/io/input_parser.cpp').read_text()+(REPO/'src/io/config_json.cpp').read_text()
 except Exception as e: fail(f'cannot read execution-source parser files: {e}'); parser_text=''
 required_exact=['total_time','bio_dt','output_interval','seed','domain_x','domain_y','domain_z','grid_dx','mucus_thickness','radial_turnover','distal_transit','peristaltic_enabled','crypts_enabled','motility.enabled','carbon_z_gradient','carbon.boundary_conc','metabolism.uptake_limit',OXYGEN_K_ROS_KEY,'dysbiosis_threshold','gpu_enabled','gpu_device_id','chemistry.toxin_evaluation','chemistry.toxin_lumping','initial_population.placement','initial_population.z_min','initial_population.z_max','fixes','hdf5_file','kd_corrinoid_btuB','kd_colicinE_btuB','b12_initial_conc',MUCIN_AMPLITUDE_KEY,'sos_basal_rate','sos_lysis_prob','grid_species','receptor_expression']
 for key in required_exact:
  if f'"{key}"' not in parser_text: fail(f'exact config key not found in execution-source parser: {key}')
 aliases=[('b12.initial_conc','b12_initial_conc','corrinoid.initial_conc','corrinoid_initial_conc'),('kd_b12_btuB','kd_corrinoid_btuB'),('hdf5_file','hdf5.file')]
 for stage,n in expected.items():
  mp=G/f'stage_{stage}'/'manifest.json'
  if not mp.exists(): fail(f'missing {mp.relative_to(ROOT)}'); continue
  m=json.loads(mp.read_text()); image_digests.add(m.get('container_image_digest')); execution_shas.add(m.get('execution_source_sha')); baseline_shas.add(m.get('model_baseline_sha')); runs=m.get('runs',[])
  if m.get('schema_version')!=2: fail(f'Stage {stage}: manifest schema is not 2')
  if len(runs)!=n: fail(f'Stage {stage}: {len(runs)} runs != {n}')
  if m.get('model_baseline_sha')!=BASELINE: fail(f'Stage {stage}: model baseline SHA mismatch')
  expected_exec=a.execution_source_sha if a.deployment else EXEC_PLACEHOLDER
  if m.get('execution_source_sha')!=expected_exec: fail(f"Stage {stage}: execution source SHA {m.get('execution_source_sha')!r} != {expected_exec!r}")
  if m.get('mpi_ranks')!=1 or not m.get('gpu_required'): fail(f'Stage {stage}: GPU/one-rank manifest contract violated')
  if m.get('attempt_timeout_s')>7200: fail(f'Stage {stage}: timeout exceeds 2 h')
  for i,e in enumerate(runs):
   if e.get('array_index')!=i: fail(f'Stage {stage}: non-contiguous array index at {i}')
   p=ROOT/e['input_relpath']
   if not p.exists(): fail(f'missing input {p}'); continue
   raw=p.read_bytes(); sha=hashlib.sha256(raw).hexdigest()
   if sha!=e.get('input_sha256'): fail(f'{p}: digest mismatch')
   cfg=json.loads(raw); rid=e['run_id']; ids.append(rid); configs[rid]=cfg; meta=cfg.get('_campaign',{})
   if meta.get('run_id')!=rid or meta.get('stage')!=stage or meta.get('model_baseline_sha')!=BASELINE or meta.get('execution_source_sha')!=expected_exec: fail(f'{rid}: provenance metadata mismatch')
   for fam in aliases:
    present=[k for k in fam if k in cfg]
    if len(present)>1: fail(f'{rid}: alias collision {present}')
   for k,want in [('gpu_enabled',True),('gpu_device_id',0),('total_time',21600),('bio_dt',60),('domain_x',1e-4),('domain_y',1e-4),('domain_z',1e-4),('grid_dx',2e-6)]:
    got=cfg.get(k); ok=(got is want) if isinstance(want,bool) else (got is not None and approx(got,want))
    if not ok: fail(f'{rid}: {k}={got!r}, expected {want!r}')
   if cfg.get('fixes')!=['metabolism','bacteriocin','receptor','mechanics']: fail(f'{rid}: mechanism/fix axis drift')
   if cfg.get('metabolism.uptake_limit')!='delivery' or cfg.get('chemistry.toxin_evaluation')!='grid': fail(f'{rid}: unsupported GPU chemistry placement request')
   h=cfg.get('hdf5',{}); sch=h.get('schedule',{})
   for k,v in [('summary',1),('agents',10),('provenance',10),('lineage',0),('genome',0)]:
    if sch.get(k)!=v: fail(f'{rid}: HDF5 {k} schedule drift')
   if stage=='C':
    if sch.get('grid')!=60 or sch.get('grid_species')!=['bacteriocin_BtuB']: fail(f'{rid}: Stage C must write only named hourly BtuB toxin grid')
   elif sch.get('grid')!=0 or sch.get('grid_species') not in ([],None): fail(f'{rid}: heavy grid output outside Stage C')
   if len(cfg.get('initial_strains',[]))!=2 or [x.get('count') for x in cfg['initial_strains']]!=[60,60]: fail(f'{rid}: founder axis drift')
   signatures.append((stage,rid,json.dumps(strip_meta(cfg),sort_keys=True)))
 if campaign_manifest:
  for k,want in [('schema_version',2),('model_baseline_sha',BASELINE),('execution_source_sha',a.execution_source_sha if a.deployment else EXEC_PLACEHOLDER)]:
   if campaign_manifest.get(k)!=want: fail(f'campaign manifest {k} mismatch: {campaign_manifest.get(k)!r} != {want!r}')
  if len(campaign_manifest.get('runs',[]))!=65: fail('campaign manifest does not contain 65 runs')
 if len(ids)!=len(set(ids)): fail('run-id collision')
 if len(signatures)!=65: fail(f'total explicit runs {len(signatures)} != 65')
 q=[s for st,_,s in signatures if st=='Q']
 if len(set(q))!=1: fail('Q repeats are not parser-identical')
 A=[c for r,c in configs.items() if r.startswith('A_')]
 aset={(c['kd_corrinoid_btuB'],c['_campaign']['arm'],c['seed']) for c in A}; awant={(k,t,s) for k in C['axes']['A']['kd_corrinoid_btuB_mol_m3'] for t in ('producer','null') for s in C['seeds']}
 if aset!=awant: fail('Stage A explicit run set does not equal 4 Kd x 2 arms x 3 paired seeds')
 for k in C['axes']['A']['kd_corrinoid_btuB_mol_m3']:
  factor=1+1e-3/k
  if not any(approx(factor,w,1e-12) for w in (1001.0,101.0,11.0,2.0)): fail(f'derived competition factor incorrect for {k}: {factor}')
 B=[c for r,c in configs.items() if r.startswith('B_')]
 if len({c['kd_corrinoid_btuB'] for c in B})!=1 or {c['b12_initial_conc'] for c in B}!=set(C['axes']['B']['b12_initial_conc_mol_m3']): fail('Stage B fixed-Kd/supply-only axis violated')
 for c in B:
  if c['initial_strains'][0].get('receptor_expression')!={'BtuB':1.0} or c['initial_strains'][1].get('receptor_expression')!={'BtuB':0.0} or any(x['plasmids'] for x in c['initial_strains']): fail(f"{c['_campaign']['run_id']}: Stage B is not toxin-free BtuB-normal/null")
 CP=[c for r,c in configs.items() if r.startswith('C_amp')]; CN=[c for r,c in configs.items() if r.startswith('C_shared')]
 if len(CP)!=9 or len(CN)!=3 or {c[MUCIN_AMPLITUDE_KEY] for c in CP}!={0,15,60}: fail('Stage C explicit 9 producer + 3 shared-null design violated')
 D=[c for r,c in configs.items() if r.startswith('D_')]
 if len({(c['kd_corrinoid_btuB'],c['b12_initial_conc'],c[MUCIN_AMPLITUDE_KEY]) for c in D})!=1: fail('Stage D fixed upstream axes violated')
 DP=[c for c in D if c['_campaign']['arm']=='producer']; DN=[c for c in D if c['_campaign']['arm']=='plasmid_free_null']
 if len(DP)!=12 or len(DN)!=3: fail('Stage D must contain exactly 12 ColE1 carriers plus 3 plasmid-free nulls')
 if {c['seed'] for c in DN}!=set(C['seeds']): fail('Stage D plasmid-free null set must contain exactly one run per seed')
 if any(c['initial_strains'][0].get('plasmids')!=[] for c in DN): fail('Stage D formal null unexpectedly carries a plasmid')
 if any(c['initial_strains'][0].get('plasmids')!=['ColE1'] for c in DP): fail('Stage D carrier arm lost ColE1')
 if any(c['_campaign']['axes'].get('nominal_target_per_generation') is not None for c in DN): fail('Stage D null controls were accidentally crossed with the lysis axis')
 for c in DP:
  target=c['_campaign']['axes']['nominal_target_per_generation']; want=1-math.sqrt(1-target)
  if not approx(c['sos_lysis_prob'],want,1e-10): fail(f"{c['_campaign']['run_id']}: lysis transform incorrect")
  if c.get('sos_basal_rate')!=0 or c.get(OXYGEN_K_ROS_KEY)!=0: fail(f"{c['_campaign']['run_id']}: lysis attribution controls missing")
 for c in DN:
  if c.get('sos_basal_rate')!=0 or c.get('sos_lysis_prob')!=0 or c.get(OXYGEN_K_ROS_KEY)!=0: fail(f"{c['_campaign']['run_id']}: null attribution controls missing")
 if len(image_digests)!=1: fail(f'mixed image digests: {image_digests}')
 image=next(iter(image_digests),None); pattern=C['execution']['image_digest_pattern']
 if not image or not re.fullmatch(pattern,image):
  msg=f'container image digest is not pinned ({image!r}); cloud submission is blocked'; (fail if a.deployment else WARN.append)(msg)
 if a.deployment and execution_shas!={a.execution_source_sha}: fail(f'mixed/wrong execution source SHAs: {execution_shas}')
 if baseline_shas!={BASELINE}: fail(f'mixed/wrong model baseline SHAs: {baseline_shas}')
 print(json.dumps({'status':'FAIL' if ERR else 'PASS','mode':'deployment' if a.deployment else 'planning','model_baseline_sha':BASELINE,'execution_source_sha':next(iter(execution_shas),None),'git_head':head,'runs_checked':len(signatures),'errors':ERR,'warnings':WARN},indent=2))
 return 1 if ERR else 0
if __name__=='__main__': sys.exit(main())
