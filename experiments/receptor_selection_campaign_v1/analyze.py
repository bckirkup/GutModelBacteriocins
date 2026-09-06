#!/usr/bin/env python3
"""Summarize campaign HDF5 outputs into run, paired, and gate-ready CSV/JSON.
Missing or incomplete output is an explicit BLOCKED result, never an inferred zero.

Paired slopes/tails are recomputed on a shared calendar window ending at
min(t_end_treatment, t_end_control). Own-end windows remain on run rows for
diagnostics only.

realized_lysis_per_division uses a producer-exposure division denominator
(composition-weighted cumulative division increments). The all-type ratio is
retained only as realized_lysis_per_total_division and must not be used as the
Stage D lysis response.
"""
from __future__ import annotations
import argparse,csv,gzip,json,math,shutil,tempfile
from pathlib import Path
try:
 import numpy as np
 NUMPY_IMPORT_ERROR=None
except Exception as exc:
 np=None
 NUMPY_IMPORT_ERROR=f'{type(exc).__name__}: {exc}'
try:
 import h5py
 H5PY_IMPORT_ERROR=None
except Exception as exc:
 h5py=None
 H5PY_IMPORT_ERROR=f'{type(exc).__name__}: {exc}'
ROOT=Path(__file__).resolve().parent; G=ROOT/'generated'
TAIL_S=7200; SENS_S=3600
def val(x):
 a=np.asarray(x[()]); v=a.reshape(-1)[0] if a.size else None
 if isinstance(v,(bytes,np.bytes_)): return v.decode(errors='replace')
 if isinstance(v,np.generic): return v.item()
 return v
def steps(g): return sorted(g,key=lambda x:int(x.rsplit('_',1)[1]))
def open_h5(path):
 if np is None: raise RuntimeError(f'numpy unavailable: {NUMPY_IMPORT_ERROR}')
 if h5py is None: raise RuntimeError(f'h5py unavailable: {H5PY_IMPORT_ERROR}')
 if path.suffix!='.gz': return h5py.File(path,'r'),None
 tmp=tempfile.NamedTemporaryFile(suffix='.h5',delete=False); tmp.close()
 with gzip.open(path,'rb') as a,open(tmp.name,'wb') as b: shutil.copyfileobj(a,b)
 return h5py.File(tmp.name,'r'),Path(tmp.name)
def slope(ts,ys):
 if len(ts)<2 or max(ts)==min(ts): return None
 return float(np.polyfit(np.asarray(ts,dtype=float),np.asarray(ys,dtype=float),1)[0])
def window_metrics(series,t_end,width_s):
 """Slope/tail on series points with t in [t_end-width_s, t_end], clipped to available times."""
 if not series or t_end is None: return {'slope_log10_ratio_per_h':None,'tail_median_log10_ratio':None,'samples':0}
 z=[x for x in series if x[0]<=t_end+1e-12 and x[0]>=t_end-width_s-1e-12]
 s=slope([x[0] for x in z],[x[3] for x in z])
 return {'slope_log10_ratio_per_h':None if s is None else s*3600,'tail_median_log10_ratio':float(np.median([x[3] for x in z])) if z else None,'samples':len(z)}
def get_event(last,names):
 ev=last.get('events')
 if ev is None:return None
 for n in names:
  if n in ev:return val(ev[n])
 return None
def producer_division_exposure(summaries):
 """Allocate cumulative-division increments by live type-1 fraction at interval start.

 HDF5 does not emit per-type divisions. Interval `divisions` is often zero on
 sparse dumps; cumulative deltas between summary samples are the durable signal.
 Type-1 fraction at the start of each increment is the producer-exposure weight.
 """
 if len(summaries)<2: return None
 total=0.0
 for prev,cur in zip(summaries,summaries[1:]):
  ddiv=cur['cumdiv']-prev['cumdiv']
  if ddiv<=0: continue
  n1,n2=prev['n1'],prev['n2']; denom=n1+n2
  if denom<=0: continue
  total+=ddiv*(n1/denom)
 return float(total)
def output_for(stage,index,entry):
 candidates=[G/f'stage_{stage}'/'results'/str(index)/'output.h5.gz',G/f'stage_{stage}'/'results'/str(index)/'output.h5',(ROOT/entry['input_relpath']).parent/'output.h5',(ROOT/entry['input_relpath']).parent/'output.h5.gz']
 return next((p for p in candidates if p.exists()),None)
def one(entry):
 cfg=json.loads((ROOT/entry['input_relpath']).read_text()); m=cfg['_campaign']; p=output_for(m['stage'],entry['array_index'],entry)
 base={'run_id':entry['run_id'],'stage':m['stage'],'arm':m['arm'],'seed':cfg['seed'],'output_status':'missing','output_path':str(p) if p else None,'analysis_error':None,'model_baseline_sha':m['model_baseline_sha'],'execution_source_sha_expected':m['execution_source_sha']}
 for k,v in m['axes'].items(): base['axis_'+k]=v
 if not p:return base
 temp=None
 try:
  h,temp=open_h5(p)
  with h:
   for req in ('run_provenance','summary','agents'):
    if req not in h: raise ValueError(f'missing /{req}')
   rp=h['run_provenance']; source=str(val(rp['git_sha'])) if 'git_sha' in rp else None; placement=str(val(rp['chemistry_placement'])) if 'chemistry_placement' in rp else None
   term=str(val(rp['termination_cause'])) if 'termination_cause' in rp else 'missing'
   ss=steps(h['summary']); aa=steps(h['agents']); series=[]; summary_rows=[]
   for sk in ss:
    s=h['summary'][sk]; t=float(val(s['time'])) if 'time' in s else int(sk.rsplit('_',1)[1])*float(cfg['bio_dt'])
    nbt=np.asarray(s['n_by_type'][()]).ravel() if 'n_by_type' in s else None
    n1=int(nbt[1]) if nbt is not None and nbt.size>1 else None; n2=int(nbt[2]) if nbt is not None and nbt.size>2 else None
    cumdiv=get_event(s,['cumulative_divisions','divisions'])
    if n1 is not None and n2 is not None and cumdiv is not None: summary_rows.append({'t':t,'n1':n1,'n2':n2,'cumdiv':float(cumdiv)})
   # Agent dumps drive the composition time series used for slopes.
   for sk in aa:
    ag=h['agents'][sk]; typ=np.asarray(ag['type'][()]); step=int(sk.rsplit('_',1)[1]); t=step*float(cfg['bio_dt'])
    sumkey=sk if sk in h['summary'] else None
    if sumkey and 'time' in h['summary'][sumkey]: t=float(val(h['summary'][sumkey]['time']))
    n1=int((typ==1).sum());n2=int((typ==2).sum()); series.append((t,n1,n2,math.log10((n1+0.5)/(n2+0.5))))
   if not series: raise ValueError('no agent snapshots')
   tend=series[-1][0]; last=h['summary'][ss[-1]]
   w2=window_metrics(series,tend,TAIL_S); w1=window_metrics(series,tend,SENS_S)
   kills=get_event(last,['cumulative_mortality_colicin','mortality_colicin']); lys=get_event(last,['cumulative_mortality_lysis','mortality_lysis']); div=get_event(last,['cumulative_divisions','divisions']); outflow=get_event(last,['cumulative_outflow_boundary','outflow_boundary','cumulative_boundary_exports'])
   prod_div=producer_division_exposure(summary_rows)
   ag=h['agents'][aa[-1]]; typ=np.asarray(ag['type'][()]); mu=np.asarray(ag['mu_realized'][()] if 'mu_realized' in ag else ag['mu'][()]) if ('mu_realized'in ag or'mu'in ag) else np.array([])
   for t in (1,2): base[f'n_type{t}']=int((typ==t).sum()); base[f'mean_mu_type{t}']=float(np.mean(mu[typ==t])) if mu.size and np.any(typ==t) else None
   btu=[]
   if 'lineage' in h and aa[-1] in h['lineage'] and 'btuB_expression' in h['lineage'][aa[-1]]: btu=np.asarray(h['lineage'][aa[-1]]['btuB_expression'][()])
   # Own-end windows are diagnostics; paired contrasts recompute on t_common below.
   base.update({'output_status':'complete' if term=='horizon_reached' else 'terminated','execution_source_sha_observed':source,'execution_source_sha_match':source==m['execution_source_sha'],'chemistry_placement':placement,'termination_cause':term,'t_end_s':tend,**w2,'sensitivity_1h_slope_log10_ratio_per_h':w1['slope_log10_ratio_per_h'],'sensitivity_1h_tail_median_log10_ratio':w1['tail_median_log10_ratio'],'divisions':div,'producer_divisions_exposure':prod_div,'mortality_colicin':kills,'mortality_lysis':lys,'outflow_boundary':outflow,'kills_per_lysis':kills/lys if kills is not None and lys else None,'kills_per_division':kills/div if kills is not None and div else None,'realized_lysis_per_total_division':lys/div if lys is not None and div else None,'realized_lysis_per_division':lys/prod_div if lys is not None and prod_div else None,'realized_lysis_per_producer_division':lys/prod_div if lys is not None and prod_div else None,'mean_btuB_expression_final':float(np.mean(btu)) if len(btu) else None,'_series':series})
   if m['stage']=='C':
    names=[]
    if 'grid'in h:
     for sk in steps(h['grid']): names.extend(list(h['grid'][sk].keys()))
    base['btuB_grid_present']=bool(names) and set(names)=={'bacteriocin_BtuB'}; base['grid_species_observed']=';'.join(sorted(set(names)))
    base['source_centered_profile_status']='AVAILABLE_FOR_POSTPROCESSING' if base['btuB_grid_present'] and 'kill_provenance' in h else 'BLOCKED_MISSING_SOURCE_EVENT_COORDINATES'
 except Exception as e: base['output_status']='invalid';base['analysis_error']=f'{type(e).__name__}: {e}'
 finally:
  if temp: temp.unlink(missing_ok=True)
 return base
def public_row(r): return {k:v for k,v in r.items() if not k.startswith('_')}
def write_csv(path,rows):
 keys=sorted({k for r in rows for k in r});
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def paired_contrast(r,ctrl,label):
 t_common=min(r['t_end_s'],ctrl['t_end_s'])
 rw=window_metrics(r.get('_series') or [],t_common,TAIL_S); cw=window_metrics(ctrl.get('_series') or [],t_common,TAIL_S)
 return {'stage':r['stage'],'run_id':r['run_id'],'control_run_id':ctrl['run_id'],'seed':r['seed'],'contrast':label,'t_common_s':t_common,'window_s':TAIL_S,'slope_log10_ratio_per_h_treatment':rw['slope_log10_ratio_per_h'],'slope_log10_ratio_per_h_control':cw['slope_log10_ratio_per_h'],'tail_median_log10_ratio_treatment':rw['tail_median_log10_ratio'],'tail_median_log10_ratio_control':cw['tail_median_log10_ratio'],'delta_slope_log10_ratio_per_h':None if rw['slope_log10_ratio_per_h'] is None or cw['slope_log10_ratio_per_h'] is None else rw['slope_log10_ratio_per_h']-cw['slope_log10_ratio_per_h'],'delta_tail_median_log10_ratio':None if rw['tail_median_log10_ratio'] is None or cw['tail_median_log10_ratio'] is None else rw['tail_median_log10_ratio']-cw['tail_median_log10_ratio'],'samples_treatment':rw['samples'],'samples_control':cw['samples']}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--results-dir',type=Path,default=ROOT/'analysis');a=ap.parse_args();a.results_dir.mkdir(parents=True,exist_ok=True)
 cm=json.loads((G/'campaign_manifest.json').read_text()); rows=[one(e) for e in cm['runs']]; pub=[public_row(r) for r in rows]
 write_csv(a.results_dir/'run_metrics.csv',pub);(a.results_dir/'run_metrics.json').write_text(json.dumps(pub,indent=2,allow_nan=False)+'\n')
 complete=[r for r in rows if r['output_status'] in ('complete','terminated')]; missing=[public_row(r) for r in rows if r['output_status'] not in ('complete','terminated')]
 idx={(r['stage'],r['arm'],r['seed'],r.get('axis_kd_corrinoid_btuB_mol_m3'),r.get('axis_amplitude')):r for r in complete}; pairs=[]
 for r in complete:
  ctrl=None; label=None
  if r['stage']=='A' and r['arm']=='producer': ctrl=idx.get(('A','null',r['seed'],r.get('axis_kd_corrinoid_btuB_mol_m3'),None));label='producer_minus_null'
  elif r['stage']=='C' and r['arm']=='producer': ctrl=next((x for x in complete if x['stage']=='C' and x['arm']=='shared_null' and x['seed']==r['seed']),None);label='producer_minus_shared_null'
  elif r['stage']=='D': ctrl=next((x for x in complete if x['stage']=='C' and x['arm']=='shared_null' and x['seed']==r['seed']),None);label='producer_minus_plasmid_free_null'
  if ctrl: pairs.append(paired_contrast(r,ctrl,label))
 write_csv(a.results_dir/'paired_metrics.csv',pairs);(a.results_dir/'paired_metrics.json').write_text(json.dumps(pairs,indent=2)+'\n')
 stages={s:[r for r in rows if r['stage']==s] for s in 'QABCD'}; gates={}
 for s in 'QABCD':
  miss=sum(r['output_status'] not in ('complete','terminated') for r in stages[s]); gates[s]={'status':'BLOCKED_MISSING_OUTPUTS' if miss else 'READY_FOR_SCIENTIFIC_REVIEW','expected':len(stages[s]),'missing_or_invalid':miss}
 q=stages['Q']
 if not any(r['output_status'] not in ('complete','terminated') for r in q):
  if any(not r.get('execution_source_sha_match') or r.get('chemistry_placement')!='device_delivery' for r in q): gates['Q']['status']='FAIL_PROVENANCE_OR_PLACEMENT'
 cgood=[r for r in stages['C'] if r['output_status'] in ('complete','terminated')]
 if len(cgood)==12 and any(not r.get('btuB_grid_present') or r.get('source_centered_profile_status','').startswith('BLOCKED') for r in cgood): gates['C']['status']='BLOCKED_TRANSPORT_PROVENANCE'
 report={'campaign_id':cm['campaign_id'],'runs_total':len(rows),'outputs_readable':len(complete),'missing_or_invalid':len(missing),'gates':gates,'interpretation':'Sequential scientific gates require review of all seed-level rows; this analyzer never auto-promotes downstream stages. Paired deltas use a shared calendar window ending at t_common_s. realized_lysis_per_division is producer-exposure-denominator only.'}
 (a.results_dir/'missing_outputs.json').write_text(json.dumps(missing,indent=2)+'\n');(a.results_dir/'gate_status.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));return 0
if __name__=='__main__': raise SystemExit(main())
