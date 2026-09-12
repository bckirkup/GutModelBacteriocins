#!/usr/bin/env python3
"""Print (never execute) digest-pinned AWS S3 upload and Batch array commands."""
import argparse,json,re,shlex,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent; REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO/"python"))
from gut_ibm_tools.path_utils import PathValidationError,validate_input_path
ap=argparse.ArgumentParser(); ap.add_argument('--stage',choices=list('QABCD'),required=True); ap.add_argument('--execution-source-sha',required=True,help='full commit used to build image; must match git HEAD and manifests'); ap.add_argument('--image-uri',required=True,help='ECR URI ending @sha256:<64 hex>'); ap.add_argument('--input-prefix',required=True); ap.add_argument('--output-prefix',required=True); ap.add_argument('--job-name'); ap.add_argument('--job-queue',required=True); ap.add_argument('--job-definition',required=True); ap.add_argument('--approval-file',type=Path); a=ap.parse_args()
if not re.fullmatch(r'[0-9a-f]{40}',a.execution_source_sha): raise SystemExit('REFUSED: --execution-source-sha must be full lowercase 40-hex')
try:
 if subprocess.check_output(['git','-C',str(REPO),'rev-parse','--is-inside-work-tree'],text=True,stderr=subprocess.STDOUT).strip()!='true': raise ValueError('not a work tree')
 head=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True,stderr=subprocess.STDOUT).strip()
except Exception as e: raise SystemExit(f'REFUSED: command generation requires a normal git checkout: {e}')
if head!=a.execution_source_sha: raise SystemExit(f'REFUSED: --execution-source-sha differs from git HEAD ({head})')
if not re.fullmatch(r'.+@sha256:[0-9a-f]{64}',a.image_uri): raise SystemExit('REFUSED: --image-uri must be digest pinned')
m=json.loads((ROOT/'generated'/f'stage_{a.stage}'/'manifest.json').read_text())
if m.get('execution_source_sha')!=a.execution_source_sha: raise SystemExit('REFUSED: generated manifest execution_source_sha differs; regenerate deployment artifacts from this HEAD')
if m['container_image_digest']!=a.image_uri.rsplit('@',1)[1]: raise SystemExit('REFUSED: generated manifest digest differs; rerun prepare.py --deployment with this digest')
if a.stage!='Q':
 if not a.approval_file: raise SystemExit(f'REFUSED: Stage {a.stage} requires a gate approval JSON')
 try:
  approval_path=validate_input_path(a.approval_file)
  d=json.loads(approval_path.read_text()); need=m['gate_requires']
 except (OSError,json.JSONDecodeError,PathValidationError) as e:
  raise SystemExit(f'REFUSED: Stage {a.stage} approval file unreadable: {e}')
 if d.get(need) is not True or d.get('model_baseline_sha')!=m['model_baseline_sha'] or d.get('execution_source_sha')!=m['execution_source_sha'] or d.get('image_digest')!=m['container_image_digest']: raise SystemExit(f'REFUSED: approval must set {need}=true and match model baseline, execution source, and image digest')
q=shlex.quote; local=ROOT/'generated'/f'stage_{a.stage}'/'jobs'; n=len(m['runs']); name=a.job_name or f"gutibm-receptor-v1-{a.stage.lower()}"
print('# Review these commands. This program does not execute them.')
print(f'# model_baseline_sha={m["model_baseline_sha"]}')
print(f'# execution_source_sha={m["execution_source_sha"]}')
print(f'# Verify job definition {a.job_definition} resolves to exactly {a.image_uri}.')
print(f"aws s3 cp {q(str(local))} {q(a.input_prefix.rstrip('/')+'/')} --recursive --exclude '*' --include '*/input.json'")
env=[{'name':'INPUT_S3_PREFIX','value':a.input_prefix.rstrip('/')},{'name':'OUTPUT_S3_PREFIX','value':a.output_prefix.rstrip('/')},{'name':'MPI_RANKS','value':'1'},{'name':'REQUIRE_GPU','value':'1'},{'name':'GPU_DEVICE_ID','value':'0'}]
cmd=['aws','batch','submit-job','--job-name',name,'--job-queue',a.job_queue,'--job-definition',a.job_definition,'--array-properties',f'size={n}','--timeout','attemptDurationSeconds=7200','--container-overrides',json.dumps({'environment':env},separators=(',',':'))]
print(' '.join(q(x) for x in cmd))
