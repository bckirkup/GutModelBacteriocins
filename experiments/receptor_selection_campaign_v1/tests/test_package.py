import json, math, subprocess, sys, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASELINE='ad037ef96f3433a1ad4008994956f89c68a3ed39'
EXEC_PLACEHOLDER='EXECUTION_SOURCE_SHA_REQUIRED_AFTER_CAMPAIGN_COMMIT'
IMAGE_PLACEHOLDER='REQUIRED_BEFORE_SUBMISSION'
class PackageTests(unittest.TestCase):
 def test_contract_separates_design_and_execution_identity(self):
  c=json.loads((ROOT/'campaign_contract.json').read_text())
  self.assertEqual(c['schema_version'],2); self.assertEqual(c['model_baseline_sha'],BASELINE)
  self.assertEqual(c['execution_source_sha_policy']['planning_placeholder'],EXEC_PLACEHOLDER)
  self.assertNotIn('source_sha',c)
  self.assertEqual(sum(v['jobs'] for v in c['stages'].values()),62)
  for target,p in zip(c['axes']['D']['realized_lysis_target_per_generation'],c['axes']['D']['sos_lysis_prob']): self.assertTrue(math.isclose(p,1-math.sqrt(1-target),rel_tol=1e-10,abs_tol=1e-12))
 def test_planning_artifacts_record_both_identities(self):
  c=json.loads((ROOT/'campaign_contract.json').read_text())
  cm=json.loads((ROOT/'generated'/'campaign_manifest.json').read_text())
  self.assertEqual(cm['model_baseline_sha'],BASELINE); self.assertEqual(cm['execution_source_sha'],EXEC_PLACEHOLDER); self.assertEqual(cm['container_image_digest'],IMAGE_PLACEHOLDER)
  for stage,meta in c['stages'].items():
   m=json.loads((ROOT/'generated'/f'stage_{stage}'/'manifest.json').read_text())
   self.assertEqual(len(m['runs']),meta['jobs']); self.assertEqual([x['array_index'] for x in m['runs']],list(range(meta['jobs'])))
   self.assertEqual(m['model_baseline_sha'],BASELINE); self.assertEqual(m['execution_source_sha'],EXEC_PLACEHOLDER)
   for entry in m['runs']:
    cfg=json.loads((ROOT/entry['input_relpath']).read_text()); provenance=cfg['_campaign']
    self.assertEqual(provenance['model_baseline_sha'],BASELINE); self.assertEqual(provenance['execution_source_sha'],EXEC_PLACEHOLDER); self.assertNotIn('source_sha',provenance)
 def test_planning_preflight_passes_with_explicit_digest_warning(self):
  p=subprocess.run([sys.executable,str(ROOT/'preflight.py')],capture_output=True,text=True)
  self.assertEqual(p.returncode,0,p.stdout+p.stderr); report=json.loads(p.stdout); self.assertEqual(report['mode'],'planning'); self.assertTrue(any('cloud submission is blocked' in x for x in report['warnings']))
 def test_deployment_preflight_refuses_tarball_and_missing_execution_sha(self):
  p=subprocess.run([sys.executable,str(ROOT/'preflight.py'),'--deployment'],capture_output=True,text=True)
  self.assertEqual(p.returncode,1); report=json.loads(p.stdout); self.assertTrue(any('--execution-source-sha' in x for x in report['errors'])); self.assertTrue(any('git checkout' in x or 'git HEAD' in x for x in report['errors']))
 def test_deployment_generation_refuses_without_git_even_with_pinned_values(self):
  p=subprocess.run([sys.executable,str(ROOT/'prepare.py'),'--deployment','--execution-source-sha','0'*40,'--image-digest','sha256:'+'1'*64],capture_output=True,text=True)
  self.assertNotEqual(p.returncode,0); self.assertIn('normal git checkout',p.stdout+p.stderr)
 def test_runtime_analysis_uses_exact_execution_sha(self):
  text=(ROOT/'analyze.py').read_text(); self.assertIn("source==m['execution_source_sha']",text); self.assertNotIn("startswith(source)",text)
 def test_paired_windows_share_calendar_end(self):
  sys.path.insert(0,str(ROOT)); import analyze
  # Treatment continues past the null; own-end slopes would disagree with a shared end.
  null=[(t,10,10,0.0) for t in range(0,16801,600)]
  prod=[(t,10+t//600,max(1,10-t//1200),math.log10(((10+t//600)+0.5)/(max(1,10-t//1200)+0.5))) for t in range(0,20401,600)]
  t_common=16800.0
  cw=analyze.window_metrics(null,t_common,analyze.TAIL_S); tw=analyze.window_metrics(prod,t_common,analyze.TAIL_S)
  own=analyze.window_metrics(prod,20400.0,analyze.TAIL_S)
  self.assertEqual(cw['samples'],tw['samples']); self.assertGreater(cw['samples'],1)
  self.assertNotEqual(own['slope_log10_ratio_per_h'],tw['slope_log10_ratio_per_h'])
  self.assertTrue(all(x[0]<=t_common for x in null if x[0]>=t_common-analyze.TAIL_S))
 def test_producer_divisions_use_integer_by_type_counter(self):
  sys.path.insert(0,str(ROOT)); import analyze
  class FakeEv(dict):
   def get(self,k,default=None): return super().get(k,default)
  class Arr:
   def __init__(self,a): self._a=a
   def __getitem__(self,_): return self._a
  last={'events':FakeEv({'cumulative_divisions_by_type':Arr([0,17,9,0,0,0,0,0])})}
  self.assertEqual(analyze.producer_divisions_from_events(last),17)
  self.assertIsNone(analyze.producer_divisions_from_events({'events':{}}))
  text=(ROOT/'analyze.py').read_text()
  self.assertIn('cumulative_divisions_by_type',text)
  self.assertIn('realized_lysis_per_total_division',text)
  self.assertNotIn('producer_division_exposure',text)
  self.assertIn('t_common',text)
if __name__=='__main__': unittest.main(verbosity=2)
