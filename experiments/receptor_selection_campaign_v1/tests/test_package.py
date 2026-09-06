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
if __name__=='__main__': unittest.main(verbosity=2)
