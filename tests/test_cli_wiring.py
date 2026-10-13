#!/usr/bin/env python3
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class CliWiringTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.source=(ROOT/'bin'/'easypseudogene').read_text()
 def test_serial_blast_mask_routes_through_mask_capable_runner(self):
  self.assertIn('[[ "${SELF_CATALOG_BLAST_PARALLEL_CHROMOSOMES}" == "1" || "${SELF_CATALOG_BLAST_MASK_CDS}" == "1" ]]',self.source)
  self.assertIn('runner_jobs=1',self.source);self.assertIn('parallel_cmd+=(--mask-cds)',self.source)
 def test_blast_fast_validation_is_wired(self):
  block=self.source.split('# Integrated BLAST6/tblastn self-catalog fast engine.',1)[1].split('# Step 1:',1)[0]
  self.assertIn('run_fast_catalog_validation "${fast_dir}"',block)
  self.assertIn('[[ "${SELF_CATALOG_VALIDATE_DISRUPTIONS}" == "1" ]] || return 0',self.source)
  self.assertIn('validate_candidate_disruptions.py',self.source)
  self.assertIn('.fast_catalog_candidates.validated.tsv',self.source)
  self.assertIn('catalog_tiers.py',self.source)
 def test_tblastn_validated_and_candidates_only_presets_are_distinct(self):
  self.assertIn('tblastn-production|tblastn-production-candidates-only|tblastn-production-validated)',self.source)
  self.assertIn('[[ "${PRESET}" == "tblastn-production-validated" ]] && SELF_CATALOG_VALIDATE_DISRUPTIONS="1"',self.source)
if __name__=='__main__':unittest.main()
