"""Checks for demo reproducibility and rejection of incomplete/changed runs."""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEMO = REPO / 'examples' / 'demo'
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')


class DemoTests(unittest.TestCase):
    def test_generator_reproduces_bundled_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp) / 'generate_data.py'
            shutil.copy2(DEMO / 'generate_data.py', script)
            subprocess.run([sys.executable, str(script)], env=ENV, check=True,
                           capture_output=True, text=True)
            for original in (DEMO / 'data').iterdir():
                self.assertEqual(original.read_bytes(), (Path(tmp) / 'data' / original.name).read_bytes())

    def test_manifest_matches_bundled_inputs(self):
        truth = json.loads((DEMO / 'data/truth.json').read_text())
        self.assertEqual(truth['genome_bases'], {'cross': 4205, 'self': 6308})
        for name, expected in truth['files'].items():
            self.assertEqual(hashlib.sha256((DEMO / 'data' / name).read_bytes()).hexdigest(), expected)
        protein_text = (DEMO / 'data/proteins.fa').read_text()
        self.assertEqual(protein_text.count('>'), 2)
        self.assertIn('demo_A.t1', (DEMO / 'data/self_annotation.gff3').read_text())

    def test_missing_output_is_not_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run([sys.executable, str(DEMO / 'check_results.py'),
                                     '--output', str(Path(tmp) / 'absent')], env=ENV,
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('DEMO FAIL', result.stderr)
            self.assertNotIn('DEMO PASS', result.stdout)

    def test_edited_dataset_is_not_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            copied = Path(tmp) / 'demo'
            shutil.copytree(DEMO, copied)
            with (copied / 'data/proteins.fa').open('a') as handle:
                handle.write('\n')
            output = Path(tmp) / 'output'
            output.mkdir()
            result = subprocess.run([sys.executable, str(copied / 'check_results.py'),
                                     '--output', str(output)], env=ENV,
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Demo data changed', result.stderr)
            self.assertEqual(json.loads((output / 'demo_check.json').read_text())['status'], 'failed')

    def test_existing_directory_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'existing'
            output.mkdir()
            (output / 'keep.txt').write_text('preserve me\n')
            result = subprocess.run(['bash', str(DEMO / 'run_demo.sh'), str(output)],
                                    env=ENV, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('already exists', result.stderr)
            self.assertEqual([p.name for p in output.iterdir()], ['keep.txt'])
            self.assertEqual((output / 'keep.txt').read_text(), 'preserve me\n')


if __name__ == '__main__':
    unittest.main()
