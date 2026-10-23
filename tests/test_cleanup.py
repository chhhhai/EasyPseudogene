"""Cleanup regression checks using isolated copies of the synthetic snapshot."""
import contextlib
import fcntl
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'lib'))
import cleanup_intermediates as cleanup
from reference_cache import artifacts, sha


class CleanupTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.run_dir = self.directory / 'run'
        expected = ROOT / 'examples/demo/expected'
        shutil.copytree(expected, self.run_dir)
        self.work = self.run_dir / 'work'
        (self.work / 'inputs').mkdir(parents=True)
        (self.work / 'miniprot').mkdir()
        (self.work / 'miniprot/search.paf').write_text('synthetic temporary alignment\n')
        (self.work / 'genewise').mkdir()
        shutil.copy2(expected / 'genewise_candidate.txt', self.work / 'genewise/demo_A.t1.gw.txt')
        for name in ('.input.lock', '.reference.lock'):
            (self.run_dir / name).touch()
        original = self.directory / 'original.fa'; original.write_text('>synthetic\nATGC\n')
        normalized = self.work / 'inputs/target.genome.fa'
        normalized.write_bytes(original.read_bytes())
        data = dict(signature=dict(sources=dict(genome=dict(path=str(original), sha256=sha(original)))),
                    normalized_paths=dict(genome=str(normalized)))
        (self.run_dir / 'metadata/input_manifest.json').write_text(json.dumps(data))
        (self.work / 'inputs/input_manifest.json').write_text(json.dumps(data))
        manifest = dict(layout_version=4, workflow='standard', comparison_mode='cross', files={},
                        note='Isolated cleanup helper fixture; no alignment analysis executed')
        (self.run_dir / 'metadata/run_manifest.json').write_text(json.dumps(manifest))
        cache = dict(status='complete', signature=dict(inputs={str(normalized): sha(normalized)}))
        self.cache_path = self.run_dir / 'metadata/reference_cache.json'
        self.cache_path.write_text(json.dumps(cache))
        self.refresh()

    def refresh(self):
        cache = json.loads(self.cache_path.read_text())
        cache['artifacts'] = artifacts(self.run_dir)
        self.cache_path.write_text(json.dumps(cache))

    def state(self):
        return {str(p.relative_to(self.run_dir)): ('link:' + os.readlink(p) if p.is_symlink() else sha(p))
                for p in self.run_dir.rglob('*') if p.is_file() or p.is_symlink()}

    def invoke(self, *options):
        return subprocess.run([sys.executable, str(ROOT / 'lib/cleanup_intermediates.py'),
                               '--run-dir', str(self.run_dir), *options],
                              capture_output=True, text=True, timeout=20)

    def assert_fingerprints(self):
        cache = json.loads(self.cache_path.read_text())
        self.assertEqual(cache['status'], 'complete')
        self.assertEqual(cache['artifacts'], artifacts(self.run_dir))

    def test_dry_run_does_not_modify_files(self):
        before = self.state()
        result = self.invoke('--dry-run')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.state(), before)

    def test_cleanup_preserves_published_files_and_final_raw(self):
        before = self.state()
        raw = (self.work / 'genewise/demo_A.t1.gw.txt').read_bytes()
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(list(self.work.iterdir()), [])
        for name, digest in before.items():
            if name == 'report.html' or name.startswith(('results/', 'qc/')):
                self.assertEqual(sha(self.run_dir / name), digest)
        self.assertTrue((self.run_dir / 'results/final_genewise.txt').read_bytes().startswith(raw))
        self.assert_fingerprints()
        after = self.state()
        self.assertEqual(self.invoke().returncode, 0)
        self.assertEqual(self.state(), after)

    def test_both_run_locks_block_cleanup(self):
        for name in ('.input.lock', '.reference.lock'):
            with self.subTest(lock=name):
                before = self.state()
                with (self.run_dir / name).open('r+') as handle:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    self.assertNotEqual(self.invoke().returncode, 0)
                self.assertEqual(self.state(), before)

    def test_incomplete_run_is_not_cleaned(self):
        data = json.loads(self.cache_path.read_text()); data['status'] = 'running'
        self.cache_path.write_text(json.dumps(data))
        before = self.state()
        self.assertNotEqual(self.invoke().returncode, 0)
        self.assertEqual(self.state(), before)

    def test_modified_published_output_is_not_cleaned(self):
        with (self.run_dir / 'results/final_loci.tsv').open('a') as handle: handle.write('\n')
        before = self.state()
        self.assertNotEqual(self.invoke().returncode, 0)
        self.assertEqual(self.state(), before)

    def test_original_input_under_work_is_protected(self):
        original = self.work / 'miniprot/original.fa'; original.write_text('>source\nATGC\n')
        path = self.run_dir / 'metadata/input_manifest.json'
        data = json.loads(path.read_text()); data['signature']['sources']['genome']['path'] = str(original)
        path.write_text(json.dumps(data)); self.refresh()
        before = self.state()
        self.assertNotEqual(self.invoke().returncode, 0)
        self.assertEqual(self.state(), before)

    def test_explicit_reuse_input_in_staging_is_protected(self):
        original = self.work / 'inputs/user_reuse.paf'; original.write_text('user-supplied alignment\n')
        data = json.loads(self.cache_path.read_text()); data['signature']['inputs'][str(original)] = sha(original)
        self.cache_path.write_text(json.dumps(data))
        before = self.state()
        self.assertNotEqual(self.invoke().returncode, 0)
        self.assertEqual(self.state(), before)

    def test_symlink_targets_and_unknown_files_are_retained(self):
        outside = self.directory / 'external'; outside.mkdir()
        original = outside / 'keep.fa'; original.write_text('>keep\nATGC\n')
        (self.work / 'mmseqs').symlink_to(outside, target_is_directory=True)
        (self.work / 'tmp').mkdir()
        (self.work / 'tmp/alias').symlink_to(original)
        keep = self.work / 'user_notes'; keep.write_text('keep\n')
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(original.read_text(), '>keep\nATGC\n')
        self.assertEqual(keep.read_text(), 'keep\n')

    def test_missing_final_alignment_refuses_before_deletion(self):
        (self.work / 'genewise/demo_A.t1.gw.txt').unlink()
        before = self.state()
        self.assertNotEqual(self.invoke().returncode, 0)
        self.assertEqual(self.state(), before)

    def test_work_symlink_is_refused(self):
        outside = self.directory / 'work_copy'
        self.work.rename(outside)
        self.work.symlink_to(outside, target_is_directory=True)
        self.assertNotEqual(self.invoke().returncode, 0)
        self.assertTrue((outside / 'genewise/demo_A.t1.gw.txt').is_file())

    def test_partial_deletion_can_be_retried(self):
        remove = shutil.rmtree

        def fail_miniprot(path, *args, **kwargs):
            if Path(path) == self.work / 'miniprot':
                raise PermissionError('Controlled deletion failure')
            return remove(path, *args, **kwargs)

        with contextlib.redirect_stdout(io.StringIO()), patch.object(cleanup.shutil, 'rmtree', fail_miniprot):
            with self.assertRaises(PermissionError):
                cleanup.cleanup(self.run_dir)
        path = self.run_dir / 'metadata/cleanup.json'
        self.assertEqual(json.loads(path.read_text())['status'], 'partial')
        self.assert_fingerprints()
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(path.read_text())['status'], 'complete')
        self.assert_fingerprints()

    def test_zero_final_candidates_need_no_genewise_archive(self):
        table = self.run_dir / 'results/final_loci.tsv'
        table.write_text(table.read_text().splitlines()[0] + '\n')
        (self.run_dir / 'results/final_loci.bed').write_text('')
        (self.run_dir / 'results/final_loci.gff3').write_text('##gff-version 3\n')
        (self.run_dir / 'results/final_queries.txt').write_text('')
        path = self.run_dir / 'results/summary.json'
        summary = json.loads(path.read_text()); summary['final_loci'] = summary['final_queries'] = 0
        path.write_text(json.dumps(summary)); self.refresh()
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse((self.run_dir / 'results/final_genewise.txt').exists())
        self.assert_fingerprints()

    def test_mount_guard_prevents_deletion(self):
        before = self.state()
        with patch.object(cleanup.os.path, 'ismount', return_value=True):
            with self.assertRaisesRegex(ValueError, 'local work directory'):
                cleanup.plan(self.run_dir, json.loads(self.cache_path.read_text()))
        self.assertEqual(self.state(), before)


if __name__ == '__main__':
    unittest.main()
