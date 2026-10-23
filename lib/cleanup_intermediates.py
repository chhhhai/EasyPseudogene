"""Remove known work directories only after a verified, completed analysis."""
import argparse
import csv
import fcntl
import json
import os
import re
import shutil
import stat
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path

from reference_cache import artifacts, sha


WORK_STAGES = ('inputs', 'input_check', 'prefilter', 'miniprot', 'mmseqs', 'tmp',
               'regions', 'genewise', 'post_genewise_filter', 'fast_catalog',
               'disruption_validation', 'mmseqs_fast', 'blast_fast')
REQUIRED = ('report.html', 'results/final_loci.tsv', 'results/final_loci.bed',
            'results/final_loci.gff3', 'results/final_queries.txt', 'results/loci.tsv',
            'results/queries.tsv', 'results/summary.json', 'metadata/run_manifest.json',
            'metadata/input_manifest.json', 'qc/input_qc.json')


def write_json(path, data):
    temporary = path.with_name(path.name + '.partial')
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    temporary.replace(path)


def inside(path, directory):
    return path == directory or directory in path.parents


def published_hashes(run):
    result = {'report.html': sha(run / 'report.html')}
    for name in ('results', 'metadata', 'qc'):
        for path in sorted((run / name).rglob('*')):
            if path.is_symlink():
                raise ValueError('Published files must be local regular files: ' + str(path))
            if path.is_file() and path != run / 'metadata/reference_cache.json':
                result[str(path.relative_to(run))] = sha(path)
    return result


def verify(run, cache):
    if cache.get('status') != 'complete':
        raise ValueError('Analysis is not complete; intermediates were retained')
    if any((run / name).exists() for name in WORK_STAGES):
        raise ValueError('Legacy work directories exist outside work/; automatic cleanup is not supported')
    for name in REQUIRED:
        if not (run / name).is_file():
            raise ValueError('Missing published output: ' + name)
    current = published_hashes(run)
    expected = {key: value for key, value in cache.get('artifacts', {}).items()
                if key == 'report.html' or key.startswith(('results/', 'metadata/', 'qc/'))}
    if not expected or current != expected:
        raise ValueError('Published outputs differ from the completed cache; nothing was deleted')
    manifest = json.loads((run / 'metadata/run_manifest.json').read_text())
    if manifest.get('snapshot_only') or manifest.get('layout_version') != 4:
        raise ValueError('Cleanup requires a complete current-layout run, not a snapshot or legacy layout')
    with (run / 'results/final_loci.tsv').open(newline='') as handle:
        finals = list(csv.DictReader(handle, delimiter='\t'))
    summary = json.loads((run / 'results/summary.json').read_text())
    if summary.get('final_loci') != len(finals):
        raise ValueError('Final candidate table and summary counts differ')
    return manifest, finals


def inventory(path):
    """Count file content without following symlinks or crossing mount points."""
    if path.is_symlink():
        return 1, 0
    if not path.is_dir() or os.path.ismount(path):
        raise ValueError('Expected a local work directory: ' + str(path))
    count = total = 0
    for parent, directories, files in os.walk(path, followlinks=False):
        for name in directories:
            child = Path(parent) / name
            if not child.is_symlink() and os.path.ismount(child):
                raise ValueError('Work directory contains a mount point: ' + str(child))
        for name in files:
            info = (Path(parent) / name).lstat()
            count += 1
            if stat.S_ISREG(info.st_mode):
                total += info.st_size
    return count, total


def plan(run, cache):
    work = run / 'work'
    if not work.is_dir() or work.is_symlink():
        raise ValueError('work/ must be a local directory')
    targets = [work / name for name in WORK_STAGES if (work / name).exists() or (work / name).is_symlink()]
    inputs = json.loads((run / 'metadata/input_manifest.json').read_text())
    protected = [Path(value['path']).resolve() for value in inputs.get('signature', {}).get('sources', {}).values()]
    generated_inputs = {Path(value).resolve() for value in inputs.get('normalized_paths', {}).values()
                        if value and inside(Path(value).resolve(), work / 'inputs')}
    generated_inputs.add((work / 'inputs/input_manifest.json').resolve())
    # Staged inputs are generated; other explicit cached inputs can be user-supplied
    # reuse files or configuration and must never be deleted.
    for name in cache.get('signature', {}).get('inputs', {}):
        path = Path(name).resolve()
        if path not in generated_inputs:
            protected.append(path)
    for target in targets:
        if target.is_symlink():
            continue  # Unlink the link only; never remove its target.
        if any(inside(path, target) for path in protected):
            raise ValueError('A source input is inside a cleanup directory; retained: ' + str(target))
    return [dict(path=str(target.relative_to(run)), files=(size := inventory(target))[0], logical_bytes=size[1])
            for target in targets]


def preserve_final_genewise(run, manifest, finals, previous):
    """Keep original final-locus GeneWise text even for older summary reports."""
    if manifest.get('workflow') != 'standard' or not finals:
        return None
    destination = run / 'results/final_genewise.txt'
    if destination.is_symlink():
        raise ValueError('final_genewise.txt must not be a symlink')
    if destination.exists():
        if previous.get('final_genewise_sha256') != sha(destination):
            raise ValueError('An unrecognized final_genewise.txt already exists; retained all work files')
        return sha(destination)
    temporary = destination.with_name(destination.name + '.partial')
    if temporary.exists() or temporary.is_symlink():
        raise ValueError('A previous final_genewise.txt.partial needs review before cleanup')
    wanted = {row['locus_id'] for row in finals}
    if any(Path(locus).name != locus or locus in ('.', '..') for locus in wanted):
        raise ValueError('Invalid final locus ID')
    found = set()
    try:
        with temporary.open('xb') as out:
            for locus in sorted(wanted):
                source = run / 'work/genewise' / (locus + '.gw.txt')
                if not source.is_file():
                    continue
                with source.open('rb') as handle:
                    header = b''.join(handle.readline() for _ in range(32))
                    match = re.search(rb'Query protein:\s+(\S+)', header)
                    if not match or match.group(1).decode() != locus:
                        raise ValueError('GeneWise record does not match final locus ' + locus)
                    out.write(header)
                    shutil.copyfileobj(handle, out)
                    out.write(b'\n')
                found.add(locus)
            combined = run / 'work/genewise/pseudogene-genewise.txt'
            if wanted != found and combined.is_file():
                with combined.open('rb') as handle:
                    block = []
                    for line in handle:
                        block.append(line)
                        if line.strip() == b'//':
                            raw = b''.join(block)
                            match = re.search(rb'Query protein:\s+(\S+)', raw)
                            locus = match.group(1).decode() if match else None
                            if locus in wanted - found:
                                out.write(raw + b'\n')
                                found.add(locus)
                            block = []
            if wanted != found:
                raise ValueError('Cannot preserve all final GeneWise records; retained all work files')
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return sha(destination)


def cleanup(run, dry_run=False):
    cache_path = run / 'metadata/reference_cache.json'
    cache = json.loads(cache_path.read_text())
    manifest, finals = verify(run, cache)
    record_path = run / 'metadata/cleanup.json'
    previous = json.loads(record_path.read_text()) if record_path.is_file() else {}
    if previous.get('status') == 'complete':
        print('Intermediate cleanup already completed. Report and results are retained.')
        return
    items = plan(run, cache)
    print(json.dumps(dict(action='preview' if dry_run else 'cleanup', directories=items,
                          logical_bytes=sum(item['logical_bytes'] for item in items),
                          preserved=['report.html', 'results/', 'qc/', 'metadata/', 'logs/'],
                          retained_other_work=[p.name for p in (run / 'work').iterdir() if p.name not in WORK_STAGES]), indent=2))
    if dry_run:
        print('Preview only: no files were written or deleted.')
        return
    archive = preserve_final_genewise(run, manifest, finals, previous)
    record = dict(status='running', started_utc=datetime.now(timezone.utc).isoformat(),
                  directories=items, removed_directories=[], final_genewise_sha256=archive,
                  preserved=['report.html', 'results/', 'qc/', 'metadata/', 'logs/'],
                  retained_other_work=[p.name for p in (run / 'work').iterdir() if p.name not in WORK_STAGES],
                  note='File byte counts are logical sizes, not a guaranteed measure of freed disk space.')
    record['previous_attempts'] = previous.get('previous_attempts', [])
    if previous:
        record['previous_attempts'].append({key: previous.get(key) for key in
                                            ('status', 'started_utc', 'finished_utc', 'removed_directories', 'error')})
    manifest['retention'] = dict(intermediates='cleanup_running', cleanup_record='cleanup.json')
    manifest['files']['cleanup'] = 'cleanup.json'
    if archive:
        manifest['files']['final_genewise'] = '../results/final_genewise.txt'
    write_json(record_path, record)
    write_json(run / 'metadata/run_manifest.json', manifest)
    # Establish a checkpoint for safe retry after interruption; preserve existing
    # work hashes here rather than re-reading hundreds of GB before deletion.
    cache['artifacts'].update(published_hashes(run))
    cache['retention'] = dict(intermediates='cleanup_running')
    write_json(cache_path, cache)
    try:
        work_info = (run / 'work').stat()
        for item in items:
            path = run / item['path']
            current = (run / 'work').stat()
            if (run / 'work').is_symlink() or (current.st_dev, current.st_ino) != (work_info.st_dev, work_info.st_ino):
                raise ValueError('work/ changed during cleanup')
            if path.is_symlink():
                path.unlink()
            else:
                shutil.rmtree(path)
            record['removed_directories'].append(item['path'])
        record['status'] = 'complete'
    except BaseException as exc:
        record.update(status='partial', error=str(exc))
        raise
    finally:
        record['finished_utc'] = datetime.now(timezone.utc).isoformat()
        manifest['retention']['intermediates'] = 'removed' if record['status'] == 'complete' else 'cleanup_partial'
        write_json(record_path, record)
        write_json(run / 'metadata/run_manifest.json', manifest)
        cache['retention'] = manifest['retention']
        cache['artifacts'] = artifacts(run)
        write_json(cache_path, cache)
    print('Intermediate cleanup complete. report.html, results/, qc/, metadata/ and logs/ were retained.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True, help='Completed OUTPUT/NAME directory')
    parser.add_argument('--dry-run', action='store_true', help='Show cleanup scope without writing or deleting files')
    parser.add_argument('--lock-fds', nargs=2, type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    try:
        run = args.run_dir.resolve(strict=True)
        with ExitStack() as stack:
            for index, name in enumerate(('.input.lock', '.reference.lock')):
                path = run / name
                if not path.is_file() or path.is_symlink():
                    raise ValueError('Expected an existing analysis lock file: ' + name)
                if args.lock_fds:
                    fd = args.lock_fds[index]
                    if (os.fstat(fd).st_dev, os.fstat(fd).st_ino) != (path.stat().st_dev, path.stat().st_ino):
                        raise ValueError('Inherited analysis lock does not match this run')
                else:
                    fd = stack.enter_context(path.open('r+')).fileno()
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            for name in ('metadata', 'results', 'qc'):
                if (run / name).is_symlink():
                    raise ValueError('Published output directories must not be symlinks')
            if (run / 'metadata/reference_cache.json').is_symlink():
                raise ValueError('The completed cache must not be a symlink')
            cleanup(run, args.dry_run)
    except (ValueError, OSError, KeyError) as exc:
        raise SystemExit('ERROR: ' + str(exc)) from exc


if __name__ == '__main__':
    main()
