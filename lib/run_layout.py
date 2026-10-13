"""Resolve current work/cache locations and support explicit legacy publishers."""
from pathlib import Path


def work_directory(run_dir):
    run_dir = Path(run_dir)
    work = run_dir / 'work'
    if not work.is_dir():
        return run_dir
    stages = ('regions', 'miniprot', 'genewise', 'post_genewise_filter', 'fast_catalog', 'mmseqs_fast', 'blast_fast')
    # A new preflight may coexist with an explicitly opened older analysis.
    # Read that analysis's actual stage files rather than the preflight folder.
    if any((run_dir / name).is_dir() for name in stages) and not any((work / name).is_dir() for name in stages):
        return run_dir
    return work


def cache_file(run_dir):
    run_dir = Path(run_dir)
    if work_directory(run_dir) != run_dir:
        return run_dir / 'metadata/reference_cache.json'
    return run_dir / 'reference_cache.json'
