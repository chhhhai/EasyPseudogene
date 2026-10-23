"""Conservative reference cache guard: incompatible/incomplete results require a new output directory.

Never deletes previous results. Content hashes cover inputs, configuration, code,
executables and published intermediates. A process lock is held by the CLI.
"""
import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
from run_layout import cache_file, work_directory


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def artifacts(directory):
    result={}
    work = work_directory(directory)
    roots = [work / sub for sub in ['miniprot','genewise','regions','prefilter','post_genewise_filter',
                                    'fast_catalog','disruption_validation','mmseqs_fast','blast_fast']]
    roots += [directory / sub for sub in ['results','qc','metadata']]
    if work != directory:
        roots += [directory / sub for sub in ['miniprot','genewise','regions','prefilter','post_genewise_filter','fast_catalog']]
    cache = cache_file(directory)
    for root in roots:
        if root.exists():
            for p in sorted(root.rglob('*')):
                if p.is_file() and p != cache and 'tmp' not in p.relative_to(root).parts:
                    result[str(p.relative_to(directory))]=sha(p)
    report = directory / 'report.html'
    if report.is_file():result['report.html']=sha(report)
    return result


def begin(a):
    if (a.directory / 'metadata/cleanup.json').exists():
        raise ValueError('Intermediates were cleaned or cleanup was interrupted. Published results remain available; use a new output directory for another analysis')
    inputs={str(p.resolve()):sha(p) for p in a.input if p.is_file()}
    source={str(p.relative_to(a.source)):sha(p) for sub in ['bin','scripts','lib','config'] for p in sorted((a.source/sub).rglob('*')) if p.is_file() and p.suffix not in ['.pyc'] and '__pycache__' not in p.parts}
    version_file=a.source/'VERSION'
    if version_file.is_file():source['VERSION']=sha(version_file)
    tools={name:sha(shutil.which(name)) if shutil.which(name) else None for name in ['python3','miniprot','genewise','samtools','mmseqs','gffread','tblastn','makeblastdb']}
    # Retention changes do not change the biological analysis/cache identity.
    analysis_arguments=[arg for arg in a.arguments if arg != '--cleanup-intermediates']
    signature={'inputs':inputs,'source':source,'tools':tools,'arguments':analysis_arguments,'environment':{k:os.environ.get(k,'') for k in ['WISECONFIGDIR','WISEPATH','GENEWISE_BIN','FLANK','MIN_ALN_LEN','MIN_MAPQ','MMSEQS_SEARCH_TYPE','MMSEQS_SENS','MMSEQS_EVALUE','MMSEQS_MAX_SEQS']}}
    wise=Path(os.environ.get('WISECONFIGDIR','/nonexistent'))
    if wise.is_dir():signature['wise_config']={p.name:sha(p) for p in sorted(wise.iterdir()) if p.is_file()}
    path=cache_file(a.directory)
    if path != a.directory / 'reference_cache.json' and (a.directory / 'reference_cache.json').exists():
        raise ValueError('Legacy output layout exists; use a new --output directory')
    if path.exists():
        old=json.loads(path.read_text())
        if old.get('signature')!=signature:raise ValueError('Cache inputs, configuration, code or tools changed; use a new --output directory')
        if old.get('status')!='complete':raise ValueError('Previous reference run is incomplete; use a new --output directory')
        if old.get('artifacts')!=artifacts(a.directory):raise ValueError('Cached results are missing or modified; use a new --output directory')
        if getattr(a, 'reuse_complete', False):
            return True
    elif artifacts(a.directory):
        raise ValueError('Unverified legacy results exist; use a new --output directory')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({'status':'running','signature':signature},indent=2)+'\n')
    return False


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['begin','complete']);p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--source',type=Path);p.add_argument('--input',type=Path,action='append',default=[])
    # Parse invocation JSON as an ordinary option: no shell interpretation.
    p.add_argument('--invocation',default='[]')
    p.add_argument('--reuse-complete', action='store_true', help='Keep a verified complete run unchanged and report reuse')
    a=p.parse_args();a.arguments=json.loads(a.invocation)
    if a.action=='begin':
        reused = begin(a)
        if a.reuse_complete:print('reused' if reused else 'new')
    else:
        path=cache_file(a.directory);data=json.loads(path.read_text());data.update(status='complete',artifacts=artifacts(a.directory))
        tmp=path.with_suffix('.json.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)

if __name__=='__main__':
    try:main()
    except (ValueError,OSError) as e:raise SystemExit('ERROR: '+str(e))
