#!/usr/bin/env python3
"""Run resumable tblastn independently per target FASTA record and combine BLAST6."""
from __future__ import annotations
import argparse,concurrent.futures,hashlib,json,os,shutil,signal,subprocess,sys,threading,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path:sys.path.insert(0,str(HERE))
from catalog_utils import atomic_text_writer,file_sha256,read_fasta,write_fasta
from self_catalog_mmseqs_fast import hardmask_cds

def safe_name(name):
 digest=hashlib.sha256(name.encode()).hexdigest()[:8];clean=''.join(c if c.isalnum() or c in '._-' else '_' for c in name)[:80];return f'{clean}.{digest}'
def run_checked(cmd,stdout=None,stderr=None):
 process=subprocess.run(cmd,stdout=stdout,stderr=stderr,text=True)
 if process.returncode!=0:raise RuntimeError(f'command failed ({process.returncode}): {" ".join(cmd)}')
def fingerprint_text(text):return hashlib.sha256(text.encode()).hexdigest()
def tool_version(command):
 try:
  process=subprocess.run([command,'-version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=20);return process.stdout.splitlines()[0].strip()
 except Exception:return 'unknown'
def main(argv=None):
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--proteins',required=True,type=Path);ap.add_argument('--genome',required=True,type=Path);ap.add_argument('--workdir',required=True,type=Path);ap.add_argument('--out-blast6',required=True,type=Path);ap.add_argument('--summary-json',required=True,type=Path);ap.add_argument('--jobs',type=int,default=4);ap.add_argument('--threads-per-job',type=int,default=4);ap.add_argument('--task',choices=['tblastn','tblastn-fast'],default='tblastn');ap.add_argument('--evalue',type=float,default=0.1);ap.add_argument('--max-target-seqs',type=int,default=0);ap.add_argument('--gff',type=Path);ap.add_argument('--mask-cds',action='store_true');ap.add_argument('--resume',action='store_true');ap.add_argument('--retries',type=int,default=1);ap.add_argument('--min-free-gb',type=float,default=1.0);args=ap.parse_args(argv)
 if args.jobs<1 or args.threads_per_job<1 or args.retries<0:ap.error('jobs/threads must be >=1 and retries >=0')
 args.workdir.mkdir(parents=True,exist_ok=True);free=shutil.disk_usage(args.workdir).free;required=max(int(args.min_free_gb*1024**3),args.genome.stat().st_size*2)
 if free<required:raise SystemExit(f'insufficient free disk: {free/1024**3:.2f} GB available, {required/1024**3:.2f} GB required')
 target_genome=args.genome;mask_summary=None
 if args.mask_cds:
  if not args.gff:ap.error('--mask-cds requires --gff')
  target_genome=args.workdir/'target.CDS_hardmasked.fa';mask_summary=hardmask_cds(args.genome,args.gff,target_genome)
 versions={'tblastn':tool_version('tblastn'),'makeblastdb':tool_version('makeblastdb')};protein_sha=file_sha256(args.proteins);records=list(read_fasta(target_genome));
 if not records:raise SystemExit(f'no FASTA records found in {target_genome}')
 names=[record[0] for record in records]
 if len(names)!=len(set(names)):raise SystemExit('duplicate target FASTA sequence identifiers are not allowed')
 targets=[]
 for index,(name,desc,seq) in enumerate(records):
  tag=safe_name(name);fa=args.workdir/f'{tag}.fa';db=args.workdir/f'{tag}.db';out=args.workdir/f'{tag}.blast6';log=args.workdir/f'{tag}.log';write_fasta([(name,desc,seq)],fa);target_fp=fingerprint_text(json.dumps({'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'protein_sha256':protein_sha,'task':args.task,'evalue':args.evalue,'max_target_seqs':args.max_target_seqs,'tool_versions':versions},sort_keys=True));targets.append((index,name,fa,db,out,log,target_fp))
 manifest_path=args.workdir/'resume_manifest.json';manifest={'schema_version':1,'targets':{}}
 if args.resume and manifest_path.exists():
  try:manifest=json.loads(manifest_path.read_text())
  except Exception:manifest={'schema_version':1,'targets':{}}
 lock=threading.Lock();temporary_files=set()
 def save_manifest():
  with atomic_text_writer(manifest_path) as handle:handle.write(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
 def cleanup(_signum=None,_frame=None):
  for path in list(temporary_files):
   try:path.unlink()
   except FileNotFoundError:pass
  if _signum is not None:raise SystemExit(128+_signum)
 signal.signal(signal.SIGTERM,cleanup);signal.signal(signal.SIGINT,cleanup)
 def worker(item):
  index,name,fa,db,out,log,target_fp=item;t0=time.time();previous=manifest.get('targets',{}).get(name,{})
  if args.resume and previous.get('fingerprint')==target_fp and out.exists() and previous.get('output_sha256')==file_sha256(out):return {'index':index,'target':name,'seconds':0.0,'rows':sum(1 for _ in out.open()),'resumed':True,'attempts':0}
  last_error=None
  for attempt in range(1,args.retries+2):
   tmp=out.with_name(out.name+f'.tmp.{os.getpid()}.{threading.get_ident()}');temporary_files.add(tmp)
   try:
    with log.open('a' if attempt>1 else 'w') as err:
     run_checked(['makeblastdb','-in',str(fa),'-dbtype','nucl','-parse_seqids','-out',str(db)],stdout=err,stderr=err);cmd=['tblastn','-task',args.task,'-query',str(args.proteins),'-db',str(db),'-evalue',str(args.evalue),'-num_threads',str(args.threads_per_job),'-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen','-out',str(tmp)]
     if args.max_target_seqs:cmd+=['-max_target_seqs',str(args.max_target_seqs)]
     run_checked(cmd,stdout=err,stderr=err)
    os.replace(tmp,out);temporary_files.discard(tmp);result={'index':index,'target':name,'seconds':time.time()-t0,'rows':sum(1 for _ in out.open()),'resumed':False,'attempts':attempt}
    with lock:manifest.setdefault('targets',{})[name]={'fingerprint':target_fp,'output_sha256':file_sha256(out),'rows':result['rows']};save_manifest()
    return result
   except Exception as error:
    last_error=error
    try:tmp.unlink()
    except FileNotFoundError:pass
    temporary_files.discard(tmp)
  raise RuntimeError(f'{name} failed after {args.retries+1} attempts: {last_error}')
 t0=time.time();results=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
  futures=[executor.submit(worker,item) for item in targets]
  for future in concurrent.futures.as_completed(futures):results.append(future.result())
 results.sort(key=lambda x:x['index'])
 with atomic_text_writer(args.out_blast6) as combined:
  for item in sorted(targets):
   with item[4].open() as source:
    for line in source:combined.write(line)
 summary={'schema_version':2,'targets':len(targets),'jobs':args.jobs,'threads_per_job':args.threads_per_job,'task':args.task,'evalue':args.evalue,'elapsed_seconds':time.time()-t0,'combined_rows':sum(x['rows'] for x in results),'resumed_targets':sum(x['resumed'] for x in results),'per_target':results,'out_blast6':str(args.out_blast6),'out_blast6_sha256':file_sha256(args.out_blast6),'protein_sha256':protein_sha,'tool_versions':versions,'cds_hardmask':mask_summary,'disk_preflight_required_bytes':required}
 with atomic_text_writer(args.summary_json) as out:out.write(json.dumps(summary,indent=2,sort_keys=True)+'\n')
 save_manifest();print(json.dumps(summary,indent=2,sort_keys=True))
if __name__=='__main__':main()
