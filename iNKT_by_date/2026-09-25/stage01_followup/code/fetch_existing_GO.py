"""Read-only fetch of the existing GO package after user restores SSH auth."""
from pathlib import Path
import subprocess,shlex,tarfile,json,hashlib
BASE=Path(__file__).resolve().parents[2];DEST=BASE/'work/followup0925/remote_existing_GO';DEST.mkdir(exist_ok=True)
code="""import sys,tarfile
from pathlib import Path
r=Path('/home/zzz0054/bio3/iNKT_by_date/2026-09-20/package')
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as t:
 for p in ['data/current_analysis_tables','notes/GO_provenance.json','notes/input_manifest.json','notes/download_manifest.json']:
  t.add(r/p,arcname=p)
"""
args=['ssh','-S','/Users/zeruzhang/.ssh/codex-zhou2.sock','-o','BatchMode=yes','-o','ForwardX11=no','-o','ForwardAgent=no','zhou2_vpn','python3 -c '+shlex.quote(code)]
a=DEST/'existing_GO.tgz'
with a.open('xb') as f:subprocess.run(args,stdout=f,check=True,timeout=180)
with tarfile.open(a) as t:
 for m in t.getmembers():
  if m.issym() or m.islnk() or Path(m.name).is_absolute() or '..' in Path(m.name).parts: raise RuntimeError('Unsafe archive member')
 t.extractall(DEST,filter='data')
expected=json.loads((BASE/'work/remote37/file_manifest.json').read_text());checks=[]
for p in (DEST/'data/current_analysis_tables').iterdir():
 if p.is_file():
  k=str(p.relative_to(DEST));h=hashlib.sha256(p.read_bytes()).hexdigest();checks.append({'file':k,'sha256':h,'expected':expected.get(k),'matches_previous_package':h==expected.get(k)})
(DEST/'hash_checks.json').write_text(json.dumps(checks,indent=2));print('fetched',len(checks),'files; changed',sum(not x['matches_previous_package'] for x in checks))
if not all(x['matches_previous_package'] for x in checks):raise RuntimeError('Changed source snapshot: review before reuse')
