from pathlib import Path
import hashlib,json,os
r=Path('/home/zzz0054/bio3')
roots=['input/iNKT','notebooks','docs','output/iNKT_meeting_followup_20260905','output/iNKT_discovery_20260905','output/iNKT_reproduction_deck/20260830_C5_paper_Fig3DEF_followup','iNKT_by_date/2026-09-19']
files={};links=[];groups={}
for root in roots:
 group=0
 for base,dirs,names in os.walk(r/root,followlinks=False):
  dirs[:]=[d for d in dirs if d not in ['__pycache__','.ipynb_checkpoints']]
  for name in names:
   p=Path(base)/name
   if p.is_symlink():links.append({'path':str(p.relative_to(r)),'target':str(p.resolve())});continue
   if not p.is_file():continue
   rel=str(p.relative_to(r));h=hashlib.sha256()
   with p.open('rb') as f:
    for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
   n=p.stat().st_size;files[rel]={'bytes':n,'sha256':h.hexdigest()};group+=n
 groups[root]=group
for p in r.iterdir():
 if p.is_file() and (p.suffix=='.sh' or p.name in ['README.md','pyproject.toml','uv.lock']):
  files[p.name]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
print(json.dumps({'source_root':str(r),'snapshot_date':'2026-09-25','files':files,'groups':groups,'skipped_symlinks':links,'total_bytes':sum(v['bytes'] for v in files.values()),'file_count':len(files),'excluded':['.git','.venv','unrelated projects','older outputs beyond selected dependency directories']},indent=2))
