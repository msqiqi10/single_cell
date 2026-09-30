"""Create a reviewable report package; keep expression objects/model weights local."""
from pathlib import Path
import json,shutil,hashlib,re,zipfile,ast
ROOT=Path(__file__).resolve().parents[1];RUN=ROOT/'runs/meeting_followup'
DEST=ROOT.parent/'outputs/iNKT_20260925_Rob_Yue_results'
DEST.mkdir(parents=True,exist_ok=True)
validation=json.loads((RUN/'validation.json').read_text());assert validation['passed']
for sub in ['tables','network','figures','existing_evidence','gene_lists','models']:
    shutil.copytree(RUN/sub,DEST/sub,dirs_exist_ok=True)
for p in RUN.iterdir():
    if p.is_file():shutil.copy2(p,DEST/p.name)
shutil.copytree(ROOT/'src',DEST/'code',dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__'))
for p in (ROOT/'src').glob('*.py'):ast.parse(p.read_text(),str(p))
(DEST/'provenance').mkdir(exist_ok=True)
for name in ['gnn_adapter.json','gnn_adapter.patch','semantic_model.json','local_verification.json','remote_manifest.json']:
    shutil.copy2(ROOT/'provenance'/name,DEST/'provenance'/name)
shutil.copy2(ROOT/'requirements.lock.txt',DEST/'requirements.lock.txt')
shutil.copy2(ROOT/'run_local.sh',DEST/'code/run_local.sh')
sources=[(ROOT.parent/'outputs/2026-09-24_会议核对依据/candidate_coverage.csv','legacy_page26_28_candidate_coverage.csv'),
         (ROOT.parent/'outputs/iNKT_20260925_followup/tables/Rob_focus_existing_robust_pathways.csv','legacy_robust_pathway_audit.csv'),
         (ROOT.parent/'outputs/iNKT_20260925_followup/tables/sample_design_registry.csv','sample_design_registry.csv')]
for source,name in sources:shutil.copy2(source,DEST/'existing_evidence'/name)
packaged_validation=json.loads((DEST/'validation.json').read_text())
packaged_validation['source_hash_manifest']='provenance/remote_manifest.json'
packaged_validation['original_run_validation_path']=str(RUN/'validation.json')
(DEST/'validation.json').write_text(json.dumps(packaged_validation,indent=2))
# Check every relative Markdown link in the main deliverables.
links=[]
for md in [DEST/'README.md',DEST/'BRIEF_EN.md',DEST/'EXPERIMENT_CANDIDATES.md',DEST/'VALIDATION_NOTES.md']:
    for target in re.findall(r'\]\(([^)]+)\)',md.read_text()):
        if target.startswith(('https://','http://','#')):continue
        p=(md.parent/target.split('#')[0]);assert p.exists(),(md.name,target)
        links.append(target)
manifest={'validation_passed':True,'relative_markdown_links_checked':len(links),
          'scope':'Report/evidence/code snapshot; local_inkt contains the actual environment, inputs, expression objects and public model weights',
          'files':{str(p.relative_to(DEST)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
                   for p in sorted(DEST.rglob('*')) if p.is_file() and p.name!='PACKAGE_MANIFEST.json'}}
(DEST/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
archive=DEST.with_suffix('.zip')
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(DEST.rglob('*')):
        if p.is_file():z.write(p,arcname=str(Path(DEST.name)/p.relative_to(DEST)))
with zipfile.ZipFile(archive) as z:assert z.testzip() is None
print(json.dumps(dict(report=str(DEST/'README.md'),archive=str(archive),files=len(manifest['files']),bytes=archive.stat().st_size,links_checked=len(links)),indent=2))
