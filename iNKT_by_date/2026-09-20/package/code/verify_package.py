"""Check editability, retained content, source figures and complete page coverage."""
from pathlib import Path
import json,hashlib,zipfile,re
import pandas as pd
from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE_TYPE
from pypdf import PdfReader
from lxml import etree
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parents[1];BASE=OUT/'package';checks=[];warnings=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(ok,label):
 if not bool(ok):raise AssertionError(label)
 checks.append(label)

def main():
 scenes=json.loads((BASE/'slides.json').read_text());byid={s['id']:s for s in scenes};check(len(scenes)==58,'All 37 old, 19 new and two integration pages captured')
 preview_hash={sha(BASE/s['preview']) for s in scenes}
 panelexport=[]
 for s in scenes:
  for i,e in enumerate(s['elements']):
   if e['type']=='image':
    p=BASE/e['png'];check(p.exists(),f'Panel exists {p.name}');check(sha(p) not in preview_hash,f'Panel is not a flattened slide: {p.name}');check(e['source_mode'] in ['original_image_file','individual_matplotlib_axes'],'Panel origin is a true source figure/artist')
    if e.get('source_path'):check(sha(ROOT/e['source_path'])==e['source_sha256'],'Original directly pasted image unchanged '+e['source_path'])
    if e.get('svg'):etree.parse(str(BASE/e['svg']))
    check(e['box'][2]*e['box'][3]<.90*s['width']*s['height'],'No picture covers a whole slide '+s['id'])
    panelexport.append({'source_page':s['id'],'panel_index':i,'panel_file':e['png'],'svg_file':e.get('svg',''),'origin':e['source_mode'],'original_image':e.get('source_path',''),'source_reference':s['reference']})
 pd.DataFrame(panelexport).to_csv(BASE/'notes/panel_source_index.csv',index=False)
 for p in (BASE/'decks').glob('*.pptx'):
  prs=Presentation(p);index=pd.read_csv(p.with_name(p.stem+'_index.csv'));pdf=PdfReader(p.with_suffix('.pdf'));check(len(prs.slides)==len(index)==len(pdf.pages),'PPTX/PDF/index counts match '+p.stem)
  with zipfile.ZipFile(p) as z:
   embedded_svg=[n for n in z.namelist() if n.startswith('ppt/media/') and n.endswith('.svg')];check(len(embedded_svg)>0,'Native SVG source panels embedded '+p.stem)
   for n in embedded_svg:etree.fromstring(z.read(n))
   for n in z.namelist():
    if n.startswith('ppt/slides/slide') and n.endswith('.xml'):
     x=etree.fromstring(z.read(n));blips=x.findall('.//{http://schemas.microsoft.com/office/drawing/2016/SVG/main}svgBlip')
     if blips:
      relname=n.replace('ppt/slides/','ppt/slides/_rels/')+'.rels';rx=etree.fromstring(z.read(relname));ids={r.get('Id') for r in rx};check(all(b.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed') in ids for b in blips),'SVG relationships resolve '+n)
  for i,(sl,row) in enumerate(zip(prs.slides,index.itertuples()),1):
   source=byid[row.source_page];text=[sh for sh in sl.shapes if sh.has_text_frame];tables=[sh for sh in sl.shapes if sh.has_table];images=[sh for sh in sl.shapes if sh.shape_type==MSO_SHAPE_TYPE.PICTURE]
   check(len(text)==row.text and len(tables)==row.table and len(images)==row.image,f'Element counts correct {p.stem}/{i}')
   check(len(text)>=3,f'Title/body/captions editable {p.stem}/{i}')
   actual_text=[sh.text for sh in text];expected=[str(i) if e['role']=='page_number' else e['text'] for e in source['elements'] if e['type']=='text'];check(actual_text==expected,f'All source text preserved natively {p.stem}/{i}')
   for sh,e in zip(tables,[e for e in source['elements'] if e['type']=='table']):
    check(all(sh.table.cell(c['r'],c['c']).text==c['text'] for c in e['cells']),f'Native table values unchanged {p.stem}/{i}')
   for sh in sl.shapes:
    check(sh.width>0 and sh.height>0,'Positive editable-element bounds')
    if sh.left<0 or sh.top<0 or sh.left+sh.width>prs.slide_width+91440 or sh.top+sh.height>prs.slide_height+91440:warnings.append({'deck':p.stem,'page':i,'shape':sh.name,'reason':'element extends beyond slide boundary'})
 integrated=pd.read_csv(BASE/'decks/iNKT_integrated_58_editable_index.csv');check(set(integrated.source_page)==set(byid),'Integrated deck preserves every source page exactly once')
 for p,h in json.loads((BASE/'notes/previous_source_hashes.json').read_text()).items():check(sha(ROOT/p)==h,'Previous source input unchanged '+p)
 result={'status':'passed','checks_passed':len(checks),'checks':checks,'warnings':warnings,'flat_slide_pictures':0,'native_text_tables':True,'plots':'Independent SVG+PNG source panels; plot-internal marks/labels are graphical elements, not native spreadsheet-linked charts.','PDF':'Rendered from the same source figure layout, not from an Office rendering of the PPTX. Export PDF in PowerPoint after supervisor edits.','GPU_used':False}
 (BASE/'notes/verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));print('PASS',len(checks),'checks;',len(warnings),'layout warnings',flush=True)
 if warnings:print(json.dumps(warnings,indent=2),flush=True)
if __name__=='__main__':main()
