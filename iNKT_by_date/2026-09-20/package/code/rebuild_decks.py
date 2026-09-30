"""Standalone PPTX rebuild from scene JSON and original independent figure assets.
Requires python-pptx and pypdf. PDF companion uses same-source vector figure pages;
after editing the PPTX, export PDF in PowerPoint to synchronize revisions.
"""
from __future__ import annotations
import json, hashlib, csv, shutil, zipfile, io
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN, MSO_AUTO_SIZE
from pptx.opc.package import Part
from pptx.opc.constants import RELATIONSHIP_TYPE as RT
from pptx.oxml import parse_xml
from pypdf import PdfReader,PdfWriter
BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'decks';OUT.mkdir(exist_ok=True)
SVG_NS='http://schemas.microsoft.com/office/drawing/2016/SVG/main'
DML_NS='http://schemas.openxmlformats.org/drawingml/2006/main'
REL_NS='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
MAIN=['custom_01','old_02','old_03','old_05','old_06','old_07','old_08','old_09','old_10','old_22','new_02','new_03','new_04','new_05','new_06','new_07','new_08','new_09','new_10','new_11','new_12','new_17','new_13','new_14','new_15','new_18','new_16','custom_02']

def color(rgba):return RGBColor(*(int(round(max(0,min(1,c))*255)) for c in rgba[:3]))
def text_style(run,e):
 run.font.name='Arial';run.font.size=Pt(e.get('font_size',11));run.font.bold=e.get('bold',False);run.font.italic=e.get('italic',False);run.font.color.rgb=color(e.get('color',[0,0,0,1]))
def add_text(slide,e,scene,page):
 x,y,w,h=e['box'];text=str(page) if e['role']=='page_number' else e['text']
 # A modest extra width accommodates Arial vs DejaVu metrics without flattening text.
 width=min(scene['width']-max(0,x),max(.14,w*1.04+.08));height=max(.18,h+.08)
 shape=slide.shapes.add_textbox(Inches(max(0,x)),Inches(max(0,y-.01)),Inches(width),Inches(height));shape.name=e['role']+' | '+text.replace('\n',' ')[:65]
 tf=shape.text_frame;tf.clear();tf.word_wrap=False;tf.margin_left=0;tf.margin_right=0;tf.margin_top=0;tf.margin_bottom=0;tf.vertical_anchor=MSO_ANCHOR.TOP;tf.auto_size=MSO_AUTO_SIZE.NONE
 for i,line in enumerate(text.split('\n')):
  p=tf.paragraphs[0] if i==0 else tf.add_paragraph();p.space_before=Pt(0);p.space_after=Pt(0);p.line_spacing=Pt(e.get('font_size',11)*1.19)
  p.alignment={'left':PP_ALIGN.LEFT,'center':PP_ALIGN.CENTER,'right':PP_ALIGN.RIGHT}.get(e.get('horizontal_alignment'),PP_ALIGN.LEFT)
  run=p.add_run();run.text=line;text_style(run,e)
 shape.rotation=e.get('rotation',0)
 return shape

def add_table(slide,e):
 x,y,w,h=e['box'];shape=slide.shapes.add_table(e['rows'],e['cols'],Inches(x),Inches(y),Inches(w),Inches(h));shape.name='Editable numerical / explanatory table';t=shape.table
 for i,v in enumerate(e['column_widths']):t.columns[i].width=Inches(v)
 for i,v in enumerate(e['row_heights']):t.rows[i].height=Inches(v)
 for ce in e['cells']:
  cell=t.cell(ce['r'],ce['c']);cell.text=ce['text'];cell.fill.solid();cell.fill.fore_color.rgb=color(ce['fill']);cell.margin_left=Inches(.06);cell.margin_right=Inches(.04);cell.margin_top=Inches(.025);cell.margin_bottom=Inches(.025);cell.vertical_anchor=MSO_ANCHOR.MIDDLE
  for para in cell.text_frame.paragraphs:
   para.space_before=Pt(0);para.space_after=Pt(0);para.alignment={'left':PP_ALIGN.LEFT,'center':PP_ALIGN.CENTER,'right':PP_ALIGN.RIGHT}.get(ce['align'],PP_ALIGN.LEFT)
   for run in para.runs:text_style(run,ce)
  cell.text_frame.word_wrap=True
  tcpr=cell._tc.get_or_add_tcPr()
  for edge in ['lnL','lnR','lnT','lnB']:
   tcpr.append(parse_xml(f'<a:{edge} xmlns:a="{DML_NS}" w="6350"><a:solidFill><a:srgbClr val="DDDDDD"/></a:solidFill><a:prstDash val="solid"/></a:{edge}>'))
 return shape

def add_image(slide,e,cache):
 x,y,w,h=e['box'];shape=slide.shapes.add_picture(str(BASE/e['png']),Inches(x),Inches(y),width=Inches(w),height=Inches(h));shape.name=e['role']+' | '+Path(e['png']).stem
 if e.get('svg'):
  svg_path=BASE/e['svg'];package=slide.part.package
  if str(svg_path) not in cache:
   part=Part(package.next_partname('/ppt/media/panel%d.svg'),'image/svg+xml',package,svg_path.read_bytes());cache[str(svg_path)]=part
  rid=slide.part.relate_to(cache[str(svg_path)],RT.IMAGE)
  ext=parse_xml(f'<a:extLst xmlns:a="{DML_NS}" xmlns:asvg="{SVG_NS}" xmlns:r="{REL_NS}"><a:ext uri="{{96DAC541-7B7A-43D3-8B79-37D633B846F1}}"><asvg:svgBlip r:embed="{rid}"/></a:ext></a:extLst>')
  shape._pic.blipFill.blip.append(ext)
 return shape

def build(name,scenes,main_count=None):
 prs=Presentation();prs.slide_width=Inches(16);prs.slide_height=Inches(9);cache={};rows=[]
 prs.core_properties.title='iNKT — editable source-element presentation';prs.core_properties.subject='Native text and tables, individual SVG/PNG figure panels';prs.core_properties.author='iNKT analysis project';prs.core_properties.comments='Full-slide images are excluded. All figure panels have source assets and provenance.'
 for number,s in enumerate(scenes,1):
  slide=prs.slides.add_slide(prs.slide_layouts[6]);counts={'text':0,'table':0,'image':0}
  # Graphics first, then fully editable text and tables preserve layering.
  for kind in ['image','table','text']:
   for e in s['elements']:
    if e['type']!=kind:continue
    if kind=='image':add_image(slide,e,cache)
    elif kind=='table':add_table(slide,e)
    else:add_text(slide,e,s,number)
    counts[kind]+=1
  section='Main talk' if main_count is None or number<=main_count else 'Appendix'
  slide.notes_slide.notes_text_frame.text='\n'.join([s['title'],s.get('notes_zh',''),s.get('reference',''),s.get('note',''),f'Source page: {s["id"]}; {section}',f'Independent graphics: assets/panels/{s["id"]}_*', 'Native text and tables can be edited directly. SVG/PNG panels can be moved, resized or replaced independently; numeric source tables are in data/.'])
  rows.append({'slide':number,'section':section,'source_page':s['id'],'title':s['title'],**counts})
 path=OUT/(name+'.pptx');prs.save(path)
 pdf=PdfWriter()
 for number,s in enumerate(scenes,1):
  source=PdfReader(BASE/s['pdf_page']);page=source.pages[0]
  numbers=[e for e in s['elements'] if e['type']=='text' and e['role']=='page_number']
  if numbers:
   fig=plt.figure(figsize=(16,9));fig.patch.set_alpha(0)
   for e in numbers:
    x,y,w,h=e['box'];fig.add_artist(Rectangle(((x-.03)/16,1-(y+h+.03)/9),(w+.15)/16,(h+.06)/9,transform=fig.transFigure,facecolor='white',edgecolor='none'))
    fig.text((x+w)/16,1-(y+h*.85)/9,str(number),fontsize=e['font_size'],ha='right',va='baseline',color='#555555')
   stream=io.BytesIO();fig.savefig(stream,format='pdf',transparent=True);plt.close(fig);stream.seek(0);page.merge_page(PdfReader(stream).pages[0])
  pdf.add_page(page)
 pdf.add_metadata({'/Title':'iNKT — '+name,'/Subject':'Same-source companion PDF; export from PowerPoint after editing slides.'})
 with (OUT/(name+'.pdf')).open('wb') as f:pdf.write(f)
 with (OUT/(name+'_index.csv')).open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
 print(name,len(scenes),'slides',sum(r['text'] for r in rows),'native text boxes',sum(r['table'] for r in rows),'native tables',sum(r['image'] for r in rows),'independent panels',flush=True)
 return rows

def main():
 scenes=json.loads((BASE/'slides.json').read_text());byid={s['id']:s for s in scenes}
 assert len(byid)==58
 old=[byid[f'old_{i:02}'] for i in range(1,38)];new=[byid[f'new_{i:02}'] for i in range(1,20)]
 integrated=[byid[x] for x in MAIN]+[s for s in old+new if s['id'] not in MAIN]
 assert len(integrated)==58 and len(set(s['id'] for s in integrated))==58
 build('iNKT_previous_37_editable',old)
 build('iNKT_cytotoxicity_GO_19_editable',new)
 rows=build('iNKT_integrated_58_editable',integrated,len(MAIN))
 build('iNKT_main_talk_28_editable',integrated[:len(MAIN)])
 (BASE/'notes/deck_order.json').write_text(json.dumps({'main_talk':MAIN,'integrated':[s['id'] for s in integrated]},indent=2))
 return rows
if __name__=='__main__':main()
