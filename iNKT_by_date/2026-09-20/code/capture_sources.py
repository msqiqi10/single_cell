"""Re-run only figure builders and export individual artists, never raster slide screenshots.
Scientific analyses are read-only. Native text/table scene records accompany source panels.
"""
from __future__ import annotations
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
import sys, json, shutil, inspect, hashlib, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.colors import to_rgba
from PIL import Image
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parents[1];PKG=OUT/'package';WORK=OUT/'work';SCENES=[]
sys.path.insert(0,str(ROOT/'notebooks/scripts/iNKT'))
sys.path.insert(0,str(ROOT/'iNKT_by_date/2026-09-19/code'))

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def rel(p):return str(Path(p).relative_to(PKG))
def rgba(x):
 try:return list(to_rgba(x))
 except (ValueError,TypeError):return [0,0,0,1]
def box(b,fig):
 dpi=fig.dpi;h=fig.get_figheight();return [float(b.x0/dpi),float(h-b.y1/dpi),float(b.width/dpi),float(b.height/dpi)]
def txt_record(t,fig,renderer,role):
 b=t.get_window_extent(renderer);return {'type':'text','text':t.get_text(),'box':box(b,fig),'font_size':float(t.get_fontsize()),'font_family':'Arial','bold':t.get_fontweight() in ['bold','heavy','semibold','demibold',700,800,900],'italic':t.get_fontstyle()=='italic','color':rgba(t.get_color()),'rotation':float(t.get_rotation()),'role':role,'horizontal_alignment':t.get_ha()}

def capture(fig,series,page,info,notes=''):
 ident=f'{series}_{page:02}';print(time.strftime('%H:%M:%S'),'CAPTURE',ident,info['title'],flush=True)
 fig.canvas.draw();renderer=fig.canvas.get_renderer();width,height=fig.get_size_inches();elements=[]
 # Preserve same-source PDF and previews only for PDF/QA. They are never inserted in PPTX.
 fig.savefig(PKG/'assets/pdf_pages'/f'{ident}.pdf')
 fig.savefig(PKG/'assets/page_previews'/f'{ident}.png',dpi=110)
 text_artists=[t for t in fig.texts if t.get_visible() and t.get_text()]
 for ax in fig.axes:
  text_artists += [t for t in [ax.title,ax._left_title,ax._right_title] if t.get_visible() and t.get_text()]
 for t in text_artists:
  role='plot_title' if t.axes is not None else 'figure_text'
  if t.get_text()==info.get('title'):role='title'
  if t.get_text()==str(page) and t.get_position()[0]>.9:role='page_number'
  elements.append(txt_record(t,fig,renderer,role))
 tables=[];table_axes=set()
 for ax in fig.axes:
  for table in ax.tables:
   cells=table.get_celld();nr=max(r for r,c in cells)+1;nc=max(c for r,c in cells)+1
   bb=table.get_window_extent(renderer);tb={'type':'table','box':box(bb,fig),'rows':nr,'cols':nc,'cells':[],'column_widths':[cells[(0,c)].get_window_extent(renderer).width/fig.dpi for c in range(nc)],'row_heights':[cells[(r,0)].get_window_extent(renderer).height/fig.dpi for r in range(nr)]}
   for (r,c),cell in sorted(cells.items()):
    t=cell.get_text();tb['cells'].append({'r':r,'c':c,'text':t.get_text(),'font_size':float(t.get_fontsize()),'bold':t.get_fontweight() in ['bold','semibold',700],'color':rgba(t.get_color()),'fill':list(cell.get_facecolor()),'align':t.get_ha()})
   elements.append(tb);table_axes.add(ax);tables.append(table)
   pd.DataFrame([[cells[(r,c)].get_text().get_text() for c in range(nc)] for r in range(nr)]).to_csv(PKG/'data'/f'{ident}_editable_table_{len(tables)}.csv',index=False,header=False)
 # Hide other axes and all source figure text; each panel is exported from its own artists.
 axes=list(fig.axes);vis={ax:ax.get_visible() for ax in axes};tvis={t:t.get_visible() for t in text_artists}
 for t in text_artists:t.set_visible(False)
 for ax in axes:ax.set_visible(False)
 for i,ax in enumerate(axes):
  if not vis[ax] or ax in table_axes:continue
  # Pure text-only axes become native editable text boxes (e.g. driver-gene callouts).
  if not ax.axison and not ax.images and not ax.collections and not ax.lines and not ax.patches:
   for t in ax.texts:
    if t.get_visible() and t.get_text():elements.append(txt_record(t,fig,renderer,'body_text'))
   continue
  ax.set_visible(True)
  fig.canvas.draw();renderer=fig.canvas.get_renderer()
  if hasattr(ax,'_direct_source_path'):
   p=Path(ax._direct_source_path);dest=PKG/'assets/panels'/f'{ident}_panel_{i+1:02}{p.suffix}';shutil.copy2(p,dest)
   b=ax.images[0].get_window_extent(renderer);elements.append({'type':'image','box':box(b,fig),'png':rel(dest),'source_mode':'original_image_file','source_path':str(p.relative_to(ROOT)),'source_sha256':sha(p),'role':'source_figure'})
  else:
   bb=ax.get_tightbbox(renderer)
   if bb is None:ax.set_visible(False);continue
   bb=Bbox.from_extents(max(0,bb.x0-3),max(0,bb.y0-3),min(width*fig.dpi,bb.x1+3),min(height*fig.dpi,bb.y1+3))
   if bb.width<=0 or bb.height<=0:ax.set_visible(False);continue
   png=PKG/'assets/panels'/f'{ident}_panel_{i+1:02}.png';svg=png.with_suffix('.svg')
   bbox=bb.transformed(fig.dpi_scale_trans.inverted())
   fig.savefig(png,dpi=220,bbox_inches=bbox,pad_inches=0,transparent=True)
   with matplotlib.rc_context({'svg.fonttype':'path'}):fig.savefig(svg,bbox_inches=bbox,pad_inches=0,transparent=True)
   elements.append({'type':'image','box':box(bb,fig),'png':rel(png),'svg':rel(svg),'source_mode':'individual_matplotlib_axes','axes_index':i,'role':'plot_panel'})
  ax.set_visible(False)
 for ax,v in vis.items():ax.set_visible(v)
 for t,v in tvis.items():t.set_visible(v)
 frame=next((f for f in inspect.stack() if f.function not in ['capture','finish','stop_after','old_finish','new_finish']),None)
 scene={'id':ident,'series':series,'original_page':page,'width':float(width),'height':float(height),'title':info['title'],'subtitle':info.get('subtitle',''),'reference':info.get('reference',''),'note':info.get('note',''),'notes_zh':notes,'elements':elements,'pdf_page':f'assets/pdf_pages/{ident}.pdf','preview':f'assets/page_previews/{ident}.png','source_builder_function':frame.function if frame else ''}
 SCENES.append(scene)
 (PKG/'slides.json').write_text(json.dumps(SCENES,ensure_ascii=False,indent=2))
 return scene

def old_figures():
 import build_inkt_reference_aligned_presentation as m
 m.OUT=WORK/'old';m.PAGES=[];m.SOURCES={};m.EXPORTS=[];m.CHECKS=[]
 def finish(fig,pdf):
  i=len(m.PAGES)+1;info=fig._page_info;capture(fig,'old',i,info,'原综合报告第'+str(i)+'页。'+info['note'])
  pdf.savefig(fig);m.PAGES.append({'page':i,**info});plt.close(fig)
 def image_page(path,title_text,subtitle,reference,pdf,note='Additional analysis; retain its original metric and interpretation.'):
  f=m.new(title_text,subtitle,reference,note);ax=f.add_axes([.035,.13,.93,.72]);ax.imshow(Image.open(m.record(path)));ax.axis('off');ax._direct_source_path=str(path);finish(f,pdf)
 m.finish=finish;m.image_page=image_page;m.package=lambda:None;m.write_audit=lambda:None
 m.main()
 shutil.copytree(m.OUT/'tables',PKG/'data/previous_display_tables',dirs_exist_ok=True)
 (PKG/'notes/previous_source_hashes.json').write_text(json.dumps(m.SOURCES,indent=2))
 assert len([s for s in SCENES if s['series']=='old'])==37

def new_figures():
 import presentation as m
 base=ROOT/'iNKT_by_date/2026-09-19';w=WORK/'new'
 for d in ['results/tables','results/figures','results/objects','presentation','notes']:(w/d).mkdir(parents=True,exist_ok=True)
 shutil.copytree(base/'results/tables',w/'results/tables',dirs_exist_ok=True);shutil.copytree(base/'results/objects',w/'results/objects',dirs_exist_ok=True)
 m.OUT=w;m.RES=w/'results';m.PRES=w/'presentation';m.FIG=w/'results/figures';m.SRC=base/'sources';m.PAGES=[]
 def finish(fig,pdf,notes=''):
  i=len(m.PAGES)+1;capture(fig,'new',i,fig._info,notes);m.PAGES.append({'page':i,**fig._info,'notes_zh':notes,'figure':str((PKG/'assets/page_previews'/f'new_{i:02}.png').relative_to(OUT))});plt.close(fig)
 m.finish=finish
 # Source main() assembles an old image-only deck after plotting; stop immediately after the final figure.
 class FiguresCaptured(Exception):pass
 orig_finish=m.finish
 def stop_after(fig,pdf,notes=''):
  orig_finish(fig,pdf,notes)
  if len(m.PAGES)==19:raise FiguresCaptured()
 m.finish=stop_after
 try:m.main()
 except FiguresCaptured:pass
 assert len([s for s in SCENES if s['series']=='new'])==19
 shutil.copytree(base/'results/tables',PKG/'data/current_analysis_tables',dirs_exist_ok=True)
 for name in ['input_manifest.json','GO_provenance.json','download_manifest.json']:
  shutil.copy2(base/'sources'/name,PKG/'notes'/name)


def custom_pages():
 def new(title,subtitle):
  f=plt.figure(figsize=(16,9));f.text(.045,.94,title,fontsize=25,weight='bold',va='top');f.text(.045,.875,subtitle,fontsize=13,va='top',color='#555');f.text(.045,.045,'20 September 2026 | Existing iNKT analyses; no new statistical analysis in this packaging step.',fontsize=9,color='#555');return f
 f=new('iNKT: cell states, cytotoxic programs and functional changes','Integrated editable presentation | main talk + complete source-page appendix')
 f.text(.055,.74,'Research question',fontsize=19,weight='bold');f.text(.055,.67,'How does the T2-associated response differ across tissues and iNKT states?',fontsize=17)
 f.text(.055,.55,'Presentation sequence',fontsize=19,weight='bold');f.text(.055,.47,'1  Cohort, cell-state markers and tissue composition\n\n2  Differential expression, cytotoxic modules and GO functions\n\n3  Local state ordering, robustness and bounded interpretation',fontsize=17,va='top')
 f.text(.055,.19,'15,532 cells | 3 tissues | Ctrl / T2 | 9 original refined clusters\nOne sample label per tissue x condition: all biological conclusions remain exploratory.',fontsize=13)
 capture(f,'custom',1,{'title':'iNKT: cell states, cytotoxic programs and functional changes','subtitle':'Integrated editable presentation','reference':'Blood Fig. 3; PAGER-scFGA Fig. 5; Borra et al. Fig. 4c','note':'Main talk followed by an appendix preserving every previous and current source page.'},'主报告按数据与细胞状态、差异与功能、轨迹与解释组织。原37页和新增19页均保留，部分进入主报告，其余进入附录。');plt.close(f)
 f=new('Interpretation: localized effector changes with shared metabolic signals','Separate direct observations, working hypotheses and unresolved evidence')
 rows=[['Bone marrow C4','Core score +0.086; q=0.025. Gzmb detection 29% to 44%.','Selective Gzmb-associated effector transcriptional shift.'],['Spleen C3','Core score -0.023; q=0.0015. Prf1 and Nkg7 decrease.','Small reduction in effector-associated transcripts.'],['Both populations','T2-up genes enriched for oxidative phosphorylation.','Respiration-related gene expression shifts; no metabolic flux measured.'],['Limits','Broader E1 signature and killing-related GO are not significant.','No demonstrated change in killing capacity or animal-level treatment effect.']]
 ax=f.add_axes([.045,.30,.91,.46]);ax.axis('off');tb=ax.table(cellText=rows,colLabels=['Population / scope','Observation','Interpretation'],colWidths=[.17,.42,.41],cellLoc='left',colLoc='left',bbox=[0,0,1,1]);tb.auto_set_font_size(False);tb.set_fontsize(11)
 for (r,c),cell in tb.get_celld().items():cell.set_facecolor('#eaf0f5' if r==0 else 'white');cell.set_edgecolor('#ccc');cell.get_text().set_wrap(True)
 f.text(.055,.20,'Both core effects persist under QC matching and overlapping stable clusters.\nThese same-data checks support candidates; independent animals and functional assays remain necessary.',fontsize=13)
 capture(f,'custom',2,{'title':'Interpretation: localized effector changes with shared metabolic signals','subtitle':'Separate direct observations, working hypotheses and unresolved evidence','reference':'20260919 numeric results; gene-level follow-up reviewed 20260919','note':'No exhaustion, lineage direction or killing-function claim is established.'},'骨髓C4主要是Gzmb相关变化；脾脏C3效应较小。两群均有氧化磷酸化基因富集，但细胞毒性分数方向不同。扩展E1与杀伤GO不显著，因此限制在转录状态变化。');plt.close(f)


def main():
 old_figures();new_figures();custom_pages()
 refs=['docs/pager-scFGA.pdf','docs/Borra_et_al_2026_GAFA_CML_NK.pdf','docs/blooda_adv-2024-014592-main.pdf']
 for p in refs:shutil.copy2(ROOT/p,PKG/'references'/Path(p).name)
 for p in [ROOT/'notebooks/scripts/iNKT/build_inkt_reference_aligned_presentation.py',ROOT/'iNKT_by_date/2026-09-19/code/presentation.py']:
  name='source_builder_previous.py' if 'reference_aligned' in p.name else 'source_builder_current.py';shutil.copy2(p,PKG/'code'/name)
 (OUT/'logs/capture.complete.json').write_text(json.dumps({'pages':len(SCENES),'panels':sum(e['type']=='image' for s in SCENES for e in s['elements']),'native_texts':sum(e['type']=='text' for s in SCENES for e in s['elements']),'native_tables':sum(e['type']=='table' for s in SCENES for e in s['elements'])},indent=2))
 print('CAPTURE COMPLETE',len(SCENES),flush=True)
if __name__=='__main__':main()
