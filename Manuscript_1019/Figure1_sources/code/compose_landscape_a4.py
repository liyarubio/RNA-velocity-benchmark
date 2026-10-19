"""Recompose existing vector Figure 1 and caption on landscape A4.
No scientific artwork is redrawn. Source and superseded deliverables are retained.
"""
from pathlib import Path
import sys,json,hashlib,shutil,re,base64
from collections import Counter
sys.path.insert(0,'/data1/liyaru/RNA_velocity_materials/Manuscript_1010/Code/vendor')
import pymupdf as fitz
from bs4 import BeautifulSoup
B=Path(__file__).resolve().parents[2]
R=B/'Reports/figure1_landscape_20261019';BEFORE=R/'before'
BEFORE.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
targets=['submission_figures/Figure01.pdf','submission_figures/Main_Figures_1-7_combined.pdf','figures/main/Figure01.pdf','figures/main/Figure01.png','figures/main/Figure01.svg','release_status.json','manifest.json','checksums.json','Reports/editorial_1019/validation.json']
for rel in targets:
 dst=BEFORE/rel
 if not dst.exists():dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(B/rel,dst)
source=BEFORE/'submission_figures/Figure01.pdf'
original=fitz.open(source);sp=original[0]
caption_blocks=[b for b in sp.get_text('blocks') if b[1]>330 and 'Fig.' in b[4]]
assert len(caption_blocks)==1
caption=' '.join(caption_blocks[0][4].replace('\xad','-').split())
title='Fig. 1 | Chronological landscape of RNA velocity algorithms.'
assert caption.startswith(title)
clip=fitz.Rect(32,54,563.2756,324)
W,H=297*72/25.4,210*72/25.4;margin=32.;top=28.
scale=(W-2*margin)/clip.width
art=fitz.Rect(margin,top,W-margin,top+clip.height*scale)
out=fitz.open();p=out.new_page(width=W,height=H)
p.show_pdf_page(art,original,0,clip=clip)
# Clip embedded objects explicitly, preserving all objects within the art region.
xref=p.get_contents()[-1]
out.update_stream(xref,f'q {art.x0} {H-art.y1} {art.width} {art.height} re W n\n'.encode()+out.xref_stream(xref)+b'\nQ\n')
fonts={False:fitz.Font(fontfile='/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'),True:fitz.Font(fontfile='/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf')}
for bold,name in [(False,'CaptionArial'),(True,'CaptionArialBold')]:p.insert_font(fontname=name,fontbuffer=fonts[bold].buffer)
size=10;tokens=[(w,True) for w in title.split()]+[(w,False) for w in caption[len(title):].split()]
lines=[];line=[];width=0;space=fonts[False].text_length(' ',fontsize=size)
for word,bold in tokens:
 length=fonts[bold].text_length(word,fontsize=size)
 if line and width+space+length>W-2*margin:lines.append(line);line=[];width=0
 width+=length+(space if line else 0);line.append((word,bold,length))
if line:lines.append(line)
y=art.y1+23
for i,line in enumerate(lines):
 gap=((W-2*margin)-sum(t[2] for t in line))/(len(line)-1) if i<len(lines)-1 and len(line)>1 else space
 x=margin
 for word,bold,length in line:
  p.insert_text((x,y),word,fontname='CaptionArialBold' if bold else 'CaptionArial',fontsize=size)
  x+=length+gap
 y+=13
assert y<H-28
out.set_metadata({'title':'Figure 1 | Chronological landscape of RNA velocity algorithms','subject':'Manuscript 1019; original vector artwork recomposed on landscape A4'})
for rel in ['submission_figures/Figure01.pdf','figures/main/Figure01.pdf']:
 temp=B/(rel+'.tmp');out.save(temp,garbage=3,deflate=True);temp.replace(B/rel)
# Headerless artwork assets for manuscript readers: same visible art, larger render.
artdoc=fitz.open();ap=artdoc.new_page(width=art.width,height=art.height);ap.show_pdf_page(ap.rect,original,0,clip=clip)
ap.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False).save(str(B/'figures/main/Figure01.png'))
(B/'figures/main/Figure01.svg').write_text(ap.get_svg_image())
p.get_pixmap(matrix=fitz.Matrix(1.25,1.25),alpha=False).save(str(R/'Figure01_landscape_A4_preview.png'))
# All forms of the full manuscript use the same replacement artwork.
stand=B/'Manuscript_1019_NM_bilingual_standalone.html'
soup=BeautifulSoup(stand.read_text(),'html.parser')
base=BeautifulSoup((B/'Manuscript_1019_NM_bilingual.html').read_text(),'html.parser')
which=[i for i,im in enumerate(base.select('img')) if im.get('src')=='figures/main/Figure01.png']
assert len(which)==1
soup.select('img')[which[0]]['src']='data:image/png;base64,'+base64.b64encode((B/'figures/main/Figure01.png').read_bytes()).decode()
stand.write_text(str(soup))
# Replace page 1 only. Keep all nine subsequent pages and their vector content.
old=fitz.open(BEFORE/'submission_figures/Main_Figures_1-7_combined.pdf')
merged=fitz.open();merged.insert_pdf(out);merged.insert_pdf(old,from_page=1)
merged.set_toc(old.get_toc());merged.set_metadata({'title':'Manuscript 1019 — Main Figures 1–7 with legends','subject':'Figure 1 landscape A4; remaining pages preserved'})
combined=B/'submission_figures/Main_Figures_1-7_combined.pdf'
temp=combined.with_suffix('.tmp.pdf');merged.save(temp,garbage=3,deflate=True);temp.replace(combined)
new=fitz.open(combined)
unchanged=[]
for i in range(1,len(old)):
 assert old[i].get_pixmap(alpha=False).samples==new[i].get_pixmap(alpha=False).samples,i
 unchanged.append(i+1)
assert len(new)==10
published=fitz.open(B/'submission_figures/Figure01.pdf')
assert new[0].get_pixmap(alpha=False).samples==published[0].get_pixmap(alpha=False).samples
srcwords=Counter(w[4] for w in sp.get_text('words',clip=clip))
newwords=Counter(w[4] for w in p.get_text('words',clip=art))
assert srcwords==newwords,(srcwords-newwords,newwords-srcwords)
after_caption=' '.join(p.get_text(clip=fitz.Rect(0,art.y1+1,W,H)).split())
assert after_caption==caption,(after_caption,caption)
result={'status':'pass','source':str(source.relative_to(B)),'source_sha256':sha(source),'output':'submission_figures/Figure01.pdf','sha256':sha(B/'submission_figures/Figure01.pdf'),'page_size_pt':[W,H],'page_size_mm':[297,210],'orientation':'landscape','pages':1,'artwork_scale_relative_to_portrait':scale,'artwork_crop_pt':list(clip),'artwork_target_pt':list(art),'caption_lines':len(lines),'caption_font_pt':size,'caption_text_preserved':True,'artwork_word_multiset_preserved':True,'artwork_word_count':sum(srcwords.values()),'artwork_rasterized':False,'scientific_data_changed':False,'combined_pages':len(new),'combined_pages_pixel_identical_to_previous':unchanged,'combined_first_page_pixel_identical':True,'preview':'Reports/figure1_landscape_20261019/Figure01_landscape_A4_preview.png'}
dump(R/'validation.json',result)
result['outputs_sha256']={rel:sha(B/rel) for rel in targets[:5]}
dump(R/'validation.json',result)
status=json.loads((B/'release_status.json').read_text());status.update(scientific_artwork_changed=False,figure_layout_changed=True,figure1_layout='landscape_A4',figure1_validation='Reports/figure1_landscape_20261019/validation.json');dump(B/'release_status.json',status)
print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
