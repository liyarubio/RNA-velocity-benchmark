"""One-page ED16 using the ED15 layout and preserved vector method fields."""
from pathlib import Path
import sys,json,hashlib,csv,math,subprocess
B=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(B/'Code'),str(B/'Code/vendor')]
import pymupdf as fitz
import openpyxl
from ed16_completion import results

W,H=595.276,841.89
SOURCE=B/'Code/inputs/ed16_single_page/method_fields_before_layout.pdf'
FONT=Path('/usr/share/fonts/truetype/msttcorefonts')
COLORS=['#D62728','#CAB2D6','#6A3D9A','#33A02C','#FFB547','#B2DF8A','#E6CF5C','#0000CD','#FFB6C1','#DB7093','#20B2AA','#87CEFA','#7B4173','#9ACD32','#FF8C00','#1E90FF','#FF6347','#4682B4','#CD5C5C']
def rgb(s):return tuple(int(s[i:i+2],16)/255 for i in (1,3,5))

def render():
    w=openpyxl.load_workbook(B/'Code/inputs/extended_v2/Extended_Data_Fig18.xlsx',read_only=True,data_only=True)
    cells=list(w['Cells and coordinates'].values)[1:]
    palette={r[0]:r[2] for r in list(w['Cell-type palette'].values)[1:]}
    edges=list(w['Direct evidence edges'].values)[1:]
    scores=dict(list(w['GDC direct edges'].values)[1:])
    for method,r in results().items():scores[method]=r['GDC']
    original_order=list(scores);order=sorted(scores,key=lambda m:-scores[m]);colors=dict(zip(original_order,COLORS))
    doc=fitz.open();p=doc.new_page(width=W,height=H)
    def text(x,y,s,size=7,bold=False,center=False,right=False,color=(0,0,0),rotate=0,measure_as=None):
        ff=FONT/('Arial_Bold.ttf' if bold else 'Arial.ttf')
        length=fitz.Font(fontfile=str(ff)).text_length(measure_as or s,fontsize=size)
        if center:x-=length/2
        if right:x-=length
        p.insert_text((x,y),s,fontsize=size,fontname='KPTitle' if bold else 'KPText',fontfile=str(ff),color=color,rotate=rotate)
    def arrow(a,b):
        a,b=fitz.Point(a),fitz.Point(b);v=(b-a)/abs(b-a);q=b-v*3.5;n=fitz.Point(-v.y,v.x)*1.5
        p.draw_line(a,b,color=(.12,.12,.12),width=.65)
        p.draw_polyline([q+n,b,q-n],color=(.12,.12,.12),width=.65)
    text(24,20,'EXTENDED DATA FIGURE 16',8,True)
    text(30,48,'a',9,True);text(46,48,'KPTracer distribution and reference',9,True)
    text(108,66,'UMAP: 754 cells',7,center=True)
    text(258,66,'Differentiation schematic',7,center=True)
    xs=[r[1] for r in cells];ys=[r[2] for r in cells]
    scale=min(124/(max(xs)-min(xs)),190/(max(ys)-min(ys)))
    cx=(max(xs)+min(xs))/2;cy=(max(ys)+min(ys))/2
    px,py=108,177
    for _,x,y,cat in cells:p.draw_circle((px+(x-cx)*scale,py-(y-cy)*scale),1.10,color=None,fill=rgb(palette[cat]))
    arrow((42,272),(42,246));arrow((42,272),(67,272))
    text(54.5,283,'UMAP1',6.5,center=True);text(34,273,'UMAP2',6.5,rotate=90)
    from plots.ed16_literature_schematic import draw as draw_literature
    nodes,literature=draw_literature(p,text,palette,edges)
    for i,(name,color) in enumerate(palette.items()):
        col,row=divmod(i,4);x=42+col*153;y=298+row*9
        p.draw_circle((x,y-2),1.6,color=None,fill=rgb(color));text(x+6,y,name,6.3)
    text(359,48,'b',9,True);text(375,48,'GDC',9,True)
    x0,x1=455,558;zero=(x0+x1)/2
    bars=[]
    for i,m in enumerate(order):
        y=70+i*12
        # Retain the published alignment after the prior asterisk removal.
        text(x0-6,y+2.4,m,6.6,right=True,measure_as='TopicVelo*' if m=='TopicVelo' else None)
        end=zero+scores[m]*(x1-x0)/2
        rect=fitz.Rect(min(zero,end),y-4,max(zero,end),y+4)
        p.draw_rect(rect,color=None,fill=rgb(colors[m]));bars.append({'method':m,'GDC':scores[m],'rect':list(rect)})
    p.draw_line((zero,59),(zero,301),color=(.4,.4,.4),width=.45)
    p.draw_line((x0,301),(x1,301),color=(0,0,0),width=.5)
    for x,value in [(x0,'-1'),(zero,'0'),(x1,'1')]:
        p.draw_line((x,301),(x,304),color=(0,0,0),width=.5);text(x,312,value,6.5,center=True)
    text(362.8628,325,'Native-transition displacement proxy',5.8)
    text(30,345,'c',9,True);text(46,345,'Method-specific velocity fields',9,True)
    # Keep a clean header for exact coordinate/score validation, before adding fields.
    header=B/'Reports/ed16_single_page_header.pdf';doc.save(header,garbage=4,deflate=True)
    source=fitz.open(SOURCE);placements=[]
    from plots.ed16_aligned_fields import prepare,render as align_field
    prepared=prepare(source[0],cells,palette)
    field_dir=B/'Reports/ed16_aligned_fields';field_dir.mkdir(exist_ok=True)
    for i,m in enumerate(order):
        j=original_order.index(m);sr,sc=divmod(j,4)
        clip=fitz.Rect(51+126.2*sc,52+114.76*sr,173+126.2*sc,152+114.76*sr)
        if m in ('TFvelo','TopicVelo'):clip=fitz.Rect(51+126*sc,392,173+126*sc,498)
        row,col=divmod(i,5);x=36+col*111;y=357+row*112
        text(x+52.5,y+7,m,7.2,True,center=True,measure_as='TopicVelo*' if m=='TopicVelo' else None)
        frame=fitz.Rect(x,y+14,x+105,y+112)
        aligned,details=align_field(clip,cells,prepared)
        panel_path=field_dir/(str(i+1).zfill(2)+'.pdf')
        aligned.save(panel_path,garbage=4,deflate=True)
        p.show_pdf_page(frame,aligned,0)
        placements.append({'method':m,'source_index':j,'source_clip':list(clip),'target_frame':list(frame),
                           'aligned_panel':str(panel_path.relative_to(B)),**details})
        aligned.close()
    dest=B/'Assets/arial/Extended_Data_Fig16.pdf';doc.subset_fonts();doc.save(dest,garbage=4,deflate=True)
    p.get_pixmap(matrix=fitz.Matrix(2.5,2.5)).save(dest.with_suffix('.png'))
    subprocess.run(['pdftocairo','-svg',str(dest),str(dest.with_suffix('.svg'))],check=True)
    assert not p.get_images() and '<image' not in dest.with_suffix('.svg').read_text()
    record={'pages':1,'layout_reference':'Extended Data Fig. 15','cells':754,'umap_point_radius_pt':1.10,'reference_edges':edges,
            'umap_transform':{'scale':scale,'center_data':[cx,cy],'center_page':[px,py]},'diagram_nodes':nodes,'literature_schematic':literature,
            'GDC_display_order':order,'field_display_order':order,'bars':bars,'field_placements':placements,
            'field_source':str(SOURCE.relative_to(B)),'field_source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'fields_preserved_as_affine_vector_paths':19,'field_cell_bounds_identical':True,'scores_unchanged':True,'models_rerun':False,'embedded_raster_images':0}
    (B/'Provenance/ed16_single_page.json').write_text(json.dumps(record,indent=2))
    with (B/'Source_data/kptracer_completion/GDC_display_order.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['rank','method','GDC'])
        writer.writerows((i+1,m,scores[m]) for i,m in enumerate(order))
    print(json.dumps({'pages':1,'methods':19,'order':order,'vector_only':True},indent=2))

if __name__=='__main__':render()
