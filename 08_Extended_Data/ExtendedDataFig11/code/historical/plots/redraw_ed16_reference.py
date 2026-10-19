"""ED16a: fixed UMAP and a separate Fig. 4a-style reference schematic."""
from pathlib import Path
import sys, json, hashlib
B = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(B/'Code'), str(B/'Code/vendor')]
import pymupdf as fitz
import openpyxl
FONTS = Path('/usr/share/fonts/truetype/msttcorefonts')
def rgb(s): return tuple(int(s[i:i+2],16)/255 for i in (1,3,5))

def render(include_completion=True):
    from ed16_completion import ready
    if include_completion and ready():
        from plots.complete_ed16_fields import render as complete
        return complete()
    from style_extended_typography import style_one, SOURCE
    stem='Extended_Data_Fig16_overview'
    labels=json.loads((B/'Code/inputs/typography_labels.json').read_text())
    style_one(SOURCE/(stem+'.pdf'),labels.get(stem+'.pdf',[]))
    target=B/'Assets/arial'/(stem+'.pdf');original=fitz.open(target)
    source=B/'Code/inputs/extended_v2/Extended_Data_Fig18.xlsx'
    w=openpyxl.load_workbook(source,read_only=True,data_only=True)
    cells=list(w['Cells and coordinates'].values)[1:]
    edges=list(w['Direct evidence edges'].values)[1:]
    palette={r[0]:r[2] for r in list(w['Cell-type palette'].values)[1:]}
    assert len(cells)==754 and len(edges)==4 and len(palette)==8
    d=fitz.open();p=d.new_page(width=595.276,height=841.89)
    clip=fitz.Rect(0,372,595.276,841.89)
    p.show_pdf_page(clip,original,0,clip=clip)
    def text(x,y,s,size=9,bold=False,center=False,rotate=0):
        ff=FONTS/('Arial_Bold.ttf' if bold else 'Arial.ttf')
        if center:x-=fitz.Font(fontfile=str(ff)).text_length(s,fontsize=size)/2
        p.insert_text((x,y),s,fontname='ED16Bold' if bold else 'ED16Regular',fontfile=str(ff),fontsize=size,rotate=rotate)
    def arrow(start,end):
        start,end=fitz.Point(start),fitz.Point(end)
        v=(end-start)/abs(end-start);normal=fitz.Point(-v.y,v.x)
        p.draw_line(start,end,color=(.12,.12,.12),width=1)
        shape=p.new_shape();shape.draw_polyline([end-v*4.5+normal*2,end,end-v*4.5-normal*2])
        shape.finish(color=(.12,.12,.12),width=.8,closePath=False);shape.commit()
    text(24,20,'EXTENDED DATA FIGURE 16',8,True)
    text(58,50,'a',8,True)
    text(175,50,'KPTracer UMAP',10,True,True)
    text(175,65,'754 cells | 8 cell states',8,center=True)
    text(431,50,'Differentiation schematic',10,True,True)
    text(431,65,'Four reference relations',8,center=True)
    xs=[r[1] for r in cells];ys=[r[2] for r in cells]
    scale=min(205/(max(xs)-min(xs)),228/(max(ys)-min(ys)))
    cx=(min(xs)+max(xs))/2;cy=(min(ys)+max(ys))/2
    for _,x,y,cat in cells:
        p.draw_circle((178+(x-cx)*scale,190-(y-cy)*scale),1.05,color=None,fill=rgb(palette[cat]),fill_opacity=.9)
    arrow((68,307),(100,307));arrow((68,307),(68,275))
    text(84,320,'UMAP1',8,center=True);text(60,309,'UMAP2',8,rotate=90)
    # Fig. 4a-style colored circular nodes and labels; exactly four source edges.
    nodes={'AT2-like':(421,106),'Gastric-like':(355,165),'Lung progenitor-like':(486,165),
           'Endoderm-like':(393,227),'Pre-EMT':(393,266),'Early EMT-1':(393,305)}
    for src,dst in edges:
        a,b=fitz.Point(nodes[src]),fitz.Point(nodes[dst]);v=(b-a)/abs(b-a)
        arrow(a+v*7,b-v*7)
    for name,point in nodes.items():p.draw_circle(point,3.5,color=None,fill=rgb(palette[name]))
    text(431,109,'AT2-like');text(355,183,'Gastric-like',center=True)
    text(486,183,'Lung progenitor-like',center=True)
    for name in ['Endoderm-like','Pre-EMT','Early EMT-1']:
        x,y=nodes[name];text(x+11,y+3,name)
    for i,(name,color) in enumerate(palette.items()):
        col,row=divmod(i,3);x=72+163*col;y=339+11*row
        p.draw_circle((x,y-2.5),2,color=None,fill=rgb(color));text(x+7,y,name,7.4)
    d.subset_fonts();assert len(p.get_images())==0
    for block in p.get_text('dict')['blocks']:
        for line in block.get('lines',[]):
            for span in line['spans']:assert p.rect.contains(fitz.Rect(span['bbox'])),span['text']
    assert p.get_pixmap(clip=clip).samples==original[0].get_pixmap(clip=clip).samples
    tmp=target.with_suffix('.tmp.pdf');d.save(tmp,garbage=4,deflate=True)
    p.get_pixmap(matrix=fitz.Matrix(2.5,2.5)).save(target.with_suffix('.png'))
    target.with_suffix('.svg').write_text(p.get_svg_image(text_as_path=False))
    original.close();d.close();tmp.replace(target)
    report={'figure':16,'cells':len(cells),'cell_types':len(palette),'reference_edges':edges,
            'diagram_nodes':nodes,'umap_transform':{'scale':scale,'center_data':[cx,cy],'center_page':[178,190]},
            'umap_axis_labels':['UMAP1','UMAP2'],'umap_centroid_arrows':False,
            'schematic_style_reference':'Figures/Figure04.pdf, panel a','panel_b_pixel_identical':True,
            'embedded_raster_images':0,'source':str(source.relative_to(B)),
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'model_rerun':False,'scores_recomputed':False}
    (B/'Provenance/ed16_reference_revision.json').write_text(json.dumps(report,indent=2))
    print('ED16 UMAP axes and separate schematic rendered; GDC panel preserved.')

if __name__=='__main__':render()
