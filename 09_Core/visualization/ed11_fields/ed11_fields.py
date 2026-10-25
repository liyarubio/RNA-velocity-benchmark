"""Align KPTracer field widths and arrange methods by descending displayed GDC.

Native vector paths are transformed with their 754 fixed cell coordinates.
PhyloVelo retains backward arrows and its non-comparable adapted-score status.
"""
from pathlib import Path
import argparse,hashlib,json,sys,subprocess,tempfile
import numpy as np,pandas as pd,pymupdf as fitz
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'mean_top3'))
from mean_top3 import visible_drawings,original_font,text_spans

def root_paths():
    for root in HERE.parents:
        if (root/'12_Provenance/current_sources.json').is_file():
            return root,root/'10_Shared/data',root/'08_Extended_Data/ExtendedDataFig11/code/historical/plots'
        if (root/'source_data/Extended_Data_Fig17.xlsx').is_file():return root,root/'data',root/'upstream/plots'
    raise RuntimeError('Cannot locate saved data')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_inputs():
    root,data,helper=root_paths();sys.path.insert(0,str(helper))
    paths={'scores':data/'figure5_1021/KPTracer_GDC_with_PhyloVelo.csv','cells':data/'kptracer_phylovelo_20261021/results/phylovelo_directions.csv','palette':data/'kptracer_phylovelo_20261021/inputs/cell_type_palette.csv'}
    scores=pd.read_csv(paths['scores']);cells=pd.read_csv(paths['cells']);pal=pd.read_csv(paths['palette'])
    assert len(scores)==20 and scores.method.is_unique and len(cells)==754 and cells.cell_id.is_unique
    assert pd.isna(scores.set_index('method').loc['PhyloVelo','rank'])
    order=scores.sort_values('GDC',ascending=False,kind='stable').method.tolist()
    return scores,cells,dict(zip(pal.cell_type,pal.color)),order,{k:sha(p) for k,p in paths.items()}
def render(source=None):
    original=source if isinstance(source,fitz.Document) else fitz.open(source or HERE/'assets/Extended_Data_Fig11_before.pdf')
    p=original[0];scores,cells,palette,order,hashes=read_inputs()
    from ed16_aligned_fields import render as align_field,CELL_BOUNDS,rgb
    spans={s['text'].replace('\u00a0',' '):(s,d) for s,d in text_spans(p) if s['bbox'][1]>330}
    titles={m:spans[m if m!='PhyloVelo' else 'PhyloVelo (backward)'][0] for m in order}
    old_order=sorted(order,key=lambda m:(round(titles[m]['bbox'][1]/10),titles[m]['bbox'][0]))
    assert old_order[-1]=='PhyloVelo'
    old_centers=[(titles[m]['bbox'][0]+titles[m]['bbox'][2])/2 for m in old_order[:5]]
    step=float(np.median(np.diff(old_centers)));row_step=titles[old_order[5]]['bbox'][1]-titles[old_order[0]]['bbox'][1]
    x_centers=[old_centers[0]+j*step for j in range(5)]
    clips={m:fitz.Rect(x_centers[i%5]-step/2,titles[m]['bbox'][3]+.5,x_centers[i%5]+step/2,titles[m]['bbox'][3]+row_step-1) for i,m in enumerate(old_order)}
    top=min(r.y0 for r in clips.values())-15;bottom=max(r.y1 for r in clips.values())+1
    area=fitz.Rect(30,top,566,min(bottom,p.rect.height-8))
    drawings=[d for d in visible_drawings(p) if area.contains(d['rect'])]
    centers=np.array([[(d['rect'].x0+d['rect'].x1)/2,(d['rect'].y0+d['rect'].y1)/2] for d in drawings])
    colors={k:rgb(v) for k,v in palette.items()};fills=np.array([d.get('fill') or (-2,-2,-2) for d in drawings])
    small=np.array([d['rect'].width<4 and d['rect'].height<4 for d in drawings])
    is_cell=small&np.any(np.max(np.abs(fills[:,None,:]-np.array(list(colors.values()))[None,:,:]),axis=2)<.005,axis=1)
    xy=cells[['x','y']].to_numpy();lo=xy.min(0);span=np.ptp(xy,axis=0);l,t,r,b=CELL_BOUNDS
    positions=np.c_[l+(xy[:,0]-lo[0])/span[0]*(r-l),b-(xy[:,1]-lo[1])/span[1]*(b-t)]
    rows=list(cells[['cell_id','x','y','cell_type']].itertuples(index=False,name=None));prepared=drawings,centers,is_cell,positions,colors
    # The existing 19 methods define the common plotting scale and geometry.
    first=clips[old_order[0]];mask=(centers[:,0]>first.x0)&(centers[:,0]<first.x1)&(centers[:,1]>first.y0)&(centers[:,1]<first.y1)&is_cell
    assert mask.sum()==754;bounds=np.r_[centers[mask].min(0),centers[mask].max(0)]
    scale=(bounds[2]-bounds[0])/(r-l);assert abs(scale-(bounds[3]-bounds[1])/(b-t))<1e-5
    left=bounds[0]-l*scale;field_top=bounds[1]-t*scale
    doc=fitz.open();page=doc.new_page(width=p.rect.width,height=p.rect.height)
    # Import only the visible area outside b; clipping also excludes old titles
    # from extracted text. No full-page PDF redaction or metric rerun is needed.
    for rect in [fitz.Rect(0,0,p.rect.width,area.y0),fitz.Rect(0,area.y1,p.rect.width,p.rect.height)]:
        if not rect.is_empty:page.show_pdf_page(rect,original,0,clip=rect,keep_proportion=False)
    font=original_font(original,p,titles[old_order[0]]['font'],required=''.join(order)+'backward()')
    page.insert_font(fontname='AlignedFieldLabel',fontbuffer=font.buffer)
    title_size=titles[old_order[0]]['size'];title_y=titles[old_order[0]]['origin'][1]
    placements=[]
    for i,method in enumerate(order):
        tile,detail=align_field(clips[method],rows,prepared)
        x=left+(i%5)*step;y=field_top+(i//5)*row_step;frame=fitz.Rect(x,y,x+105*scale,y+98*scale)
        page.show_pdf_page(frame,tile,0,keep_proportion=False)
        label=method if method!='PhyloVelo' else 'PhyloVelo (backward)'
        words=label.split();gap=title_size*.278
        width=sum(font.text_length(w,fontsize=title_size) for w in words)+gap*(len(words)-1)
        at=x_centers[i%5]-width/2
        for word in words:
            page.insert_text((at,title_y+(i//5)*row_step),word,fontname='AlignedFieldLabel',fontsize=title_size)
            at+=font.text_length(word,fontsize=title_size)+gap
        placements.append(dict(method=method,GDC=float(scores.set_index('method').loc[method,'GDC']),destination=list(frame),**detail))
    saved=fitz.open(stream=doc.tobytes(),filetype='pdf');zoom=2
    pix=[d[0].get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False) for d in [original,saved]]
    aa,bb=[np.frombuffer(v.samples,dtype=np.uint8).reshape(v.height,v.width,3) for v in pix];outside=np.ones(aa.shape[:2],bool)
    outside[max(0,int(area.y0*zoom)-2):int(area.y1*zoom)+3,:]=False
    delta=np.abs(aa.astype(np.int16)-bb.astype(np.int16))[outside];assert delta.max()==0
    extracted=[s['text'].replace('\u00a0',' ') for s,_ in text_spans(saved[0]) if area.y0<s['bbox'][1]<area.y1]
    # Every method title occurs once and is located in its requested grid slot.
    for i,m in enumerate(order):
        expected=m if m!='PhyloVelo' else 'PhyloVelo (backward)';assert extracted.count(expected)==1,(expected,extracted)
    assert all(np.allclose(e['target_cell_bounds'],CELL_BOUNDS) for e in placements)
    audit=dict(status='pass',panel='KPTracer method fields (local ED11b; requested Ext Fig15b)',order=order,
        ordering='Descending displayed GDC; stable ties; display order is not comparable rank',input_sha256=hashes,
        PhyloVelo=dict(direction='native backward retained',adapted_GDC=1.,comparable_rank=None),placements=placements,
        common_frame_size_pt=[105*scale,98*scale],common_cell_bounds_pt=list(CELL_BOUNDS),methods=20,cells_per_method=754,
        changed_region=list(area),unchanged_area_max_pixel_difference=int(delta.max()),models_rerun=False,numerical_inputs_changed=False)
    return doc,audit

def save(doc,stem,panel=False):
    stem=Path(stem);stem.parent.mkdir(parents=True,exist_ok=True)
    if panel:
        crop=fitz.Rect(25,350,568,801);d=fitz.open();p=d.new_page(width=crop.width,height=crop.height);p.show_pdf_page(p.rect,doc,0,clip=crop);doc=d
    doc.save(stem.with_suffix('.pdf'),garbage=4,deflate=True)
    preview=doc
    if not panel:
        d=fitz.open();p=d.new_page(width=595.276,height=785);crop=fitz.Rect(0,(doc[0].rect.height-785)/2,doc[0].rect.width,(doc[0].rect.height+785)/2);p.show_pdf_page(p.rect,doc,0,clip=crop);preview=d
    preview[0].get_pixmap(matrix=fitz.Matrix(2,2),alpha=False).save(stem.with_suffix('.png'))
    # Poppler omits the clipped-away source fields from the SVG, avoiding
    # very large hidden copies while retaining vector paths and font outlines.
    with tempfile.TemporaryDirectory(prefix='ed11-svg-') as tmp:
        source=Path(tmp)/'preview.pdf';preview.save(source,garbage=4,deflate=True)
        subprocess.run(['pdftocairo','-svg',str(source),str(stem.with_suffix('.svg'))],check=True)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=Path('demo_output/ed11b'));p.add_argument('--full-figure',action='store_true');args=p.parse_args()
    doc,audit=render();save(doc,args.out/('Extended_Data_Fig11' if args.full_figure else 'ED11b'),panel=not args.full_figure)
    (args.out/'validation.json').write_text(json.dumps(audit,indent=2)+'\n');pd.DataFrame({'method':audit['order'],'displayed_GDC':[x['GDC'] for x in audit['placements']]}).to_csv(args.out/'display_order.csv',index=False)
    print('20 fields aligned to common geometry; descending displayed GDC verified; panel a unchanged.')
if __name__=='__main__':main()
