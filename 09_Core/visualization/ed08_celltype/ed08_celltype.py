"""Recolor ED8a by the exact ED8c cell types; preserve point geometry and panels b-d."""
from pathlib import Path
import argparse,hashlib,json,sys
import numpy as np
import pandas as pd
import pymupdf as fitz
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'mean_top3'))
from mean_top3 import visible_drawings
PALETTE={'Baso/Eos':'#279e68','Ery/Meg':'#d62728','Lym/pDC':'#1f77b4','Mono':'#aa40fc','Neutro':'#aec7e8','Progenitors':'#ff7f0e','Unknown':'#808080'}
EXTRA_PALETTE={'Other':'#b0b0b0'}
DISPLAY_LABELS={'Mast':'Other','Ccr7_DC':'Other'}
DISPLAY_PALETTE={**PALETTE,**EXTRA_PALETTE}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def root_data():
    for root in HERE.parents:
        if (root/'12_Provenance/current_sources.json').exists():return root/'08_Extended_Data/ExtendedDataFig08/data/reference_restoration'
        if (root/'source_data/Extended_Data_Fig14.xlsx').exists():return root/'data/reference_restoration'
    raise RuntimeError('Cannot locate saved inputs')
def read_data():
    d=root_data();full=pd.read_csv(d/'scLT_timepoints.csv');labels=pd.read_csv(d/'scLT_cell_types.csv');ref=pd.read_csv(d/'scLT_reference_cells.csv')
    full=full.merge(labels,on='cell_id',how='left',validate='one_to_one',sort=False)
    assert len(full)==3186 and full.cell_type.notna().all() and set(full.cell_type)==set(PALETTE)|set(DISPLAY_LABELS)
    assert full.set_index('cell_id').loc[ref.cell_id,'cell_type'].tolist()==ref.cell_type.tolist()
    assert np.allclose(full.set_index('cell_id').loc[ref.cell_id,['x','y']],ref[['x','y']])
    return full,{n:sha(d/n) for n in ['scLT_timepoints.csv','scLT_cell_types.csv','scLT_reference_cells.csv']}
def recolor(source=None):
    original=source if isinstance(source,fitz.Document) else fitz.open(source or HERE/'assets/Extended_Data_Fig08_before.pdf');doc=fitz.open();doc.insert_pdf(original,from_page=0,to_page=0);page=doc[0]
    data,hashes=read_data();draws=list(visible_drawings(original[0]))
    points=[d for d in draws if d.get('fill') and abs(d['fill_opacity']-.8)<1e-6 and 1.5<d['rect'].width<2 and d['rect'].y1<230 and d['rect'].x1<310]
    assert len(points)==3186,len(points)
    axes=[d['rect'] for d in draws if d.get('fill')==(1.,1.,1.) and 100<d['rect'].width<115 and 110<d['rect'].height<130 and d['rect'].x1<310 and d['rect'].y1<230]
    axes=sorted(axes,key=lambda r:r.x0);assert len(axes)==2,axes
    changes=[fitz.Rect(r.x0+.3,r.y0+.3,r.x1-.3,r.y1-.3) for r in axes]
    batches=[];max_error=0
    for j,samples in enumerate([['d2_5'],['d5r1','d5r2']]):
        q=pd.concat([data[data['sample']==s] for s in samples],ignore_index=True)
        paths=[p for p in points if axes[j].contains(p['rect'])];assert len(paths)==len(q)
        centers=np.array([[(p['rect'].x0+p['rect'].x1)/2,(p['rect'].y0+p['rect'].y1)/2] for p in paths])
        xy=np.c_[q[['x','y']].to_numpy(),np.ones(len(q))];transform=np.linalg.lstsq(xy,centers,rcond=None)[0]
        err=float(np.max(np.abs(xy@transform-centers)));assert err<.002,err;max_error=max(max_error,err)
        batches.append((q,paths))
    # Remove only the three literal UTF-16 legend strings in the original
    # Matplotlib form. Avoid rewriting the other panels' complex PDF paths.
    removed={label:0 for label in ['d2_5','d5r1','d5r2']}
    for xref in {x[0] for x in page.get_xobjects()}:
        stream=doc.xref_stream(xref)
        if 'Day 2.5'.encode('utf-16-be') not in stream:continue
        for label in removed:
            token=b'[ ('+label.encode('utf-16-be')+b') ] TJ'
            removed[label]+=stream.count(token);stream=stream.replace(token,b'[] TJ')
        doc.update_stream(xref,stream)
    assert all(n==1 for n in removed.values()),removed
    for region in changes:page.draw_rect(region,color=None,fill=(1,1,1))
    shape=page.new_shape();expected=[]
    for q,paths in batches:
        for label,d in zip(q.cell_type,paths):
            color=tuple(int(DISPLAY_PALETTE[DISPLAY_LABELS.get(label,label)][i:i+2],16)/255 for i in (1,3,5));expected.append(color)
            for item in d['items']:
                assert item[0]=='c';shape.draw_bezier(*item[1:])
            shape.finish(color=None,fill=color,fill_opacity=d['fill_opacity'],closePath=d['closePath'],even_odd=d['even_odd'])
    shape.commit()
    scale=axes[0].width/110
    legend=fitz.Rect(axes[0].x0,axes[0].y1+15*scale,axes[1].x1,axes[0].y1+43*scale)
    page.insert_font(fontname='helv')
    for i,(label,color) in enumerate(DISPLAY_PALETTE.items()):
        x=legend.x0+(5+(i%3)*80)*scale;y=legend.y0+(6+(i//3)*8)*scale
        rgb=tuple(int(color[k:k+2],16)/255 for k in (1,3,5))
        page.draw_circle((x,y-1.3*scale),1.4*scale,color=None,fill=rgb)
        page.insert_text((x+4*scale,y),label,fontname='helv',fontsize=5.3*scale)
    changes.append(legend)
    saved=fitz.open(stream=doc.tobytes(),filetype='pdf');actual=[d for d in visible_drawings(saved[0]) if d.get('fill') and abs(d['fill_opacity']-.8)<1e-6 and 1.5<d['rect'].width<2 and d['rect'].y1<230 and d['rect'].x1<310]
    # The opaque white panel interiors cover the retained original vector dots.
    # Validate the final 3186 visible replacement dots in their drawing order.
    actual=actual[-3186:]
    assert len(actual)==3186
    assert np.allclose([d['fill'] for d in actual],expected,atol=1e-6)
    assert np.allclose([list(d['rect']) for d in actual],[list(d['rect']) for _,paths in batches for d in paths],atol=.002)
    for label in ['d2_5','d5r1','d5r2']:assert label not in saved[0].get_text(clip=fitz.Rect(axes[0].x0,axes[0].y0,axes[1].x1,axes[1].y1))
    zoom=3;pix=[d[0].get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False) for d in [original,saved]]
    a,b=[np.frombuffer(p.samples,dtype=np.uint8).reshape(p.height,p.width,3) for p in pix];mask=np.ones(a.shape[:2],bool)
    for r in changes:mask[int(r.y0*zoom)-2:int(r.y1*zoom)+3,int(r.x0*zoom)-2:int(r.x1*zoom)+3]=False
    delta=np.abs(a.astype(np.int16)-b.astype(np.int16))[mask];assert delta.max()==0
    return doc,dict(status='pass',panel='ED08a',cells=3186,reference_labels_matched=1756,cell_type_palette=PALETTE,additional_palette=EXTRA_PALETTE,display_label_mapping=DISPLAY_LABELS,input_sha256=hashes,
        point_mapping_max_error_pt=max_error,changed_regions=[list(r) for r in changes],colors_verified=3186,point_geometry_preserved=True,
        point_opacity=.8,display_counts=data.cell_type.replace(DISPLAY_LABELS).value_counts().to_dict(),sample_counts=data['sample'].value_counts().to_dict(),unchanged_area_max_pixel_difference=int(delta.max()),
        original_sample_legends_removed=True,legend='Panel a shares seven colors with panel c; Mast and Ccr7_DC are displayed in grey as Other, with original annotations retained',models_rerun=False,numerical_inputs_changed=False)
def save(doc,stem,clip=None):
    stem=Path(stem);stem.parent.mkdir(parents=True,exist_ok=True)
    if clip:
        d=fitz.open();p=d.new_page(width=clip.width,height=clip.height);p.show_pdf_page(p.rect,doc,0,clip=clip);doc=d
    doc.save(stem.with_suffix('.pdf'),garbage=4,deflate=True)
    with fitz.open(stem.with_suffix('.pdf')) as saved:
        saved[0].get_pixmap(matrix=fitz.Matrix(3.5 if clip else 2.5,3.5 if clip else 2.5),alpha=False).save(stem.with_suffix('.png'))
        stem.with_suffix('.svg').write_text(saved[0].get_svg_image(text_as_path=False))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=Path('demo_output/ed08a'));p.add_argument('--full-figure',action='store_true');args=p.parse_args()
    doc,audit=recolor();save(doc,args.out/('Extended_Data_Fig08' if args.full_figure else 'ED08a'),None if args.full_figure else fitz.Rect(40,60,300,245))
    (args.out/'validation.json').write_text(json.dumps(audit,indent=2)+'\n');print('ED08a: 3186 cell colors and coordinates verified; other panels unchanged.')
if __name__=='__main__':main()
