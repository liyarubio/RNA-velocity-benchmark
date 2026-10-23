"""Order formal Figure 3h method columns by descending available-case mean GDC.

The retained vector template supplies the exact cells and colors. Current CSV
values determine the order; missing results remain missing. This is a
presentation change, not a model rerun. The best-performance envelope stays last.
"""
from pathlib import Path
import argparse, hashlib, json, sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'mean_top3'))
from mean_top3 import fitz, np, pd, visible_drawings, text_spans, original_font, save_artwork
HERE=Path(__file__).resolve().parent
ORDER=['velocyto','scVelo(sto)','scVelo(dyn)','VeloAE','VeloVAE','dynamo','UniTVelo',
       'DeepVelo(VAE)','DeepVelo(GCN)','LatentVelo','cellDancer','veloVI','STT',
       'TFvelo','biVI','DeepKINET','TopicVelo','SDEvelo','cell2fate']
CLIP=fitz.Rect(280,590,558,749)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read_scores():
    for root in HERE.parents:
        if (root/'12_Provenance/current_sources.json').is_file():
            path=root/'03_Figure03/data/figure3_detail/Figure3h_dataset10_global_performance.csv';break
        if (root/'source_data/Figure03.xlsx').is_file():
            path=root/'data/figure3_detail/Figure3h_dataset10_global_performance.csv';break
    else:raise RuntimeError('Current data root not found')
    expected=json.loads((HERE/'assets/input_hashes.json').read_text())['csv_sha256']
    assert sha(path)==expected,'Numerical inputs changed: refresh the complete heatmap template first'
    raw=pd.read_csv(path,index_col=0).reindex(ORDER)
    data=raw.drop(columns=['Mean','Median']);assert data.shape==(19,10)
    means=data.mean(axis=1,skipna=True)
    assert np.allclose(means,raw['Mean'],atol=1e-6,rtol=0)
    order=means.sort_values(ascending=False,kind='stable',na_position='last').index.tolist()
    return path,data,means,order

def reorder(doc=None):
    if doc is None:doc=fitz.open(HERE/'assets/Figure03_before.pdf')
    original=doc;out=fitz.open();out.insert_pdf(original);page=out[0]
    path,data,means,order=read_scores();draws=list(visible_drawings(original[0]))
    strips=sorted([d for d in draws if d.get('fill') and 375<d['rect'].x0<535
                   and 709<d['rect'].y0<711 and 3<d['rect'].height<4],key=lambda d:d['rect'].x0)
    assert len(strips)==19,len(strips)
    spans={s['text']:(s,d) for s,d in text_spans(original[0])
           if s['text'] in ORDER+['Best performance'] and s['bbox'][0]>370 and s['bbox'][1]>710}
    assert set(spans)==set(ORDER+['Best performance'])
    # Infer each method's currently plotted column from its right-aligned 50-degree label.
    # The template must still have the original order (double sorting would be ambiguous).
    endpoints=[]
    for method in ORDER:
        s,direction=spans[method];font=original_font(original,original[0],s['font'],required=method.replace(" ",""))
        endpoints.append(s['origin'][0]+font.text_length(method,fontsize=s['size'])*direction[0])
    assert all(a<b for a,b in zip(endpoints,endpoints[1:])), 'Template is already reordered'
    cells=[d for d in draws if d.get('fill') and 378<d['rect'].x0<532
           and 618<d['rect'].y0<708 and 7<d['rect'].width<9 and 8<d['rect'].height<10]
    assert len(cells)==int(data.notna().sum().sum()),len(cells)
    y0=min(d['rect'].y0 for d in cells);y1=max(d['rect'].y1 for d in strips)
    matrix=fitz.Rect(strips[0]['rect'].x0-.2,y0-.2,strips[-1]['rect'].x1+.2,y1+.2)
    labels=fitz.Rect(min(s['bbox'][0] for s,_ in spans.values())-.3,
                     min(s['bbox'][1] for s,_ in spans.values())-.3,
                     max(s['bbox'][2] for s,_ in spans.values())+.3,
                     max(s['bbox'][3] for s,_ in spans.values())+.3)
    for region in [matrix,labels]:page.add_redact_annot(region,fill=(1,1,1),cross_out=False)
    page.apply_redactions(images=0,graphics=1,text=0)
    # Isolate matrix paths before importing columns. Reusing the full source
    # page 19 times would inflate the exported SVG with hidden copies of a-g.
    vectors=fitz.open();layer=vectors.new_page(width=page.rect.width,height=page.rect.height)
    for drawing in draws:
        if not matrix.contains(drawing['rect']):continue
        shape=layer.new_shape()
        for item in drawing['items']:
            if item[0]=='re':shape.draw_rect(item[1])
            elif item[0]=='l':shape.draw_line(item[1],item[2])
            else:raise AssertionError('Unexpected heatmap path: '+item[0])
        shape.finish(color=drawing.get('color'),fill=drawing.get('fill'),
                     width=drawing.get('width') or 0,
                     closePath=bool(drawing.get('closePath')),
                     even_odd=bool(drawing.get('even_odd')))
        shape.commit()
    # Permute whole vector cell columns together with their method-color strips.
    for new_index,method in enumerate(order):
        old_index=ORDER.index(method);old=strips[old_index]['rect'];new=strips[new_index]['rect']
        clip=fitz.Rect(old.x0,y0-.15,old.x1,y1+.15)
        dest=fitz.Rect(new.x0,y0-.15,new.x1,y1+.15)
        page.show_pdf_page(dest,vectors,0,clip=clip,keep_proportion=False)
    label_origins={}
    for method in order+['Best performance']:
        span,direction=spans[method];font=original_font(original,original[0],span['font'],required=method.replace(" ",""))
        name='SortedH'+str(order.index(method) if method in order else 19)
        page.insert_font(fontname=name,fontbuffer=font.buffer)
        dx=0
        if method in order:
            old=strips[ORDER.index(method)]['rect'];new=strips[order.index(method)]['rect']
            dx=(new.x0+new.x1-old.x0-old.x1)/2
        origin=fitz.Point(span['origin'][0]+dx,span['origin'][1])
        angle=np.degrees(np.arctan2(-direction[1],direction[0]))
        if method == 'Best performance':
            # The original bold subset omits a space glyph; preserve the PDF's
            # explicit inter-word positioning instead of substituting a font.
            raw=next(s for b in original[0].get_text('rawdict')['blocks']
                     for line in b.get('lines',[]) for s in line['spans']
                     if ''.join(c['c'] for c in s['chars'])==method and s['bbox'][0]>370)
            for word,index in [('Best',0),('performance',5)]:
                at=fitz.Point(raw['chars'][index]['origin'])
                page.insert_text(at,word,fontname=name,fontsize=span['size'],morph=(at,fitz.Matrix(angle)))
        else:
            page.insert_text(origin,method,fontname=name,fontsize=span['size'],morph=(origin,fitz.Matrix(angle)))
        label_origins[method]=list(origin)
    assert sorted(s['text'] for s,_ in text_spans(page) if s['text'] in spans and s['bbox'][0]>350 and s['bbox'][1]>710)==sorted(spans)
    zoom=3;a=original[0].get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False);b=page.get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False)
    before=np.frombuffer(a.samples,dtype=np.uint8).reshape(a.height,a.width,3)
    after=np.frombuffer(b.samples,dtype=np.uint8).reshape(b.height,b.width,3)
    mask=np.ones(before.shape[:2],dtype=bool)
    for r in [matrix,labels]:mask[int(r.y0*zoom)-2:int(r.y1*zoom)+3,int(r.x0*zoom)-2:int(r.x1*zoom)+3]=False
    delta=np.abs(before.astype(np.int16)-after.astype(np.int16))[mask]
    assert np.max(delta)==0,'Pixels changed outside panel h method columns/labels'
    # Verify the reordered vector fills against the original column-by-column colors.
    # Each imported column carries its original colors; pixel centers also check clipping/order.
    source_pix=original[0].get_pixmap(matrix=fitz.Matrix(4,4),alpha=False)
    result_pix=page.get_pixmap(matrix=fitz.Matrix(4,4),alpha=False)
    aa=np.frombuffer(source_pix.samples,dtype=np.uint8).reshape(source_pix.height,source_pix.width,3)
    bb=np.frombuffer(result_pix.samples,dtype=np.uint8).reshape(result_pix.height,result_pix.width,3)
    ys=sorted(set(round((d['rect'].y0+d['rect'].y1)/2,2) for d in cells))
    assert len(ys)==10
    for j,method in enumerate(order):
        old=strips[ORDER.index(method)]['rect'];new=strips[j]['rect']
        for y in ys+[(y1+strips[0]['rect'].y0)/2]:
            assert np.array_equal(aa[int(y*4),int((old.x0+old.x1)*2)],bb[int(y*4),int((new.x0+new.x1)*2)])
    audit=dict(status='pass',panel='Figure03h',input_sha256=sha(path),order=order,
        means=means.loc[order].to_dict(),valid_datasets=data.count(axis=1).loc[order].to_dict(),
        missing_policy='Mean of available values; no zero imputation',ties='Stable previous display order',
        best_performance='Unchanged per-dataset envelope, excluded from ranking and retained at the right',
        numeric_inputs_changed=False,models_rerun=False,changed_regions=[list(matrix),list(labels)],
        unchanged_area_max_pixel_difference=int(delta.max()),permuted_cell_centers_checked=190,
        method_color_strips_checked=19,missing_cells_retained=int(data.isna().sum().sum()),label_origins=label_origins)
    return out,audit

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,default=Path('demo_output/figure3h_mean'))
    parser.add_argument('--full-figure',action='store_true');args=parser.parse_args()
    doc,audit=reorder();name='Figure03' if args.full_figure else 'Figure03h'
    save_artwork(doc,args.out/name,clip=None if args.full_figure else CLIP)
    (args.out/'validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    pd.DataFrame({'Method':audit['order'],'Mean':[audit['means'][m] for m in audit['order']],
                  'Valid datasets':[audit['valid_datasets'][m] for m in audit['order']]}).to_csv(args.out/'method_order.csv',index=False)
    print('Figure03h: 19 methods sorted by descending mean; vector/color and unchanged-area checks passed')

if __name__=='__main__':main()
