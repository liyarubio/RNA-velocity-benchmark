"""Order Figure 3h by available-case mean GDC and display a method Mean row.

The retained vector template supplies the exact cells and colors. Current CSV
values determine the order; missing results remain missing. This is a
presentation change, not a model rerun. The best-performance envelope stays last.
"""
from pathlib import Path
import argparse, hashlib, json, re, sys
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
    audit['mean_row']=add_mean_row(out,means.loc[order],strips)
    audit['changed_regions'].append(audit['mean_row']['changed_region'])
    return out,audit

def add_mean_row(doc,means,strips):
    """Make room inside panel h, retaining font sizes and the column labels.

    The ten dataset rows become slightly shorter. The Mean row uses the same
    RdYlBu_r palette and fixed [-1, 1] GDC scale, with no numeric annotations.
    Best performance is an envelope, not a method, so has no summary cell.
    """
    before=fitz.open(stream=doc.tobytes(),filetype='pdf');page=doc[0]
    region=fitz.Rect(321.8,618.4,545.7,707.7)
    top,bottom=618.623291015625,707.454345703125
    summary_height=8.5;factor=(bottom-top-summary_height)/(bottom-top)
    ymap=lambda y:top+(y-top)*factor
    point=lambda p:fitz.Point(p.x,ymap(p.y))
    drawings=[d for d in visible_drawings(before[0]) if region.contains(d['rect'])]
    labels=[s for s,_ in text_spans(before[0]) if region.contains(fitz.Rect(s['bbox']))]
    raw_labels={''.join(c['c'] for c in s['chars']):s for b in before[0].get_text('rawdict')['blocks']
                for line in b.get('lines',[]) for s in line['spans'] if region.contains(fitz.Rect(s['bbox']))}
    assert len(labels)==10
    font=original_font(before,before[0],labels[0]['font'],required='Mean0123456789.-')
    page.add_redact_annot(region,fill=(1,1,1),cross_out=False)
    page.apply_redactions(images=0,graphics=1,text=0)
    for drawing in drawings:
        shape=page.new_shape()
        for item in drawing['items']:
            if item[0]=='re':
                r=item[1];shape.draw_rect(fitz.Rect(r.x0,ymap(r.y0),r.x1,ymap(r.y1)))
            elif item[0]=='l':shape.draw_line(point(item[1]),point(item[2]))
            else:raise AssertionError('Unexpected heatmap path: '+item[0])
        shape.finish(color=drawing.get('color'),fill=drawing.get('fill'),
                     width=drawing.get('width') or 0,closePath=bool(drawing.get('closePath')),
                     even_odd=bool(drawing.get('even_odd')))
        shape.commit()
    page.insert_font(fontname='MeanRowFont',fontbuffer=font.buffer)
    for s in labels:
        center=(s['bbox'][1]+s['bbox'][3])/2
        # The embedded subset has positioned word gaps but no space glyph.
        # Preserve those original word positions instead of inserting tofu.
        chars=raw_labels[s['text']]['chars']
        for word in re.finditer(r'\S+',s['text']):
            at=chars[word.start()]['origin']
            page.insert_text((at[0],at[1]+ymap(center)-center),word.group(),
                             fontname='MeanRowFont',fontsize=s['size'])
    palette_path=HERE/'assets/gdc_palette.csv'
    palette=pd.read_csv(palette_path)
    assert palette['index'].tolist()==list(range(256))
    rgb=palette[['r','g','b']].to_numpy()/255
    color_index=lambda value:min(255,max(0,int((value+1)/2*256)))
    # Confirm that the retained palette matches every finite source cell.
    _,data,_,order=read_scores()
    source_cells=[d for d in drawings if d.get('fill') and 378<d['rect'].x0<532
                  and 7<d['rect'].width<9 and 8<d['rect'].height<10]
    assert len(source_cells)==178
    for cell in source_cells:
        r=cell['rect'];j=int(round((r.y0-top)/((bottom-top)/10)))
        i=min(range(19),key=lambda k:abs(strips[k]['rect'].x0-r.x0))
        value=data.loc[order[i]].iloc[j]
        assert np.allclose(cell['fill'],rgb[color_index(value)],atol=1e-6,rtol=0)
    center=bottom-summary_height/2+.4
    label_size=labels[0]['size']
    label_baseline=center+(font.ascender+font.descender)*label_size/2
    page.insert_text((377.8-font.text_length('Mean',fontsize=label_size),label_baseline),
                     'Mean',fontname='MeanRowFont',fontsize=label_size)
    boxes=[]
    for method,value,strip in zip(means.index,means,strips):
        r=strip['rect'];box=fitz.Rect(r.x0,bottom-summary_height+1,r.x1,bottom)
        index=color_index(value)
        page.draw_rect(box,color=(.6,.6,.6),fill=tuple(rgb[index]),width=.2665)
        boxes.append(dict(method=method,mean=float(value),palette_index=index,
                          color=palette.iloc[index]['hex'],bbox=list(box)))
    # Compare serialized documents, including all unchanged panels and labels.
    after=fitz.open(stream=doc.tobytes(),filetype='pdf');zoom=4
    assert all(s['text'] in after[0].get_text() for s in labels), 'Dataset labels changed'
    pix=[p[0].get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False) for p in [before,after]]
    a,b=[np.frombuffer(p.samples,dtype=np.uint8).reshape(p.height,p.width,3) for p in pix]
    mask=np.ones(a.shape[:2],dtype=bool)
    mask[int(region.y0*zoom)-2:int(region.y1*zoom)+3,int(region.x0*zoom)-2:int(region.x1*zoom)+3]=False
    delta=np.abs(a.astype(np.int16)-b.astype(np.int16))[mask]
    assert delta.max()==0,'Mean row changed artwork outside the heatmap body'
    centers=[(s['rect'].x0+s['rect'].x1)/2 for s in strips]+[541.4231]
    for x in centers:
        for i in range(10):
            y=top+(i+.5)*(bottom-top)/10
            assert np.array_equal(a[int(y*zoom),int(x*zoom)],b[int(ymap(y)*zoom),int(x*zoom)])
    row=[dict(text=w[4],bbox=w[:4]) for w in after[0].get_text('words')
         if w[0]>378 and w[1]>bottom-summary_height+1 and w[3]<bottom+1]
    assert not row,'Mean row must not contain numerical labels'
    for expected in boxes:
        r=fitz.Rect(expected['bbox']);x,y=(r.x0+r.x1)/2,(r.y0+r.y1)/2
        actual=b[int(y*zoom),int(x*zoom)].astype(float)
        assert np.max(np.abs(actual-rgb[expected['palette_index']]*255))<=1
    return dict(status='pass',label='Mean',numeric_annotations=False,
        method_count=len(means),colormap='RdYlBu_r',color_limits=[-1,1],
        palette_sha256=sha(palette_path),source_cell_colors_checked=len(source_cells),
        mean_cell_colors_checked=len(boxes),alignment='Aligned to the method heatmap columns',
        cells=boxes,best_performance_summary=None,changed_region=list(region),
        unchanged_area_max_pixel_difference=int(delta.max()),dataset_cell_centers_checked=200,
        dataset_row_height_factor=factor)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,default=Path('demo_output/figure3h_mean'))
    parser.add_argument('--full-figure',action='store_true');args=parser.parse_args()
    doc,audit=reorder();name='Figure03' if args.full_figure else 'Figure03h'
    save_artwork(doc,args.out/name,clip=None if args.full_figure else CLIP)
    (args.out/'validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    pd.DataFrame({'Method':audit['order'],'Mean':[audit['means'][m] for m in audit['order']],
                  'Valid datasets':[audit['valid_datasets'][m] for m in audit['order']]}).to_csv(args.out/'method_order.csv',index=False)
    print('Figure03h: 19 methods sorted by descending mean with a Mean row; vector/color and unchanged-area checks passed')

if __name__=='__main__':main()
