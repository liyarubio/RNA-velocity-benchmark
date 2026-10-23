"""Reflow the published vector panels to Mean-only, with aligned top-three ranks.

This presentation renderer uses a frozen vector template for the unchanged
matrix/annotations. Means and ranks are read from the current numerical inputs;
bar widths are independently checked against their original joint Mean/Median
scale. It does not fit models. Median is retained only for scale verification.
"""
from pathlib import Path
import argparse, hashlib, json, re, io
import numpy as np
import pandas as pd
import pymupdf as fitz

HERE = Path(__file__).resolve().parent
ORDER = ['velocyto','scVelo(sto)','scVelo(dyn)','dynamo','VeloVAE','UniTVelo',
         'DeepVelo(VAE)','cellDancer','veloVI','LatentVelo','DeepVelo(GCN)',
         'STT','biVI','DeepKINET','TopicVelo','SDEvelo','cell2fate']
ALIASES = dict(zip(['celldancer','scvelostochastic','scvelodynamical','stt','unitvelo',
 'deepvelosa','deepvelogb','velovi','latentvelo','dynamostochastic','dynamodeterministic',
 'velocyto','velovae','bivi','cell2fate','deepkinet','sdevelo','topicvelo'],
 ['cellDancer','scVelo(sto)','scVelo(dyn)','STT','UniTVelo','DeepVelo(VAE)',
 'DeepVelo(GCN)','veloVI','LatentVelo','dynamo','Dynamo_deterministic','velocyto',
 'VeloVAE','biVI','cell2fate','DeepKINET','SDEvelo','TopicVelo']))
SPECS = {
 'Figure03a': {'base':'Figure03_artwork.pdf','band':(85,204),'clip':(32,54,560,280)},
 'ED04a': {'base':'Extended_Data_Fig04.pdf','band':(110,258),'clip':(39,39,551,290)},
 'ED04b': {'base':'Extended_Data_Fig04.pdf','band':(357,505),'clip':(39,290,551,534)},
}

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def repository_root():
    for p in HERE.parents:
        if (p/'12_Provenance/current_sources.json').is_file(): return p, True
        if (p/'source_data/Figure03.xlsx').is_file(): return p, False
    raise RuntimeError('Cannot locate the current numerical inputs')

def data_for(panel):
    root, public = repository_root()
    if panel == 'Figure03a':
        path = root / ('03_Figure03/data/figure3_detail' if public else 'data/figure3_detail') / 'Figure3a_high_dim_cosine_similarity.csv'
        table = pd.read_csv(path,index_col=0)
        table.index = [ALIASES[re.sub('[^a-z0-9]','',re.sub(r'^\d+_','',s).lower())] for s in table.index]
        table = table.reindex(ORDER)
    else:
        path = root / ('08_Extended_Data/ExtendedDataFig04/data/source_data' if public else 'source_data') / 'Extended_Data_Fig05.xlsx'
        table = pd.read_excel(path,sheet_name='Gene-level cosine' if panel=='ED04a' else 'Cell-level PCA2 cosine',index_col=0)
    expected=json.loads((HERE/'assets/input_hashes.json').read_text())[panel]
    assert sha(path) in ([expected] if isinstance(expected,str) else expected), 'Numerical input changed: rebuild the full matrix and template first'
    assert table.shape == (17,42), table.shape
    return path, table

def visible_drawings(page):
    """Respect nested PDF clipping; get_drawings() alone includes hidden copies."""
    stack = {}
    for item in page.get_drawings(extended=True):
        level = item['level']; stack = {k:v for k,v in stack.items() if k<level}
        if item['type']=='clip': stack[level]=item['scissor']; continue
        if item['type']=='group': continue
        clip = fitz.Rect(page.rect)
        for rect in stack.values(): clip &= rect
        if clip.contains(item['rect']): yield item

def text_spans(page):
    for block in page.get_text('dict')['blocks']:
        for line in block.get('lines',[]):
            for span in line['spans']: yield span, line['dir']

def original_font(doc,page,name):
    normalized=lambda s: re.sub('[^a-z]','',s.split('+')[-1].lower().replace('#20',''))
    for entry in page.get_fonts(full=True):
        if normalized(entry[3]) == normalized(name):
            data=doc.extract_font(entry[0])[3]
            if data:
                font = fitz.Font(fontbuffer=data)
                if all(font.has_glyph(ord(c)) for c in 'Mean123'): return font
                # Some published CID subsets omit the TrueType Unicode cmap.
                # Restore it from the PDF ToUnicode map, preserving glyph outlines.
                from fontTools.ttLib import TTFont, newTable
                from fontTools.ttLib.tables._c_m_a_p import CmapSubtable
                tt = TTFont(io.BytesIO(data)); glyphs = tt.getGlyphOrder()
                ref = int(doc.xref_get_key(entry[0], 'ToUnicode')[1].split()[0])
                mapping = doc.xref_stream(ref).decode('ascii'); cmap = {}
                for block in re.findall(r'beginbfrange(.*?)endbfrange', mapping, re.S):
                    for a,b,c in re.findall(r'<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>',block):
                        start,end,code = int(a,16),int(b,16),int(c,16)
                        for gid in range(start,min(end+1,len(glyphs))): cmap[code+gid-start] = glyphs[gid]
                for block in re.findall(r'beginbfchar(.*?)endbfchar', mapping, re.S):
                    for a,b in re.findall(r'<([0-9A-Fa-f]+)>\s*<([0-9A-Fa-f]+)>',block):
                        if int(a,16)<len(glyphs): cmap[int(b,16)] = glyphs[int(a,16)]
                sub = CmapSubtable.newSubtable(4);sub.platformID=3;sub.platEncID=1;sub.language=0
                sub.cmap={k:v for k,v in cmap.items() if k<65536}
                tt['cmap']=newTable('cmap');tt['cmap'].tableVersion=0;tt['cmap'].tables=[sub]
                stream=io.BytesIO();tt.save(stream);font=fitz.Font(fontbuffer=stream.getvalue())
                assert all(font.has_glyph(ord(c)) for c in 'Mean123')
                return font
    raise RuntimeError('Published font not found: '+name)

def render_figure(panels, template=None):
    spec=SPECS[panels[0]]; source=Path(template) if template else HERE/'assets'/spec['base']
    original=fitz.open(source); out=fitz.open();out.insert_pdf(original)
    page=out[0]; audits=[]; changes=[]
    for panel in panels:
        cfg=SPECS[panel]; low,high=cfg['band']; data_path,table=data_for(panel)
        mean=table.mean(axis=1); median=table.median(axis=1)
        lo=min(mean.min(),median.min())-.01; hi=max(mean.max(),median.max())
        normalized=(mean-lo)/(hi-lo); ranks=mean.rank(ascending=False,method='min')
        draws=list(visible_drawings(original[0]))
        headers={s['text']:(s,direction) for s,direction in text_spans(original[0])
                 if s['text'] in ['Mean','Median'] and s['bbox'][0]>490 and low<s['origin'][1]<high}
        assert set(headers)=={'Mean','Median'}
        bars=[d for d in draws if d.get('fill') and d.get('color') and d['rect'].x0>490
              and low<d['rect'].y0<high and 5<d['rect'].height<8
              and all(x[0] in ['l','re'] for x in d['items'])]
        # Rectangles form exactly two 17-row summary columns.
        assert len(bars)==34,(panel,len(bars))
        xs=sorted(set(round(d['rect'].x0,3) for d in bars));assert len(xs)==2
        medbars=sorted([d for d in bars if abs(d['rect'].x0-xs[0])<.01],key=lambda d:d['rect'].y0)
        meanbars=sorted([d for d in bars if abs(d['rect'].x0-xs[1])<.01],key=lambda d:d['rect'].y0)
        dx=meanbars[0]['rect'].x0-medbars[0]['rect'].x0
        scale=np.median([d['rect'].width/v for d,v in zip(meanbars,normalized)])
        errors=[abs(d['rect'].width-scale*v) for d,v in zip(meanbars,normalized)]
        assert max(errors)<.008,(panel,max(errors))
        stripes=[d for d in draws if d.get('fill') and d['rect'].x0<100
                 and d['rect'].x1>meanbars[-1]['rect'].x1 and low<d['rect'].y0<high
                 and 5<d['rect'].height<8]
        assert len(stripes)==9,(panel,len(stripes))
        left=medbars[0]['rect'].x0-.25; right=max(d['rect'].x1 for d in stripes)+.4
        region=fitz.Rect(left,min(headers[k][0]['bbox'][1] for k in headers)-.1,
                         right,max(d['rect'].y1 for d in stripes)+.25)
        # Physically remove summary text/contained paths, then rebuild only this region.
        page.add_redact_annot(region,fill=(1,1,1),cross_out=False)
        page.apply_redactions(images=0,graphics=1,text=0)
        for d in stripes:
            r=fitz.Rect(d['rect']);r.x0=left;r.x1-=dx
            page.draw_rect(r,color=None,fill=d['fill'],overlay=True)
        for d in meanbars:
            r=fitz.Rect(d['rect']);r.x0-=dx;r.x1-=dx
            page.draw_rect(r,color=d['color'],fill=d['fill'],width=d['width'],overlay=True)
        # Use the embedded original font; no system-font substitution is needed.
        span,direction=headers['Mean']; font=original_font(original,original[0],span['font'])
        fontname='MeanTop3'+panel; page.insert_font(fontname=fontname,fontbuffer=font.buffer)
        origin=fitz.Point(span['origin'][0]-dx,span['origin'][1])
        angle=np.degrees(np.arctan2(-direction[1],direction[0]))
        page.insert_text(origin,'Mean',fontname=fontname,fontsize=span['size'],
                         morph=(origin,fitz.Matrix(angle)))
        old_labels=[s for s,_ in text_spans(original[0]) if s['text'] in ['1','2','3','4','5','6']
                    and s['bbox'][0]>xs[1] and low<s['origin'][1]<high]
        assert len(old_labels)>=6
        size=old_labels[0]['size']; rank_font=original_font(original,original[0],old_labels[0]['font'])
        rankname='RankTop3'+panel;page.insert_font(fontname=rankname,fontbuffer=rank_font.buffer)
        first=old_labels[0]; nearest=min(meanbars,key=lambda d:abs((d['rect'].y0+d['rect'].y1)/2-first['origin'][1]))
        baseline_offset=first['origin'][1]-(nearest['rect'].y0+nearest['rect'].y1)/2
        rank_x=medbars[0]['rect'].x0+.4*scale
        labels=[]
        for (method,value),rank,d in zip(mean.items(),ranks,meanbars):
            if rank>3: continue
            text=str(int(rank));r=d['rect'];y=(r.y0+r.y1)/2
            half=rank_font.text_length(text,fontsize=size)/2
            assert rank_x-half>r.x0-dx and rank_x+half<r.x1-dx
            page.insert_text((rank_x-half,y+baseline_offset),text,fontname=rankname,fontsize=size)
            labels.append(dict(method=method,rank=int(rank),mean=float(value),x=rank_x,y=y))
        visible=page.get_text('words',clip=region)
        assert 'Median' not in [w[4] for w in visible]
        assert sorted(w[4] for w in visible if w[4].isdigit())==sorted(str(v['rank']) for v in labels),(panel,visible,labels)
        audits.append(dict(panel=panel,input_sha256=sha(data_path),baseline_sha256=sha(source),
            mean=mean.to_dict(),rank=ranks.astype(int).to_dict(),method_order=list(table.index),
            valid_datasets=table.count(axis=1).to_dict(),normalization_min=float(lo),normalization_max=float(hi),
            normalized_mean=normalized.to_dict(),bar_width_check_max_error_pt=max(errors),
            labels=labels,changed_rect=list(region),median_displayed=False,rank_limit=3,
            ties='minimum rank, retained',template='Frozen published vectors; original Mean fills and scale retained'))
        changes.append(region)
    # Compare all pixels outside the edited summary rectangles, including every other panel.
    zoom=3; before=original[0].get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False)
    after=page.get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False)
    a=np.frombuffer(before.samples,dtype=np.uint8).reshape(before.height,before.width,3)
    b=np.frombuffer(after.samples,dtype=np.uint8).reshape(after.height,after.width,3)
    mask=np.ones(a.shape[:2],dtype=bool)
    for r in changes:
        mask[max(0,int(r.y0*zoom)-2):int(r.y1*zoom)+3,max(0,int(r.x0*zoom)-2):int(r.x1*zoom)+3]=False
    delta=np.abs(a.astype(np.int16)-b.astype(np.int16))[mask]
    result=dict(status='pass',panels=audits,unchanged_area_mean_absolute_error=float(delta.mean()),
                unchanged_area_changed_fraction=float(np.mean(delta>0)),models_rerun=False,numeric_inputs_changed=False)
    assert result['unchanged_area_mean_absolute_error']<.01,result
    assert result['unchanged_area_changed_fraction']<.001,result
    return out,result

def save_artwork(doc,stem,clip=None):
    stem=Path(stem);stem.parent.mkdir(parents=True,exist_ok=True)
    if clip:
        cropped=fitz.open();page=cropped.new_page(width=clip.width,height=clip.height)
        page.show_pdf_page(page.rect,doc,0,clip=clip);doc=cropped
    doc.save(stem.with_suffix('.pdf'),garbage=4,deflate=True)
    preview=doc; zoom=3.5
    if stem.name == 'Figure03' and clip is None:
        # Retain the existing 600 x 780 pt web-preview frame and 300 dpi PNG.
        marker=next(s for s,_ in text_spans(doc[0]) if s['text']=='a')
        scale=marker['size']/9
        x=marker['origin'][0]-5*scale; y=marker['origin'][1]-11*scale
        frame=fitz.Rect(x,y,x+600*scale,y+780*scale)
        preview=fitz.open();p=preview.new_page(width=600,height=780)
        p.show_pdf_page(p.rect,doc,0,clip=frame); zoom=300/72
    preview[0].get_pixmap(matrix=fitz.Matrix(zoom,zoom)).save(stem.with_suffix('.png'))
    stem.with_suffix('.svg').write_text(preview[0].get_svg_image(text_as_path=False))

def main(panel):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('demo_output/mean_top3'))
    args=parser.parse_args();doc,audit=render_figure([panel])
    save_artwork(doc,args.out/panel,fitz.Rect(SPECS[panel]['clip']))
    (args.out/(panel+'_validation.json')).write_text(json.dumps(audit,indent=2)+'\n')
    print(panel+': Mean only; top-three ranks aligned; validation passed')
