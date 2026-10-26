"""Reweight the retained 19-configuration composite; preserve all endpoint inputs."""
from pathlib import Path
import argparse, json, sys
import numpy as np
import pandas as pd
import pymupdf as fitz
import matplotlib
matplotlib.use('Agg')
from matplotlib import colormaps

HERE = Path(__file__).resolve().parent
WEIGHTS = {'Accuracy':.5, 'Scalability':.25, 'Stability':.25}
TITLE = 'Weighted composite (Accuracy 0.50, Scalability 0.25, Stability 0.25)'

def paths():
    for root in HERE.parents:
        if (root/'12_Provenance/current_sources.json').exists():
            return root/'06_Figure06/data/source_data/Figure06.xlsx', root/'09_Core/visualization/mean_top3'
        if (root/'source_data/Figure06.xlsx').exists():
            return root/'source_data/Figure06.xlsx', root/'Code/mean_top3'
    raise RuntimeError('Cannot locate workbook and vector helper')

def recalculate(summary):
    result = summary.copy()
    assert len(result)==19 and result.index.is_unique
    ranks = result[list(WEIGHTS)].rank(method='average')
    assert ranks.notna().all().all()
    for name in WEIGHTS:
        assert np.allclose(result[name+'_rank_index'],ranks[name]/len(result),atol=1e-12,rtol=0)
    # Integer arithmetic preserves exact mathematical ties.
    result['Overall_Score'] = (2*ranks.Accuracy+ranks.Scalability+ranks.Stability)/(4*len(result))
    result['Overall_Score_rank'] = result.Overall_Score.rank(ascending=False,method='min')
    return result

def pixels(page):
    pix=page.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False)
    return np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,3)

def render(source=None):
    workbook,helper=paths();sys.path.insert(0,str(helper))
    from mean_top3 import visible_drawings
    saved=pd.read_excel(workbook,sheet_name='6c Summary').set_index('Method')
    old=pd.read_csv(HERE/'assets/summary_before.csv',index_col='Method',float_precision='round_trip')
    current=recalculate(saved)
    # Only the derived total and its rank may differ from the retained inputs.
    columns=[c for c in old if c not in ['Overall_Score','Overall_Score_rank']]
    assert saved.index.equals(old.index) and np.allclose(saved[columns],old[columns],atol=1e-12,rtol=0,equal_nan=True)
    source=source or fitz.open(HERE/'assets/Figure06_before.pdf')
    doc=fitz.open();doc.insert_pdf(source,from_page=0,to_page=0);page=doc[0]
    bars=sorted([d for d in visible_drawings(page) if 72<d['rect'].x0<74 and 490<d['rect'].y0<700 and d.get('fill')],key=lambda d:d['rect'].y0)
    assert len(bars)==19
    spans=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
    labels=sorted([s for s in spans if s['text'] in old.index and 24<s['bbox'][0]<26 and 489<s['bbox'][1]<700],key=lambda s:s['bbox'][1])
    order=[s['text'] for s in labels];assert len(order)==19 and set(order)==set(old.index)
    widths=np.array([d['rect'].width for d in bars])/old.reindex(order).Overall_Score.to_numpy()
    assert np.ptp(widths)<.004  # Original serialized PDF coordinates are rounded.
    fullwidth=float(np.mean(widths));left=bars[0]['rect'].x0
    title=next(s for s in spans if s['text']=='Historical composite (19 configurations)')
    font=fitz.Font(fontname='helv');size=title['size']
    titlebox=fitz.Rect(title['bbox']);titlebox.x1=max(titlebox.x1,title['origin'][0]+font.text_length(TITLE,fontsize=size))
    titlebox += (-1,-1,1,1)
    column=fitz.Rect(left-.4,bars[0]['rect'].y0-.4,left+fullwidth+.4,bars[-1]['rect'].y1+.4)
    for box in [column,titlebox]:page.add_redact_annot(box,fill=(1,1,1))
    page.apply_redactions(images=0,graphics=0)
    page.insert_font(fontname='CompositeLabel',fontbuffer=font.buffer)
    page.insert_text(title['origin'],TITLE,fontname='CompositeLabel',fontsize=size,color=(.4,.26,.15))
    plotted=current.reindex(order)
    colors=colormaps['Reds_r']((plotted.Overall_Score.rank(method='average')-1)/(len(current)-1))
    for i,((method,row),oldbar,color) in enumerate(zip(plotted.iterrows(),bars,colors)):
        rect=fitz.Rect(oldbar['rect']);center=(rect.y0+rect.y1)/2
        # The original alternating backgrounds extend through this column.
        background=(214/255,)*3 if i%2==0 else (1,1,1)
        page.draw_rect(fitz.Rect(column.x0,rect.y0,column.x1,rect.y1),fill=background,color=None)
        rect.x1=left+fullwidth*row.Overall_Score
        page.draw_rect(rect,color=oldbar['color'],fill=tuple(color[:3]),width=oldbar['width'])
        if row.Overall_Score_rank<=6:
            label=str(int(row.Overall_Score_rank));fs=4.813009262084961
            x=(rect.x0+rect.x1-font.text_length(label,fontsize=fs))/2
            page.insert_text((x,center+fs*.187),label,fontname='CompositeLabel',fontsize=fs)
    doc=fitz.open(stream=doc.tobytes(garbage=4,deflate=True),filetype='pdf')
    a,b=pixels(source[0]),pixels(doc[0]);mask=np.ones(a.shape[:2],bool)
    for box in [column,titlebox]:
        x0,y0,x1,y1=box;mask[int(y0*2)-3:int(y1*2)+4,int(x0*2)-3:int(x1*2)+4]=False
    assert np.array_equal(a[mask],b[mask]), 'Artwork changed outside overall-score column and title'
    ranking=current.sort_values('Overall_Score',ascending=False,kind='stable')
    audit=dict(status='pass',weights=WEIGHTS,methods=19,top3=ranking.index[:3].tolist(),
               rank_normalization='Ascending average rank / 19',tie_rule='Descending minimum rank',
               unchanged_endpoint_inputs=True,unchanged_other_artwork_pixels=True,
               changed_regions=[list(column),list(titlebox)],models_rerun=False,
               caption_update_pending=True)
    return doc,current,audit

def export(doc,stem,crop=None):
    stem=Path(stem);stem.parent.mkdir(parents=True,exist_ok=True)
    if crop:
        other=fitz.open();page=other.new_page(width=crop.width,height=crop.height);page.show_pdf_page(page.rect,doc,0,clip=crop);doc=other
    doc.save(stem.with_suffix('.pdf'),garbage=4,deflate=True)
    doc[0].get_pixmap(matrix=fitz.Matrix(3,3)).save(stem.with_suffix('.png'))
    stem.with_suffix('.svg').write_text(doc[0].get_svg_image(text_as_path=True))

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();doc,current,audit=render()
    export(doc,args.out/'Figure06c_weighted',fitz.Rect(18,414,313,731))
    current.to_csv(args.out/'Figure06c_weighted_summary.csv',float_format='%.17g')
    (args.out/'validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))

if __name__=='__main__':main()
