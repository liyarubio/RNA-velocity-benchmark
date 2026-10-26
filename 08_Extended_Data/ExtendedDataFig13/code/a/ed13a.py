"""Extended Data 13: saved GDC over read depth, with Figure 6b method colors."""
from pathlib import Path
import argparse, hashlib, json, sys
import numpy as np
import pandas as pd
import pymupdf as fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
ALIASES={
 'scVelo (stochastic)':'scVelo(sto)','scvelo_stochastic':'scVelo(sto)',
 'scVelo (dynamic)':'scVelo(dyn)','scvelo_dynamical':'scVelo(dyn)',
 'veloAE':'VeloAE','Dynamo':'dynamo','veloVAE':'VeloVAE','uniTvelo':'UniTVelo',
 'DeepVelo (VAE-based)':'DeepVelo(VAE)','DeepVelo_SA':'DeepVelo(VAE)',
 'DeepVelo (GCN-based)':'DeepVelo(GCN)','DeepVelo_GB':'DeepVelo(GCN)',
 'latentvelo':'LatentVelo','stt':'STT','multivelo':'MultiVelo','phylovelo':'PhyloVelo',
 'topicvelo':'TopicVelo','bivi':'biVI','deepkinet':'DeepKINET','sdevelo':'SDEvelo'}

def locations():
    for root in HERE.parents:
        if (root/'12_Provenance/current_sources.json').exists():
            return root/'08_Extended_Data/ExtendedDataFig13/data/read_depth_gdc_20261026',root/'06_Figure06/figures/Figure06.pdf',root/'09_Core/visualization/mean_top3'
        if (root/'source_data/Figure06.xlsx').exists():
            return root/'data/read_depth_gdc_20261026',root/'figures/main/Figure06.pdf',root/'Code/mean_top3'
    raise RuntimeError('Cannot locate release inputs')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def palette_from_main():
    _,pdf,helper=locations();sys.path.insert(0,str(helper))
    from mean_top3 import visible_drawings
    page=fitz.open(pdf)[0]
    labels=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']
            if 325<s['bbox'][1]<400 and s['bbox'][0]>310]
    assert len(labels)==21 and len({s['text'] for s in labels})==21
    labels.sort(key=lambda s:(round(s['bbox'][0]),s['bbox'][1]))
    lines=[d for d in visible_drawings(page) if d.get('color') and d['rect'].height<.01 and 5<d['rect'].width<10]
    palette=[]
    for s in labels:
        x0,y0,x1,y1=s['bbox'];candidates=[d for d in lines if x0-15<d['rect'].x0<x0 and y0<d['rect'].y0<y1]
        assert len(candidates)==1,(s['text'],len(candidates))
        d=candidates[0];rgb=list(d['color'])
        palette.append(dict(Method=s['text'],r=rgb[0],g=rgb[1],b=rgb[2],
                            hex='#'+''.join(f'{round(x*255):02x}' for x in rgb),
                            main_legend_bounds=list(d['rect'])))
    return palette,sha(pdf)

def data_from_saved():
    data,_,_=locations();provenance=json.loads((data/'provenance.json').read_text());frames=[]
    for tag in ['0_2','0_4','0_6','0_8','1_0']:
        path=data/f'inputs/Eval_{tag}.csv';assert sha(path)==provenance[path.name]['sha256']
        df=pd.read_csv(path,float_precision='round_trip');df['Source_method']=df.Method
        df.Method=df.Method.replace(ALIASES);df['Read_fraction']=float(tag.replace('_','.'));df['Source_file']=path.name
        assert len(df)==21 and not df.Method.duplicated().any()
        frames.append(df)
    frame=pd.concat(frames,ignore_index=True)
    assert len(frame)==105 and np.isfinite(frame.GDC).all() and frame.GDC.between(-1,1).all()
    frame['Interpretation']='Saved GDC score'
    frame.loc[frame.Method.eq('TopicVelo'),'Interpretation']='Archived base-scVelo proxy; not native integrated transitions'
    return frame

def render(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    data=data_from_saved();palette,main_hash=palette_from_main();order=[x['Method'] for x in palette]
    assert set(data.Method)==set(order)
    pivot=data.pivot(index='Read_fraction',columns='Method',values='GDC')
    assert np.array_equal(pivot.TopicVelo,pivot['scVelo(sto)'])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(210/25.4,297/25.4));ax=fig.add_axes([.14,.40,.81,.46])
    fig.text(.5,.935,'Extended Data Fig. 13 | Read-depth sensitivity of GDC',ha='center',fontsize=12,weight='bold')
    fig.text(.5,.90,'ReDeeM young2 HSPC',ha='center',fontsize=11)
    plot_lines=[]
    for color in palette:
        method=color['Method'];proxy=method=='TopicVelo';rgb=[color[k] for k in ['r','g','b']]
        line,=ax.plot(pivot.index,pivot[method],color=rgb,marker='s' if proxy else 'o',
                      markersize=5 if proxy else 4,markerfacecolor='none' if proxy else rgb,
                      linestyle='--' if proxy else '-',lw=1.3,
                      label=method+(' (base-scVelo proxy)' if proxy else ''),zorder=5 if proxy else 3)
        assert np.array_equal(np.asarray(line.get_color()),rgb)
        assert np.array_equal(line.get_ydata(),pivot[method].to_numpy())
        plot_lines.append(line)
    ax.set(xlabel='Read-depth fraction',ylabel='Global Direction-weighted Coherence (GDC)',xlim=(.17,1.03),ylim=(-1,1))
    ax.set_xticks([.2,.4,.6,.8,1.]);ax.set_yticks(np.arange(-1,1.01,.25))
    ax.axhline(0,color='#777777',lw=.7,zorder=0);ax.grid(alpha=.22,lw=.6);ax.spines[['top','right']].set_visible(False)
    legend=ax.legend(loc='upper center',bbox_to_anchor=(.5,-.16),ncol=3,frameon=False,fontsize=8.5,handlelength=2.2,columnspacing=1.25)
    fig.text(.14,.105,'TopicVelo: archived base-scVelo proxy; native integrated transitions were not evaluated here.\n'
             'Full-depth GDC is an evaluated score, not a self-comparison fixed to 1.',fontsize=8,color='#555555',linespacing=1.6)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    for artist in [legend,*fig.texts,ax.xaxis.label,ax.yaxis.label]:
        bounds=artist.get_window_extent(renderer);assert fig.bbox.contains(bounds.x0,bounds.y0) and fig.bbox.contains(bounds.x1,bounds.y1)
    for ext in ['pdf','svg','png']:fig.savefig(out/f'Extended_Data_Fig13.{ext}',dpi=180)
    plt.close(fig)
    pdf=fitz.open(out/'Extended_Data_Fig13.pdf');page=pdf[0]
    assert np.allclose([page.rect.width,page.rect.height],[210*72/25.4,297*72/25.4],atol=.01)
    # Compare serialized 8-bit RGB; PDF DeviceGray conversion can shift floats.
    vectors=[d for d in page.get_drawings() if d.get('color') and d['rect'].width>300]
    for color in palette:
        rgb=np.array([color[k] for k in ['r','g','b']]);assert any(np.array_equal(np.rint(np.asarray(d['color'])*255),np.rint(rgb*255)) for d in vectors),color['Method']
    data.to_csv(out/'plotted_GDC.csv',index=False,float_format='%.17g')
    pd.DataFrame(palette).drop(columns='main_legend_bounds').to_csv(out/'method_palette.csv',index=False,float_format='%.17g')
    report=dict(status='pass',figure='ED13',points=105,methods=21,read_fractions=list(pivot.index),
                method_order=order,palette=palette,main_figure_sha256=main_hash,
                colors_match_main_fig6b=True,serialized_pdf_8bit_rgb_verified=True,
                numerical_inputs_unchanged=True,models_rerun=False,metrics_recomputed=False,
                page_size='A4 portrait',topicvelo='Explicitly labeled archived base-scVelo proxy',
                caption_author_managed=True)
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n');return report

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();r=render(args.out);print(json.dumps({k:r[k] for k in ['status','figure','points','methods','colors_match_main_fig6b']},indent=2))

if __name__=='__main__':main()
