"""Portable numerical-panel redraws for checking the retained figures.

These plots use the saved figure-level data. They intentionally do not replace
the published layouts, manual schematics, historical fields or original images.
Usage: python 09_Core/visualization/plot_saved_panels.py --out panel_checks
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from release_paths import REPO_ROOT as ROOT, source_path

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'panel_checks')
    args=parser.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none','font.size':7})
    log=[]
    def read(file,sheet):return pd.read_excel(source_path('source_data/' + file),sheet_name=sheet)
    def matrix(file,sheet):
        d=read(file,sheet);return d.set_index(d.columns[0]).apply(pd.to_numeric,errors='coerce')
    def save(fig,name,source):
        fig.savefig(out/(name+'.pdf'),bbox_inches='tight')
        fig.savefig(out/(name+'.svg'),bbox_inches='tight')
        plt.close(fig);log.append(dict(output=name,source=source,scope='saved-data redraw; not original layout'))
    def heat(data,name,title,source):
        fig,ax=plt.subplots(figsize=(min(15,max(6,data.shape[1]*.17)),max(4,data.shape[0]*.19)),constrained_layout=True)
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#eeeeee')
        im=ax.imshow(data.astype(float),aspect='auto',cmap=cmap,vmin=-1,vmax=1,interpolation='nearest')
        ax.set_yticks(range(len(data)));ax.set_yticklabels(data.index)
        if data.shape[1]<=45:
            ax.set_xticks(range(len(data.columns)));ax.set_xticklabels(data.columns,rotation=90,fontsize=5)
        else:ax.set_xlabel('Dataset order in the source workbook')
        ax.set_title(title);fig.colorbar(im,ax=ax,label='Signed score; grey = unavailable',shrink=.7)
        save(fig,name,source)
    def bar(values,name,title,source,xlabel):
        values=values.sort_values()
        fig,ax=plt.subplots(figsize=(7,max(4,len(values)*.22)),constrained_layout=True)
        ax.barh(values.index,values,color='#40769c');ax.set_title(title);ax.set_xlabel(xlabel)
        save(fig,name,source)

    heat(matrix('Figure03.xlsx','3a High dim'),'Figure03a','Simulation: high-dimensional cosine','Figure03.xlsx / 3a High dim')
    f4=read('Figure04.xlsx','4a GDC').set_index('Method')
    bar(f4.GDC,'Figure04a_GDC','ReDeeM: saved GDC','Figure04.xlsx / 4a GDC','GDC')
    f5=read('Figure05.xlsx','5a GDC').set_index('Method')
    bar(f5.GDC,'Figure05a_GDC','ST_small: 608 spots; current native TopicVelo','Figure05.xlsx / 5a GDC','GDC')
    rt=read('Figure05.xlsx','5b Runtime seconds').set_index('Cell')
    fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
    for c in rt:ax.plot(rt.index,rt[c],label=c,lw=.8)
    ax.set(xlabel='Cells',ylabel='Seconds',yscale='log',title='Figure 5b: recorded runtime');ax.legend(fontsize=5,ncol=3)
    save(fig,'Figure05b_runtime','Figure05.xlsx / 5b Runtime seconds')
    reads=read('Figure05.xlsx','5c Plotted vs full')
    fig,ax=plt.subplots(figsize=(9,5),constrained_layout=True)
    for m,g in reads.groupby('Method'):
        g=g.sort_values('Ratio_i');ax.plot(g.Ratio_i,g.CosineSim,marker='.',lw=.8,label=m)
    ax.set(xlabel='Read-depth fraction',ylabel='Cosine versus full depth',title='Figure 5c: consistency; TopicVelo uses a transition proxy');ax.legend(fontsize=5,ncol=3)
    save(fig,'Figure05c_stability','Figure05.xlsx / 5c Plotted vs full')
    f6=read('Figure06.xlsx','6a Summary').set_index('Method')
    bar(f6.Overall_Score,'Figure06a_frozen','Frozen composite; stability excluded','Figure06.xlsx / 6a Summary','0.6 accuracy rank + 0.4 scalability rank')
    for sheet,label in [('Gene-level cosine','a'),('Cell-level PCA2 cosine','b')]:
        heat(matrix('Extended_Data_Fig05.xlsx',sheet),'ED05'+label,'dyngen: '+sheet,'Extended_Data_Fig05.xlsx / '+sheet)
    pairs=read('Extended_Data_Fig06.xlsx','Paired cell-level values')
    fig,ax=plt.subplots(figsize=(6,5),constrained_layout=True)
    ax.scatter(pairs.high,pairs.low,s=4,alpha=.35);ax.plot([-1,1],[-1,1],c='grey',lw=.7)
    ax.set(xlabel='High-dimensional cosine',ylabel='PCA2 cosine',title='ED6: 698 paired method–dataset observations')
    save(fig,'ED06_paired','Extended_Data_Fig06.xlsx / Paired cell-level values')
    heat(matrix('Extended_Data_Fig07.xlsx','HuBMAP raw GDC').T,'ED07_HuBMAP','HuBMAP: 176 datasets; raw source order','Extended_Data_Fig07.xlsx / HuBMAP raw GDC')
    for n,name in [(13,'scNT'),(14,'scLT')]:
        frame=read(f'Extended_Data_Fig{n}.xlsx','Method summary').set_index('method')
        bar(frame['median'],f'ED{n}_median',name+': median per-cell reference cosine',f'Extended_Data_Fig{n}.xlsx / Method summary','Median cosine; saved pipelines')
    for n,sheet in [(15,'GDC source values'),(17,'GDC direct edges')]:
        frame=read(f'Extended_Data_Fig{n}.xlsx',sheet)
        name='method' if 'method' in frame else 'Method'
        bar(frame.set_index(name).GDC,f'ED{n}_GDC','ST_raw (1,292 spots)' if n==15 else 'KPTracer (754 cells; denominator 706)',f'Extended_Data_Fig{n}.xlsx / '+sheet,'GDC')
    paired=pd.read_csv(source_path('data/spatial_context_consistency/paired_spots.csv'))
    bar(paired.groupby('method').spatial_cosine.median(),'ED16_spatial','Same 608 target spots in two saved workflows','data/spatial_context_consistency/paired_spots.csv','Median spatial direction cosine')
    for n in [100,200]:
        members=pd.read_csv(source_path(f'data/metascape/top{n}/mapped_membership.csv'))
        counts=pd.read_csv(source_path(f'data/metascape/top{n}/list_counts.csv'))
        binary=members[counts.method].astype(int);overlap=binary.T @ binary
        fig,ax=plt.subplots(figsize=(6,5),constrained_layout=True)
        im=ax.imshow(overlap,cmap='Blues');ax.set_xticks(range(9));ax.set_xticklabels(overlap.columns,rotation=65,ha='right');ax.set_yticks(range(9));ax.set_yticklabels(overlap.index)
        for i in range(9):
            for j in range(9):ax.text(j,i,str(overlap.iloc[i,j]),ha='center',va='center',fontsize=6,color='white' if overlap.iloc[i,j]>n*.6 else 'black')
        ax.set_title(f'ED18: top{n} mapped-gene intersections');fig.colorbar(im,ax=ax,label='Shared mapped genes')
        save(fig,f'ED18_top{n}_overlap',f'data/metascape/top{n}/mapped_membership.csv')
    (out/'plot_manifest.json').write_text(json.dumps(dict(status='complete',plots=log,models_rerun=False,publication_artwork_modified=False),indent=2)+'\n')
    print(f'{len(log)} numerical-panel checks written to {out}')

if __name__=='__main__':main()
