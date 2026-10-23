"""Shared saved-data readers and plotting primitives; no panel execution on import."""
from pathlib import Path
import argparse, json, sys, os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/rna-velocity-panel-mpl")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from release_paths import REPO_ROOT as ROOT, source_path

def read(file,sheet):return pd.read_excel(source_path('source_data/' + file),sheet_name=sheet)

def matrix(file,sheet):
        d=read(file,sheet);return d.set_index(d.columns[0]).apply(pd.to_numeric,errors='coerce')

def save(fig,name,source,out,log):
        fig.savefig(out/(name+'.pdf'),bbox_inches='tight')
        fig.savefig(out/(name+'.svg'),bbox_inches='tight')
        plt.close(fig);log.append(dict(output=name,source=source,scope='saved-data redraw; not original layout'))

def heat(data,name,title,source,out,log):
        fig,ax=plt.subplots(figsize=(min(15,max(6,data.shape[1]*.17)),max(4,data.shape[0]*.19)),constrained_layout=True)
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#eeeeee')
        im=ax.imshow(data.astype(float),aspect='auto',cmap=cmap,vmin=-1,vmax=1,interpolation='nearest')
        ax.set_yticks(range(len(data)));ax.set_yticklabels(data.index)
        if data.shape[1]<=45:
            ax.set_xticks(range(len(data.columns)));ax.set_xticklabels(data.columns,rotation=90,fontsize=5)
        else:ax.set_xlabel('Dataset order in the source workbook')
        ax.set_title(title);fig.colorbar(im,ax=ax,label='Signed score; grey = unavailable',shrink=.7)
        save(fig,name,source,out,log)

def bar(values,name,title,source,xlabel,out,log):
        values=values.sort_values()
        fig,ax=plt.subplots(figsize=(7,max(4,len(values)*.22)),constrained_layout=True)
        ax.barh(values.index,values,color='#40769c');ax.set_title(title);ax.set_xlabel(xlabel)
        save(fig,name,source,out,log)

def run_panel(render):
    parser=argparse.ArgumentParser(description=render.__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'panel_checks')
    args=parser.parse_args();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    log=[]
    plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none','font.size':7})
    render(out,log)
    registry=json.loads(Path(__file__).with_name('panel_registry.json').read_text())
    entry=next(e for e in registry if (ROOT/e['path']).resolve()==Path(render.__code__.co_filename).resolve())
    for row in log:
        row.update(figure=entry['figure'],panel=entry['panel'],entrypoint=entry['path'],scope=entry['scope'])
    for row in log:
        (out/(row['output']+'.json')).write_text(json.dumps(dict(**row,models_rerun=False,publication_artwork_modified=False),indent=2)+'\n')
    print(json.dumps(log,indent=2))
