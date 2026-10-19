"""Fill ED16 TFvelo/TopicVelo slots and score bars from completed reruns."""
from pathlib import Path
import os,sys,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/kptracer-completion-mpl')
B=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(B/'Code'),str(B/'Code/vendor')]
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pymupdf as fitz
import openpyxl
from ed16_completion import ready,results,D
for name in ['Arial.ttf','Arial_Bold.ttf']:font_manager.fontManager.addfont('/usr/share/fonts/truetype/msttcorefonts/'+name)
plt.rcParams.update({'font.family':'Arial','font.size':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.5})
def save(doc,target):
    doc.subset_fonts();tmp=target.with_suffix('.tmp.pdf');doc.save(tmp,garbage=4,deflate=True)
    doc[0].get_pixmap(matrix=fitz.Matrix(2.5,2.5)).save(target.with_suffix('.png'))
    target.with_suffix('.svg').write_text(doc[0].get_svg_image(text_as_path=False));doc.close();tmp.replace(target)
def render():
    assert ready()
    w=openpyxl.load_workbook(B/'Code/inputs/extended_v2/Extended_Data_Fig18.xlsx',read_only=True,data_only=True)
    palette={r[0]:r[2] for r in list(w['Cell-type palette'].values)[1:]}
    scores=dict(list(w['GDC direct edges'].values)[1:]);rr=results()
    for name,r in rr.items():scores[name]=r['GDC']
    # Regenerate UMAP/schematic first, so repeat calls never accumulate overlays.
    from plots.redraw_ed16_reference import render as reference
    reference(include_completion=False)
    target=B/'Assets/arial/Extended_Data_Fig16_overview.pdf';old=fitz.open(target)
    d=fitz.open();p=d.new_page(width=595.276,height=841.89)
    p.show_pdf_page(fitz.Rect(0,0,595.276,372),old,0,clip=fitz.Rect(0,0,595.276,372))
    fig,ax=plt.subplots(figsize=(7.0,3.2));fig.subplots_adjust(left=.09,right=.995,bottom=.35,top=.89)
    colors=['#D62728','#CAB2D6','#6A3D9A','#33A02C','#FFB547','#B2DF8A','#E6CF5C','#0000CD','#FFB6C1','#DB7093','#20B2AA','#87CEFA','#7B4173','#9ACD32','#FF8C00','#1E90FF','#FF6347','#4682B4','#CD5C5C']
    x=np.arange(len(scores));bars=ax.bar(x,list(scores.values()),color=colors,width=.7)
    for i in [12,13]:bars[i].set_edgecolor('#333');bars[i].set_linewidth(.6)
    ax.axhline(0,color='#555',lw=.5);ax.set_ylim(-1,1);ax.set_ylabel('GDC')
    ax.set_xticks(x);ax.set_xticklabels([n for n in scores],rotation=50,ha='right',fontsize=6.4)
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    ax.set_title('Four direct-evidence relations; 706 reference cells',loc='left',fontsize=8)
    fig.text(.09,.025,'Native-transition displacement proxy. Outlined bars: rerun methods.',fontsize=6.5)
    tmp=B/'Reports/ed16_completed_b.pdf';fig.savefig(tmp);plt.close(fig)
    with fitz.open(tmp) as plot:p.show_pdf_page(fitz.Rect(65,378,562,618),plot,0)
    p.insert_text((58,387),'b',fontname='RerunLabel',fontfile='/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf',fontsize=8)
    top=fitz.Rect(0,0,595.276,372);assert p.get_pixmap(clip=top).samples==old[0].get_pixmap(clip=top).samples
    assert not p.get_images();old.close();save(d,target)
    # Reuse the verified vector page for all 17 historical fields.
    from style_extended_typography import style_one,SOURCE
    stem='Extended_Data_Fig16_method_fields';labels=json.loads((B/'Code/inputs/typography_labels.json').read_text())
    style_one(SOURCE/(stem+'.pdf'),labels.get(stem+'.pdf',[]))
    target=B/'Assets/arial'/(stem+'.pdf');d=fitz.open(target);p=d[0]
    rectangles={'tfvelo':fitz.Rect(51,379,173,498),'topicvelo':fitz.Rect(177,379,299,498)}
    from scipy.spatial import cKDTree
    for key,rect in rectangles.items():
        data=pd.read_csv(D/(key+'_directions.csv'));xy=data[['x','y']].to_numpy();v=data[['vx','vy']].to_numpy()
        fig=plt.figure(figsize=(rect.width/72,rect.height/72));ax=fig.add_axes([.025,.025,.95,.84])
        for name,color in palette.items():
            mask=data.cell_type==name;ax.scatter(xy[mask,0],xy[mask,1],s=2,c=color,edgecolors='none',alpha=.9,rasterized=False)
        span=np.ptp(xy,axis=0);lo=xy.min(0)-span*.05;hi=xy.max(0)+span*.05
        grid=[np.linspace(lo[i],hi[i],35) for i in range(2)];xx,yy=np.meshgrid(*grid);points=np.c_[xx.ravel(),yy.ravel()]
        valid=np.isfinite(v).all(1)&(np.linalg.norm(v,axis=1)>0)
        dist,idx=cKDTree(xy[valid]).query(points,k=30);bw=np.mean(span)/35*.75
        weights=np.exp(-.5*(dist/bw)**2);mass=weights.sum(1)
        field=(v[valid][idx]*weights[:,:,None]).sum(1)/np.maximum(mass[:,None],1e-10);field[mass<.2]=np.nan
        ax.streamplot(*grid,np.ma.masked_invalid(field[:,0].reshape(xx.shape)),np.ma.masked_invalid(field[:,1].reshape(xx.shape)),density=1,color='#222',linewidth=.4,arrowsize=.48)
        ax.set(xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]));ax.set_aspect('auto');ax.axis('off')
        fig.text(.5,.955,'TFvelo' if key=='tfvelo' else 'TopicVelo',ha='center',va='top',fontsize=7,fontweight='bold')
        file=B/'Reports'/('ed16_'+key+'_field.pdf');fig.savefig(file);plt.close(fig)
        p.add_redact_annot(rect,fill=(1,1,1));p.apply_redactions(images=2,graphics=1,text=0)
        with fitz.open(file) as fieldpdf:p.show_pdf_page(rect,fieldpdf,0)
    assert not p.get_images();save(d,target)
    report={'status':'complete','umap_and_schematic_pixel_identical':True,'cells_per_new_panel':754,
            'new_fields':['TFvelo','TopicVelo'],'existing_17_scores_unchanged':True,
            'new_scores':{k:v['GDC'] for k,v in rr.items()},'TopicVelo_score_is_transition_proxy':True,'embedded_raster_images':0}
    (B/'Provenance/ed16_method_completion.json').write_text(json.dumps(report,indent=2))
    refpath=B/'Provenance/ed16_reference_revision.json'
    refreport=json.loads(refpath.read_text())
    refreport.update(panel_b_pixel_identical=False,model_rerun=True,scores_recomputed=True,
                     umap_schematic_preserved=True,method_completion='Provenance/ed16_method_completion.json')
    refpath.write_text(json.dumps(refreport,indent=2))
    print(json.dumps(report,indent=2))
if __name__=='__main__':render()
