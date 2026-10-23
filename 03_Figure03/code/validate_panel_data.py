"""Recalculate Figure 3d/f from current saved inputs; no fitting or layout rendering."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
DATA=Path(__file__).resolve().parents[1]/'data/figure3_detail'
def cosine(x,y):
    keep=np.isfinite(x)&np.isfinite(y);x=x[keep];y=y[keep]
    den=np.linalg.norm(x)*np.linalg.norm(y)
    return float(x@y/den) if den else np.nan

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('/tmp/figure3-panel-validation.json'))
    args=parser.parse_args();checks=[]
    def check(name,ok,details=None):checks.append(dict(check=name,passed=bool(ok),details=details))
    def read(name):return pd.read_csv(DATA/name)
    high=read('Figure3a_high_dim_cosine_similarity.csv');high=high.set_index(high.columns[0])
    pca=read('Figure3d_pca2_to_30_cosine_similarity.csv')
    out=[];all_x={i:[] for i in range(2,31)};all_y={i:[] for i in range(2,31)}
    for dataset,part in pca.groupby('Dataset',sort=False):
        methods=[c for c in part if c not in ['Dataset','Dim'] and 'dynamo_deterministic' not in c.lower()]
        y=high.reindex(methods)[dataset].to_numpy(float)
        for _,row in part.iterrows():
            dim=int(str(row.Dim).removeprefix('dim_'));x=row[methods].to_numpy(float)
            out.append(dict(Dataset_type=dataset.rsplit('_seed',1)[0],Dim=dim,Cosine_similarity=cosine(x,y)))
            all_x[dim].extend(x);all_y[dim].extend(y)
    actual=pd.DataFrame(out).groupby(['Dataset_type','Dim']).Cosine_similarity.mean()
    totals=pd.Series({('total',i):cosine(np.asarray(all_x[i]),np.asarray(all_y[i])) for i in range(2,31)})
    actual=pd.concat([actual,totals])
    expected=read('Figure3d_plot_data.csv').set_index(['Dataset_type','Dim']).Cosine_similarity
    check('panel_d_435_dimension_summaries',len(actual)==435 and np.allclose(actual.reindex(expected.index),expected,atol=1e-12,rtol=0,equal_nan=True))
    files={'Global':'Figure3f_global_metric.csv','Incluster':'Figure3f_incluster_metric.csv','Direction':'Figure3f_direction_metric.csv','velocity':'Figure3a_high_dim_cosine_similarity.csv','velocity_pca':'Figure3bc_low_dim_cosine_similarity.csv','trans_matrix':'Figure3f_transition_matrix_cosine_similarity.csv'}
    frames=[]
    for metric,file in files.items():
        x=read(file);x=x.rename(columns={x.columns[0]:'Method'});x=x[~x.Method.str.contains('Dynamo_deterministic',case=False)]
        frames.append(x.melt(id_vars='Method',var_name='Dataset',value_name=metric).set_index(['Method','Dataset']))
    actual=pd.concat(frames,axis=1).sort_index();expected=read('Figure3f_aligned_metrics.csv').set_index(['Method','Dataset']).sort_index()
    check('panel_f_798_aligned_units',len(actual)==798 and actual.index.equals(expected.index) and np.allclose(actual[expected.columns],expected,atol=1e-12,rtol=0,equal_nan=True))
    corr=actual.corr();saved=read('Figure3f_correlation.csv').set_index('Metric')
    check('panel_f_pairwise_Pearson_correlations',np.allclose(corr.loc[saved.index,saved.columns],saved,atol=1e-12,rtol=0,equal_nan=True))
    report=dict(status='pass' if all(c['passed'] for c in checks) else 'fail',checks=checks,models_rerun=False,R_plotting_executed=False,scope='Independent Python checks of saved d/f numerical derivations; not an R renderer execution or a formal-layout reconstruction.')
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));return report['status']!='pass'
if __name__=='__main__':raise SystemExit(main())
