"""Merge independently validated reads results, retaining all historical values."""
from pathlib import Path
import json
import hashlib
import shutil
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
R=HERE.parent
M=R.parent
RUN=Path('/data1/liyaru/software/velogif/velogif_0315/stability_completion_20261008')
source=M/'Figure5_sources/data/Figure5c_reads_velocity_cosine_similarity.csv'
before=HERE/'before/Figure5_sources/data/Figure5c_reads_velocity_cosine_similarity.csv'
old=pd.read_csv(before,float_precision='round_trip')
assert len(old)==475
new=[];protocol=[];coverage=[]
for method in ['cell2fate','topicvelo']:
    folder=RUN/'evaluation'/method
    report=json.loads((folder/'validation.json').read_text())
    assert report['status']=='complete'
    frame=pd.read_csv(folder/'pairwise_cosine.csv',float_precision='round_trip')
    assert len(frame)==25 and frame.CosineSim.between(-1,1).all()
    square=frame.pivot(index='Ratio_i',columns='Ratio_j',values='CosineSim')
    np.testing.assert_allclose(square,square.T,atol=1e-12,rtol=0)
    np.testing.assert_allclose(np.diag(square),1,atol=1e-12,rtol=0)
    for key,n in report['valid_cells_by_comparison'].items():
        assert n>0
    new.append(frame)
    protocol.append({'Method':report['method'],'Metric':report['metric'],
                     'Common_cells':report['common_cells'],'Common_genes':report['common_genes'],
                     'Mean_stability':report['mean_stability'],'Run_directory':str(RUN),
                     'Validation':str(folder/'validation.json')})
    for i,ri in enumerate([.2,.4,.6,.8,1.]):
        for j,rj in enumerate([.2,.4,.6,.8,1.]):
            valid=report['valid_cells_by_comparison'][str(i*5+j)]
            coverage.append({'Method':report['method'],'Ratio_i':ri,'Ratio_j':rj,
                'common_cells':report['common_cells'],'common_genes':report['common_genes'],
                'valid_cells':valid,'invalid_cells':report['common_cells']-valid,
                'valid_fraction':valid/report['common_cells'],'metric':report['metric']})
    shutil.copy2(folder/'validation.json',HERE/(method+'_validation.json'))
    shutil.copy2(folder/'pairwise_cosine.csv',HERE/(method+'_pairwise_cosine.csv'))
    shutil.copy2(folder/'per_cell_cosine.csv.gz',HERE/(method+'_per_cell_cosine.csv.gz'))
shutil.copy2(RUN/'evaluation/topicvelo/native_transition_row_cosine.csv',HERE/'topicvelo_native_transition_row_cosine.csv')
combined=pd.concat([old]+new,ignore_index=True)
assert len(combined)==525 and combined.Method.nunique()==21
pd.testing.assert_frame_equal(combined.iloc[:475],old)
# Preserve the historical 475-row block verbatim, including original precision.
payload=before.read_bytes()
if not payload.endswith(b'\n'):payload+=b'\n'
source.write_bytes(payload+pd.concat(new,ignore_index=True).to_csv(index=False,header=False).encode())
actual=pd.read_csv(source,float_precision='round_trip')
pd.testing.assert_frame_equal(actual.iloc[:475],old)
for method, frame in zip(['cell2fate','TopicVelo'],new):
    pd.testing.assert_frame_equal(actual.loc[actual.Method==method].reset_index(drop=True),frame)
plot=actual.loc[actual.Ratio_j==1].copy()
plot.to_csv(M/'Figure5_sources/reproduced/Figure5c_reads_velocity_stability_plot_data.csv',index=False)
pd.DataFrame(protocol).to_csv(HERE/'added_methods_protocol.csv',index=False)
pd.DataFrame(protocol).to_csv(M/'Figure5_sources/data/Figure5c_added_methods_protocol.csv',index=False)
pd.DataFrame(coverage).to_csv(HERE/'added_methods_coverage.csv',index=False)
pd.DataFrame(coverage).to_csv(M/'Figure5_sources/data/Figure5c_added_methods_coverage.csv',index=False)
manifest=M/'Figure5_sources/data/data_manifest.csv'
d=pd.read_csv(HERE/'before/Figure5_sources/data/data_manifest.csv')
row=d.Local_file.eq(source.name) if 'Local_file' in d else d.iloc[:,1].eq(source.name)
for col in d.columns:
    if col=='Local_MD5':d.loc[row,col]=hashlib.md5(source.read_bytes()).hexdigest()
    if 'size' in col.lower():d.loc[row,col]=source.stat().st_size
    if col=='Meaning':d.loc[row,col]='21-method 5x5 reads cosine; TopicVelo uses fixed-reference transition projection'
    if col=='Original_source':d.loc[row,col]=d.loc[row,col].astype(str)+'; '+str(HERE/'added_methods_protocol.csv')
    if col=='Transformation':d.loc[row,col]='Historical 475 rows preserved; 50 validated rows appended; Original_MD5 refers to historical source'
d.to_csv(manifest,index=False)
result={'status':'pass','historical_rows_preserved':475,'added_rows':50,'methods_figure5':21,
        'methods_figure6':19,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'protocol':protocol,'topicvelo_interpretation':'Native integrated transition displacement in fixed full-depth common-gene Ms; not the base scVelo velocity. Historical one-topic stochastic configuration and topic fits retained.'}
(HERE/'merge_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
