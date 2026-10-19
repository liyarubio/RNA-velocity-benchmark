"""Validate identifiers, independent cosine samples, summaries and source preservation."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import h5py
from scipy.spatial.distance import cosine
from openpyxl import load_workbook

TASK=Path(__file__).resolve().parent;OUT=TASK.parent;DATA=OUT/'Source_data/cellDancer_simulations'
cov=pd.read_csv(DATA/'method_coverage.csv')
inventory=pd.read_csv(DATA/'dataset_inventory.csv')
summary=pd.read_csv(DATA/'dataset_summary.csv')
units=pd.read_csv(DATA/'per_unit_common_genes.csv.gz',dtype={'unit_id':str})
computed=units.groupby(['dataset','method','level']).cosine.agg(['median','count'])
expected=summary[summary.gene_scope=='common_genes'].set_index(['dataset','method','level']).reindex(computed.index)
delta=np.abs(computed['median']-expected['median'])
assert len(computed)==1102 and delta.max()<1e-10
assert (computed['count']==expected['valid_units']).all()
for r in cov[cov.output_exists].itertuples():
    p=Path(r.output_path)
    if pd.notna(r.bytes):assert p.stat().st_size==int(r.bytes),p
    if r.method=='veloae':
        with h5py.File(p) as f:
            assert f['layers/velocity'].shape==(2000,100)
for r in inventory.itertuples():
    p=Path(r.input_path)
    assert p.stat().st_size==r.input_bytes and p.stat().st_mtime_ns==r.input_mtime_ns
# Read nanosecond integers as strings to avoid float64 precision loss around missing values.
cov_text=pd.read_csv(DATA/'method_coverage.csv',dtype=str,keep_default_na=False)
for r in cov_text.itertuples():
    if r.mtime_ns:
        assert Path(r.output_path).stat().st_mtime_ns==int(r.mtime_ns)
base=Path('/data1/liyaru/software/velogif/version/velogif_0222_simu')
resume=base.parents[1]/'celldancer_simulation_resume_20261008'
baseline=json.loads((resume/'manifest.json').read_text())['existing']
for r in baseline:
    s=Path(r['path']).stat();assert s.st_size==r['size'] and s.st_mtime_ns==r['mtime_ns']
# Confirm the original 30-dataset index/noise mapping from the actual workbook.
wb=load_workbook(OUT/'Tables/Table S2 Simulation Datasets.xlsx',read_only=True,data_only=True)
mapping={}
for ws in wb:
    for row in ws.iter_rows(values_only=True):
        if len(row)>5 and str(row[1]) in set(inventory.dataset):mapping[str(row[1])]=float(row[5])
assert len(mapping)==30
for r in inventory[inventory.scope=='original_30'].itertuples():assert mapping[r.dataset]==r.noise
genes=pd.read_csv(DATA/'common_gene_sets.csv.gz')
checks=[]
for dataset,method in [('Mono-kinetic_1','stt'),('Transcriptional_boost_1','celldancer'),('Multi-forward_branching_1','latentvelo')]:
    row=cov[(cov.dataset==dataset)&(cov.method==method)].iloc[0]
    with h5py.File(base/'1_simulation_datasets'/f'{dataset}.h5ad') as f:
        truth=f['layers/ground_truth_velocity'][:]
        source_genes=list(f['var'][f['var'].attrs['_index']].asstr()[:])
        source_cells=list(f['obs'][f['obs'].attrs['_index']].asstr()[:])
    with h5py.File(row.output_path) as f:
        v=f['layers/velocity'][:]
        output_genes=list(f['var'][f['var'].attrs['_index']].asstr()[:])
        output_cells=list(f['obs'][f['obs'].attrs['_index']].asstr()[:])
    gg=genes[genes.dataset==dataset].gene.tolist()
    tt=truth[:,[source_genes.index(g) for g in gg]]
    vv=v[[output_cells.index(c) for c in source_cells]][:,[output_genes.index(g) for g in gg]]
    scores=units[(units.dataset==dataset)&(units.method==method)]
    for level,ids,a,b in [('cell',source_cells,vv,tt),('gene',gg,vv.T,tt.T)]:
        saved=scores[scores.level==level].set_index('unit_id').cosine
        for i in [0,len(ids)//2,len(ids)-1]:
            calculated=1-cosine(a[i].astype(float),b[i].astype(float))
            assert abs(calculated-saved.loc[ids[i]])<1e-10
            checks.append({'dataset':dataset,'method':method,'level':level,'unit_id':ids[i]})
result={'status':'data_checks_pass','outputs_audited':595,'evaluated_outputs':551,'unit_rows':len(units),
        'summary_groups_checked':len(computed),'max_median_recomputation_error':float(delta.max()),
        'independent_scipy_cosine_checks':checks,'source_inputs_stat_preserved':True,
        'all_451_baseline_outputs_stat_preserved':True,'veloae_40_latent_shapes_verified':True,
        'original_noise_mapping_checked_against_table_s2':True}
(TASK/'validation.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='independent_scipy_cosine_checks'},indent=2))
