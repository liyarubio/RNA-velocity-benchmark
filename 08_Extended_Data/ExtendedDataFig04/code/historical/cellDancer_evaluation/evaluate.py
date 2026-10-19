"""Audit saved outputs and evaluate against source truth; never train or modify inputs.
Run: /data1/liyaru/miniconda/bin/python evaluate.py
"""
from pathlib import Path
import csv, gzip, hashlib, json
import numpy as np
import pandas as pd
import h5py
from scipy import sparse

TASK = Path(__file__).resolve().parent
OUT = TASK.parent
BASE = Path('/data1/liyaru/software/velogif/version/velogif_0222_simu')
RESUME = BASE.parents[1] / 'celldancer_simulation_resume_20261008'
DATA = OUT / 'Source_data/cellDancer_simulations'
DATA.mkdir(parents=True, exist_ok=True)
METHODS = ['velocyto','scvelo.sto','scvelo.dyn','dynamo','velovae','unitvelo',
           'deepvelo_vae','deepvelo_gcn','velovi','sdevelo','celldancer','latentvelo','stt','bivi']
EXCLUDED = {'veloae':'latent coordinates are not gene-space velocities',
            'tfvelo':'synthetic identifiers have no verified TF mapping; four historical outputs not validated',
            'deepkinet':'integer count input unavailable', 'cell2fate':'integer count input unavailable',
            'topicvelo':'native count-based inference not validated on continuous input',
            'topovelo':'physical spatial coordinates unavailable',
            'multivelo':'chromatin accessibility unavailable', 'phylovelo':'lineage information unavailable'}
FAMILIES = ['Mono-kinetic','Transcriptional_boost','Multi-forward_branching','Multi-backward_branching']

def names(f, key):
    g=f[key]
    return g[g.attrs.get('_index','_index')].asstr()[:].astype(str)

def matrix(node):
    if isinstance(node,h5py.Dataset): return node[:].astype(np.float64)
    enc=node.attrs['encoding-type']
    constructor=sparse.csr_matrix if enc=='csr_matrix' else sparse.csc_matrix
    assert enc in ('csr_matrix','csc_matrix'), enc
    return constructor((node['data'][:],node['indices'][:],node['indptr'][:]),shape=tuple(node.attrs['shape'])).toarray().astype(np.float64)

def cosine(a,b,axis):
    # Invalid vectors, including zero norms, are undefined rather than zero.
    assert a.shape == b.shape
    den=np.sqrt(np.sum(a*a,axis=axis)*np.sum(b*b,axis=axis))
    result=np.full(den.shape,np.nan)
    valid=np.isfinite(a).all(axis)&np.isfinite(b).all(axis)&np.isfinite(den)&(den>0)
    result[valid]=np.sum(a*b,axis=axis)[valid]/den[valid]
    return np.clip(result,-1,1)

def digest(a): return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def stats(a):
    a=a[np.isfinite(a)]
    return {'valid_units':len(a),'median':float(np.median(a)) if len(a) else None,
            'mean':float(np.mean(a)) if len(a) else None}

summary, coverage, audits, gene_sets = [], [], [], []
sample = {}
with gzip.open(DATA/'per_unit_common_genes.csv.gz','wt',newline='') as fh:
    writer=csv.writer(fh)
    writer.writerow(['dataset','method','level','unit_id','cosine'])
    for family in FAMILIES:
      for index in range(1,11):
        dataset=f'{family}_{index}'
        source=BASE/'1_simulation_datasets'/f'{dataset}.h5ad'
        with h5py.File(source) as f:
            cells,genes=names(f,'obs'),names(f,'var')
            truth=matrix(f['layers/ground_truth_velocity'])
        assert len(set(cells))==len(cells) and len(set(genes))==len(genes)
        assert truth.shape==(2000,1000) and np.isfinite(truth).all()
        lookup_c={v:i for i,v in enumerate(cells)}
        lookup_g={v:i for i,v in enumerate(genes)}
        meta={'dataset':dataset,'family':family,'series_index':index,
              'noise':round(index*.2,1) if family!=FAMILIES[-1] else None,
              'scope':'original_30' if family!=FAMILIES[-1] else 'backward_extension_10'}
        fitted={}
        for method in METHODS+list(EXCLUDED):
            old=BASE/'result_simu'/dataset/f'{method}.h5ad'
            new=RESUME/'results'/dataset/f'{method}.h5ad'
            assert not (old.exists() and new.exists()), (dataset,method,'duplicate')
            path=new if new.exists() else old
            row={**meta,'method':method,'output_exists':path.exists(),'output_path':str(path) if path.exists() else '',
                 'provenance':'recovery_20261008' if new.exists() else 'historical',
                 'exploratory':method=='bivi'}
            if method in EXCLUDED:
                row.update(status='not_evaluated',reason=EXCLUDED[method]);coverage.append(row);continue
            if not path.exists():
                state=RESUME/'state'/f'{dataset}__{method}.json'
                state_data=json.loads(state.read_text()) if state.exists() else {}
                row.update(status=state_data.get('status','missing'),reason=state_data.get('error','No saved velocity output'))
                coverage.append(row);continue
            with h5py.File(path) as f:
                oc,og=names(f,'obs'),names(f,'var')
                v=matrix(f['layers/velocity'])
            assert len(set(oc))==len(oc) and len(set(og))==len(og)
            assert set(oc)==set(cells), (dataset,method,'cell mismatch')
            assert set(og)<=set(genes), (dataset,method,'gene mismatch')
            oi={name:i for i,name in enumerate(oc)}
            v=v[[oi[name] for name in cells],:]
            gi=np.array([lookup_g[name] for name in og])
            finite=np.isfinite(v).all(0)
            gi,v=gi[finite],v[:,finite]
            assert len(gi)>0
            fitted[method]=(gi,v)
            row.update({'status':'exploratory' if method=='bivi' else 'evaluated', 'reason':'',
                       'cells':len(cells),'output_genes':len(og),'finite_genes':len(gi),
                       'bytes':path.stat().st_size,'mtime_ns':path.stat().st_mtime_ns,
                       'aligned_finite_velocity_sha256':digest(v),
                       'protocol_note':('Relaxed filtering in backward recovery; historical Dynamo uses original filtering' if method=='dynamo' else
                                        'Historical outputs retain selected genes; recovery outputs retain all genes' if method=='stt' else
                                        'Continuous input violates count observation domain; software execution only' if method=='bivi' else 'Saved method-specific preprocessing')})
            coverage.append(row)
        # Same cells and same finite gene coordinates across displayed methods in this dataset.
        # biVI does not define the comparison cohort. Its coordinates must cover that cohort.
        common=sorted(set.intersection(*(set(gi) for m,(gi,v) in fitted.items() if m!='bivi')))
        assert len(common)>0, dataset
        common=np.array(common)
        for g in common:gene_sets.append({'dataset':dataset,'gene':genes[g]})
        audits.append({**meta,'input_path':str(source),'input_bytes':source.stat().st_size,
                       'input_mtime_ns':source.stat().st_mtime_ns,'truth_sha256':digest(truth),
                       'cells':len(cells),'source_genes':len(genes),'common_genes':len(common),
                       'comparison_methods':[m for m in fitted if m!='bivi'],
                       'noise_mapping_evidence':'Table S2 cellDancer sheet' if family!=FAMILIES[-1] else 'Unverified: show series index only'})
        for method,(gi,v) in fitted.items():
            positions={g:i for i,g in enumerate(gi)}
            assert set(common)<=set(gi)
            for scope in ('method_retained_genes','common_genes'):
                gg=gi if scope=='method_retained_genes' else common
                vv=v if scope=='method_retained_genes' else v[:,[positions[g] for g in common]]
                tt=truth[:,gg]
                for level,axis,unit_ids in [('cell',1,cells),('gene',0,genes[gg])]:
                    scores=cosine(vv,tt,axis)
                    summary.append({**meta,'method':method,'gene_scope':scope,'level':level,
                                    'cells':len(cells),'genes':len(gg),'total_units':len(scores),
                                    **stats(scores),'exploratory':method=='bivi'})
                    if scope=='common_genes':
                        writer.writerows((dataset,method,level,u,'' if not np.isfinite(x) else format(x,'.12g')) for u,x in zip(unit_ids,scores))
                        if index==1 and family!=FAMILIES[-1]:sample[f'{dataset}__{method}__{level}']=scores
        print(dataset, 'evaluated',len(fitted),'common genes',len(common),flush=True)
pd.DataFrame(summary).to_csv(DATA/'dataset_summary.csv',index=False)
pd.DataFrame([{**r, 'mtime_ns':str(r['mtime_ns']) if 'mtime_ns' in r else ''} for r in coverage]).to_csv(DATA/'method_coverage.csv',index=False)
pd.DataFrame(gene_sets).to_csv(DATA/'common_gene_sets.csv.gz',index=False,compression='gzip')
pd.DataFrame(audits).drop(columns=['comparison_methods']).to_csv(DATA/'dataset_inventory.csv',index=False)
(DATA/'alignment_audit.json').write_text(json.dumps(audits,indent=2))
np.savez_compressed(TASK/'noise_02_distributions.npz',**sample)
(DATA/'evaluation_protocol.json').write_text(json.dumps({
    'truth':'source input layers/ground_truth_velocity, never transformed result-file truth',
    'estimate':'saved layers/velocity, cell/gene identifiers aligned, no scale recalibration',
    'plotted_scope':'same finite genes across available non-exploratory methods within each dataset; biVI evaluated on that set separately',
    'sensitivity_scope':'method_retained_genes in dataset_summary.csv',
    'unit':'cell across genes or gene across cells; median of finite cosines per dataset',
    'invalid':'zero norm or nonfinite vector -> missing; never assign zero',
    'replication':'one saved simulated dataset per scenario/noise condition; cells and genes are not independent simulation replicates',
    'preprocessing':'saved protocols differ across methods; coordinate alignment is not a calibration of physical velocity units',
    'original_scope':'30 datasets: 3 scenarios x 10 noise levels; FigS2 uses noise 0.2',
    'extension':'10 backward-branching datasets; noise mapping unverified and displayed as series index',
    'excluded':EXCLUDED,'exploratory':['bivi'],
    'main_figures_and_rankings':'unchanged; these results are not added to the dyngen or composite benchmark',
    'source_code':str(Path(__file__).resolve())},indent=2))
print('DONE',len(summary),'summary rows',len(coverage),'coverage rows')
