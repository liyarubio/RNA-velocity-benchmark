"""Paired velocity consistency for target spots in full versus restricted context.

Read-only historical predictions. No models are refitted. Latent axes are never
matched by column number. Undefined directions remain missing.
"""
from pathlib import Path
import json,hashlib,importlib,traceback
import numpy as np,pandas as pd
from scipy import sparse
import scanpy as sc,scvelo as scv

B=Path('/work');OUT=Path('/task/results');OUT.mkdir(exist_ok=True)
RELEASE=B/'Manuscript_1010';D=RELEASE/'Source_data/spatial_completion'
METHODS=[('velocyto','12_velocyto'),('scVelo(sto)','2_scvelo_stochastic'),('scVelo(dyn)','3_scvelo_dynamical'),('VeloAE','14_veloAE'),('dynamo','10_Dynamo_stochastic'),('VeloVAE','13_veloVAE'),('UniTVelo','5_unitvelo'),('DeepVelo(VAE)','6_DeepveloSA'),('cellDancer','1_celldancer'),('veloVI','8_veloVI'),('LatentVelo','9_latentvelo'),('DeepVelo(GCN)','7_DeepveloGB'),('TFvelo',None),('TopicVelo','20_TopicVelo'),('STT','4_STT'),('biVI','15_biVI'),('DeepKINET','17_DeepKINET'),('SDEvelo','18_SDEvelo'),('cell2fate','16_cell2fate')]
TYPES=['Inner EPC 1','Inner EPC 2','Outer EPC 1','Outer EPC 2','P-TGC (Sec. P-TGC 1)','P-TGC (Sec. P-TGC 2)']
cells=pd.read_csv(D/'spatial_cells.csv');target=cells[cells.cell_type.isin(TYPES)].copy();ids=pd.Index(target.cell_id);assert len(ids)==608
full_directions=pd.read_csv(D/'spatial_directions.csv')
target.to_csv(OUT/'paired_target_spots.csv',index=False)
def dense(x):return x.toarray() if sparse.issparse(x) else np.asarray(x)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for x in iter(lambda:f.read(8*1024**2),b''):h.update(x)
    return h.hexdigest()
def cosine(x,y):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    nx=np.linalg.norm(x,axis=1);ny=np.linalg.norm(y,axis=1)
    valid=np.isfinite(x).all(1)&np.isfinite(y).all(1)&(nx>1e-12)&(ny>1e-12)
    c=np.full(len(x),np.nan);c[valid]=np.clip(np.sum(x[valid]*y[valid],axis=1)/(nx[valid]*ny[valid]),-1,1)
    return c,valid,nx,ny
def summarize(x):
    v=np.asarray(x);v=v[np.isfinite(v)]
    if not len(v):return {'valid_pairs':0,'median':None,'q25':None,'q75':None,'mean':None,'fraction_positive':None,'fraction_negative':None,'fraction_ge_0_9':None,'median_angle_degrees':None}
    return {'valid_pairs':len(v),'median':float(np.median(v)),'q25':float(np.quantile(v,.25)),'q75':float(np.quantile(v,.75)),'mean':float(v.mean()),'fraction_positive':float((v>0).mean()),'fraction_negative':float((v<0).mean()),'fraction_ge_0_9':float((v>=.9).mean()),'median_angle_degrees':float(np.median(np.degrees(np.arccos(v))))}
def input_path(scope,folder,name):
    if name=='TFvelo':return B/'Manuscript_2607/Figure5_sources/runs/tfvelo_20260913/results'/scope/'tfvelo_visualized.h5ad'
    return (B/'ST_XC/result/velocity/MPSTA_raw' if scope=='ST_raw' else B/'velogif_XC/real_data/result/velocity/MPSTA')/folder/'anndata.h5ad'
vg=importlib.import_module('scvelo.tools.velocity_graph');original_parallel=vg.parallelize
def compatible_parallel(*a,**kw):kw['as_array']=False;return original_parallel(*a,**kw)
vg.parallelize=compatible_parallel
records=[];all_pairs=[];summaries=[];groups=[]
for method,folder in METHODS:
    record={'method':method,'status':'running','models_refit':False,'scvelo':scv.__version__}
    try:
        paths={scope:input_path(scope,folder,method) for scope in ['ST_raw','ST_small']}
        record['sources']={s:{'path':str(p),'sha256':sha(p)} for s,p in paths.items()}
        a=sc.read_h5ad(paths['ST_raw']);b=sc.read_h5ad(paths['ST_small'])
        assert a.obs_names.is_unique and b.obs_names.is_unique and a.var_names.is_unique and b.var_names.is_unique
        assert a.n_obs==1292 and b.n_obs==608 and set(b.obs_names)==set(ids)
        ia=a.obs_names.get_indexer(ids);ib=b.obs_names.get_indexer(ids);assert (ia>=0).all() and (ib>=0).all()
        record['coordinate_translation_to_full']={}
        for scope,obj,index in [('ST_raw',a,ia),('ST_small',b,ib)]:
            delta=target[['x','y']].values-np.asarray(obj.obsm['X_spatial'])[index]
            assert np.allclose(delta,delta[0],atol=1e-8), 'Coordinates differ by more than a translation'
            record['coordinate_translation_to_full'][scope]=delta[0].tolist()
            assert np.array_equal(obj.obs.clusters.astype(str).values[index],target.cell_type.values)
        record.update(full_spots=a.n_obs,restricted_spots=b.n_obs,full_features=a.n_vars,restricted_features=b.n_vars,paired_coordinates_identical_after_translation=True)
        genes=a.var_names[a.var_names.isin(b.var_names)];ja=a.var_names.get_indexer(genes);jb=b.var_names.get_indexer(genes)
        record['shared_features_before_finite_filter']=len(genes)
        # Compare preserved processed inputs, not the currently mutable input path.
        record['aligned_input_layer_differences']={}
        if method!='VeloAE':
            for key in ['matrix','counts','spliced','unspliced','Ms','Mu']:
                if key not in a.layers or key not in b.layers:continue
                x=dense(a.layers[key][ia][:,ja]).astype(float);y=dense(b.layers[key][ib][:,jb]).astype(float)
                record['aligned_input_layer_differences'][key]={'exact_equal':bool(np.array_equal(x,y)),'relative_l2_difference':float(np.linalg.norm(x-y)/max(np.linalg.norm(x),1e-12)),'max_abs_difference':float(np.nanmax(abs(x-y)))}
        if method=='VeloAE':
            gene_cos=np.full(608,np.nan);gene_valid=np.zeros(608,dtype=bool);gene_n1=gene_n2=np.full(608,np.nan)
            record.update(gene_comparison='unavailable_independently_learned_latent_axes',shared_genes_used=0)
        else:
            if method=='TopicVelo':
                root=B/'spatial_topicvelo_retry_20261010/results/all_topics_qc'
                ta=sparse.load_npz(root/'ST_raw/native_transition.npz');tb=sparse.load_npz(root/'ST_small/native_transition.npz')
                # Match T indices to the exact source H5AD ordering, then fix a
                # single full-context Ms reference for both transition operators.
                for scope,obj in [('ST_raw',a),('ST_small',b)]:
                    saved=sc.read_h5ad(root/scope/'topicvelo.h5ad')
                    assert saved.obs_names.equals(obj.obs_names);del saved
                ref=dense(a.layers['Ms'][:,ja]);small_ref=ref[a.obs_names.get_indexer(b.obs_names)]
                x=(ta@ref-ref)[ia];y=(tb@small_ref-small_ref)[ib]
                record.update(gene_comparison='native_transition_displacement_in_fixed_full_Ms',topic_counts=[11,5],topic_thresholds_and_topic_counts_differ=True)
            else:
                x=dense(a.layers['velocity'][ia][:,ja]);y=dense(b.layers['velocity'][ib][:,jb])
                record['gene_comparison']='gene_velocity_on_shared_gene_names'
            keep=np.isfinite(x).all(0)&np.isfinite(y).all(0)
            if method=='TFvelo':
                keep &= np.isfinite(a.var.fit_scaling_y.values[ja])&np.isfinite(b.var.fit_scaling_y.values[jb])&(a.var.fit_scaling_y.values[ja]>0)&(b.var.fit_scaling_y.values[jb]>0)
            assert keep.any()
            record['shared_genes_used']=int(keep.sum())
            pd.DataFrame({'gene':genes,'used':keep}).to_csv(OUT/(method.replace('/','_')+'_shared_genes.csv'),index=False)
            gene_cos,gene_valid,gene_n1,gene_n2=cosine(x[:,keep],y[:,keep])
        if method=='TopicVelo':
            dr=pd.read_csv(D/'topicvelo_retry/ST_raw_directions.csv').set_index('cell_id').loc[ids]
            ds=pd.read_csv(D/'topicvelo_retry/ST_small_directions.csv').set_index('cell_id').loc[ids]
            v1=dr[['vx','vy']].values;v2=ds[['vx','vy']].values
            record['spatial_projection']='saved_native_transition_projection_in_fixed_spatial_coordinates'
        else:
            dr=full_directions[full_directions.method==method].set_index('cell_id').loc[ids];v1=dr[['vx','vy']].values
            if 'velocity_spatial' in b.obsm:
                v2=np.asarray(b.obsm['velocity_spatial'])[ib];record['restricted_spatial_projection']='saved_velocity_spatial'
            else:
                if 'velocity_graph' not in b.uns:
                    np.random.seed(2024);scv.tl.velocity_graph(b,vkey='velocity',n_jobs=2)
                    record['restricted_spatial_projection']='graph_recomputed_from_saved_velocity_no_model_refit'
                else:record['restricted_spatial_projection']='saved_velocity_graph_projected_no_model_refit'
                scv.tl.velocity_embedding(b,basis='spatial',vkey='velocity',autoscale=False)
                v2=np.asarray(b.obsm['velocity_spatial'])[ib]
            record['full_spatial_projection']='exact_vectors_displayed_in_ED15'
        spatial_cos,spatial_valid,n1,n2=cosine(v1,v2)
        record.update(status='complete',gene_valid_pairs=int(gene_valid.sum()),spatial_valid_pairs=int(spatial_valid.sum()))
        table=target.copy();table['method']=method;table['gene_cosine']=gene_cos;table['spatial_cosine']=spatial_cos;table['gene_valid']=gene_valid;table['spatial_valid']=spatial_valid
        table['full_vx']=v1[:,0];table['full_vy']=v1[:,1];table['restricted_vx']=v2[:,0];table['restricted_vy']=v2[:,1]
        table['gene_full_norm']=gene_n1;table['gene_restricted_norm']=gene_n2;table['spatial_full_norm']=n1;table['spatial_restricted_norm']=n2
        table.to_csv(OUT/(method.replace('/','_')+'_paired_spots.csv'),index=False);all_pairs.append(table)
        for metric in ['gene','spatial']:
            row={'method':method,'metric':metric,'paired_spots':608,'shared_genes_used':record['shared_genes_used'],'interpretation':record['gene_comparison'] if metric=='gene' else 'fixed_spatial_direction'};row.update(summarize(table[metric+'_cosine']));summaries.append(row)
            for kind,part in table.groupby('cell_type'):
                row={'method':method,'metric':metric,'cell_type':kind,'paired_spots':len(part)};row.update(summarize(part[metric+'_cosine']));groups.append(row)
        del a,b
    except Exception:record.update(status='failed',error=traceback.format_exc())
    records.append(record)
    (OUT/'audit.json').write_text(json.dumps(records,indent=2)+'\n')
    pd.DataFrame(summaries).to_csv(OUT/'summary.csv',index=False);pd.DataFrame(groups).to_csv(OUT/'by_cell_type.csv',index=False)
    print(json.dumps({k:record[k] for k in ['method','status','full_features','restricted_features','shared_genes_used','gene_valid_pairs','spatial_valid_pairs'] if k in record}),flush=True)
    if record['status']=='failed':print(record['error'],flush=True)
pd.concat(all_pairs,ignore_index=True).to_csv(OUT/'paired_spots.csv',index=False)
assert len(records)==19 and all(r['status']=='complete' for r in records)
