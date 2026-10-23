"""Independently audit latest saved summaries; never fit or impute a model."""
from pathlib import Path
import argparse,json,sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from release_paths import source_path

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('/tmp/rna-velocity-1021'))
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    checks=[]
    def check(name,condition,details=None):
        checks.append(dict(check=name,passed=bool(condition),details=details))
    def read(rel):return pd.read_csv(source_path('data/'+rel))
    def equal(a,b):return np.allclose(a,b,atol=1e-12,rtol=0,equal_nan=True)
    plot=read('cellDancer_noise_20261020/plot_data.csv')
    original=read('cellDancer_current_30/dataset_summary.csv')
    original=original[original.gene_scope.eq('common_genes')]
    key=['dataset','method','level'];a=plot.set_index(key);b=original.set_index(key)
    check('noise_grid',len(a)==840 and a.index.is_unique and plot.dataset.nunique()==30 and plot.method.nunique()==14)
    check('noise_saved_medians',equal(a.loc[b.index,'median'],b['median']))
    check('noise_valid_units',(a.loc[b.index,'valid_units']==b.valid_units).all())
    check('noise_missing_slots',plot['median'].isna().sum()==18 and a.loc[a.index.difference(b.index),'median'].isna().all())
    group=['family','level','method'];means=read('cellDancer_noise_20261020/noise_means.csv').set_index(group)
    expected=original.groupby(group)['median'].agg(['mean','count']).reindex(means.index)
    check('noise_cross_condition_means',len(means)==84 and equal(means['mean'],expected['mean']))
    check('noise_coverage',(means.n_available==expected['count']).all() and means.coverage.value_counts().to_dict()=={'10/10':82,'1/10':2})
    check('noise_limited_groups',set(means[means.n_available.lt(10)].index)=={('Multi-forward_branching','cell','stt'),('Multi-forward_branching','gene','stt')})
    book=source_path('source_data/Extended_Data_Fig05.xlsx')
    ws=pd.read_excel(book,sheet_name='c Noise-series plot').set_index(key)
    check('noise_workbook_values',a.index.equals(ws.index) and equal(a['median'],ws['median']))
    e=json.loads(source_path('data/topovelo_historical_20261021/derived/evaluation.json').read_text())
    for r in e['inputs']:
        terms=[]
        for cluster,v in r['inner'].items():
            edges=[x for x in r['edges'] if cluster in (x['source'],x['target']) and x['n']]
            if v['n'] and edges:terms.append(v['n']/r['denominator']*v['mean']*np.mean([np.sign(x['median']) for x in edges]))
        check('TopoVelo_'+r['dataset']+'_GDC',abs(sum(terms)-r['GDC'])<1e-12,dict(recomputed=sum(terms),saved=r['GDC'],denominator=r['denominator']))
    paired=read('spatial_context_consistency/paired_spots.csv');t=paired[paired.method.eq('TopoVelo')]
    check('paired_population',len(paired)==20*608 and not paired.duplicated(['method','cell_id']).any())
    check('TopoVelo_missing_pairs',len(t)==608 and t.gene_cosine.notna().sum()==606 and t.spatial_cosine.notna().sum()==606 and set(t.loc[t.spatial_cosine.isna(),'cell_id'])==set(e['undefined_pair_ids']))
    shared=read('topovelo_historical_20261021/derived/TopoVelo_shared_genes.csv')
    check('TopoVelo_shared_genes',len(shared)==47 and shared.used.all())
    display=read('figure5_1021/KPTracer_GDC_with_PhyloVelo.csv')
    phy=json.loads(source_path('data/kptracer_phylovelo_20261021/evaluation.json').read_text())
    check('PhyloVelo_no_comparable_rank',phy['GDC_comparable'] is None and phy['GDC_rank'] is None and phy['GDC_adapted_forward']==1.0)
    check('PhyloVelo_direction_contract',phy['plot_vector_key']=='velocity_umap_backward' and phy['metric_vector_key']=='velocity_umap_forward')
    check('PhyloVelo_display_missing_rank',display.loc[display.method.eq('PhyloVelo'),'rank'].isna().all())
    f6=source_path('source_data/Figure06.xlsx');summary=pd.read_excel(f6,sheet_name='6c Summary')
    check('frozen_composite_19_configurations',len(summary)==19 and 'PhyloVelo' not in summary.Method.values)
    report=dict(status='pass' if all(x['passed'] for x in checks) else 'fail',checks=checks,number_of_checks=len(checks),models_rerun=False,scope='Saved-summary recalculation. TopoVelo gene-level input audit and model training are not rerun. PhyloVelo H5AD audit is a separate command.')
    (args.out/'latest_results_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));return report['status']!='pass'
if __name__=='__main__':raise SystemExit(main())
