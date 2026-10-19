"""Import only validated KPTracer rerun evidence into Manuscript_1010."""
from pathlib import Path
import json,shutil,hashlib
T=Path(__file__).resolve().parents[1];B=T.parent/'Manuscript_1010'
D=B/'Source_data/kptracer_completion';A=B/'Code/archive/kptracer_completion'
control=json.loads((T/'evaluation/control.json').read_text());assert control['absolute_difference']<1e-6
for method in ['tfvelo','topicvelo']:
    r=json.loads((T/'evaluation'/(method+'.json')).read_text())
    assert r['status']=='complete' and r['cells']==754 and r['reference_cell_denominator']==706
    assert r['nonzero_direction_cells']>0
assert json.loads((T/'records/topicvelo_validation.json').read_text())['status']=='pass'
assert json.loads((T/'records/tfvelo_validation.json').read_text())['status']=='pass'
D.mkdir(exist_ok=True);A.mkdir(exist_ok=True)
for p in (T/'evaluation').glob('*'):
    if p.is_file():shutil.copy2(p,D/p.name)
for p in (T/'records').glob('*.json'):shutil.copy2(p,D/p.name)
shutil.copy2(T/'records/tfvelo_gene_mapping.csv',D/'tfvelo_gene_mapping.csv')
shutil.copy2(T/'results/topicvelo/native_transition.npz',D/'native_topicvelo_transition.npz')
shutil.copy2(T/'results/topicvelo/status.json',D/'topicvelo_inference_status.json')
for p in (T/'code').glob('*.py'):shutil.copy2(p,A/p.name)
shutil.copy2(T/'logs/tfvelo.log',D/'tfvelo_fit.log')
shutil.copy2(T/'logs/topicvelo.log',D/'topicvelo_inference.log')
manifest={'run_directory':str(T),'model_outputs':{m:str(T/'results'/m/(m+'.h5ad')) for m in ['tfvelo','topicvelo']},
          'imported_files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in D.iterdir() if p.is_file()}}
(B/'Provenance/kptracer_completion_import.json').write_text(json.dumps(manifest,indent=2))
print('Validated TFvelo and native TopicVelo evidence imported.')
