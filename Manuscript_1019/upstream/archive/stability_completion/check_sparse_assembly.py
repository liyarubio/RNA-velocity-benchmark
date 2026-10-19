"""Compare literal historical transition assembly with the sparse replacement."""
import ast
import inspect
import json
from pathlib import Path
import warnings
import numpy as np
from scipy.sparse import csr_matrix
from topicvelo.transition_matrix import Combined_Topics_Transitions

HERE=Path(__file__).resolve().parent
tree=ast.parse((HERE/'complete_topicvelo.py').read_text())
replacement=next(n.value for n in ast.walk(tree) if isinstance(n,ast.Constant)
                 and isinstance(n.value,str) and n.value.startswith('        combined_TM = csr_matrix'))
original=inspect.getsource(Combined_Topics_Transitions)
start=original.index('        #construct the combined transition matrix')
end=original.index('        #save transition matrix',start)
old='\n'.join(x[8:] for x in original[start:end].splitlines())
new='\n'.join(x[8:] for x in replacement.splitlines())
rng=np.random.default_rng(2024)
checks=[]
for n,k in [(40,1),(40,3),(200,1)]:
    inds=[np.sort(rng.choice(n,int(n*.65),replace=False)) for _ in range(k)]
    mats=[]
    for idx in inds:
        a=rng.random((len(idx),len(idx)));a[a<.8]=0;a+=np.eye(len(idx));a/=a.sum(1)[:,None]
        mats.append(csr_matrix(a))
    weights=np.zeros((n,k))
    for i,idx in enumerate(inds):weights[idx,i]=rng.random(len(idx))
    norm=weights.sum(1);weights=np.divide(weights,norm[:,None],out=np.zeros_like(weights),where=norm[:,None]!=0)
    class Data: pass
    adata=Data();adata.n_obs=n;adata.obsp={'connectivities':csr_matrix(np.ones((n,n))-np.eye(n))}
    common=dict(np=np,csr_matrix=csr_matrix,n=n,topics=list(range(k)),TMs=mats,ttc_indices=inds,
                cells_weights_for_tm=weights,adata=adata)
    a=dict(common);b=dict(common)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore');exec(old,a);exec(new,b)
    A=csr_matrix(a['combined_TM']).toarray();B=csr_matrix(b['combined_TM']).toarray()
    # Historical code rounds probabilities to six decimal places. Different
    # sparse summation order can cross one rounding boundary. Bound the
    # difference by one six-decimal rounding unit plus float32 normalization.
    np.testing.assert_allclose(A,B,rtol=0,atol=1.1e-6)
    checks.append({'cells':n,'topics':k,'max_absolute_difference':float(np.max(abs(A-B)))})
# Same math, but preserve singleton gene dimension in the old GLS routine.
import importlib,types
opt=importlib.import_module('scvelo.tools.optimization');f=opt.leastsq_generalized
g=dict(f.__globals__);g['make_dense']=lambda x:x.toarray() if hasattr(x,'toarray') else np.asarray(x)
fixed=types.FunctionType(f.__code__,g,f.__name__,f.__defaults__,f.__closure__)
for genes in [1,3]:
    x=rng.random((30,genes))+.1
    z=fixed(x,2*x,3*x,6*x)
    np.testing.assert_allclose(z[2],2,atol=1e-6)
    if genes>1:
        for aa,bb in zip(z,f(x,2*x,3*x,6*x)):np.testing.assert_array_equal(aa,bb)
report={'status':'pass','assembly_absolute_tolerance':1.1e-6,'sparse_assembly_tests':checks,'singleton_GLS_known_slope':2,'multigene_GLS_unchanged':True}
(HERE/'compatibility_validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
