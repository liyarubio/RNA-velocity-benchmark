"""Historical read-depth plot only; use the current figure6b.py for the released panel."""
from pathlib import Path
import argparse,os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input',type=Path,required=True,help='Historical velocity_cosine_sim.csv')
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True);out_dir=str(args.out)
df=pd.read_csv(args.input)
ratios=['0.2','0.4','0.6','0.8','1.0']
ref_idx=ratios.index('1.0')
pairwise={m:g.pivot(index='Ratio_i',columns='Ratio_j',values='CosineSim').reindex(index=[float(r) for r in ratios],columns=[float(r) for r in ratios]).to_numpy(float) for m,g in df.groupby('Method')}
METHOD_ORDER = [
    'velocyto', 'scVelo(sto)', 'scVelo(dyn)', 'VeloAE', 'dynamo',
    'VeloVAE', 'MultiVelo', 'UniTVelo', 'DeepVelo(VAE)', 'cellDancer',
    'PhyloVelo', 'veloVI', 'LatentVelo', 'DeepVelo(GCN)', 'TFvelo', 'STT',
    'biVI', 'SDEvelo', 'DeepKINET',
]
methods_present=[m for m in METHOD_ORDER if m in pairwise]

methods_present = [m for m in METHOD_ORDER if m in pairwise]
n = len(methods_present)
ncols = 4
nrows = (n + ncols - 1) // ncols

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 3.2))
axes = axes.flatten()

ratio_labels = [f'{r}×' for r in ratios]

for idx, mname in enumerate(methods_present):
    ax = axes[idx]
    mat = pairwise[mname]
    sns.heatmap(
        mat,
        ax=ax,
        cmap='RdYlGn',
        vmin=-0.2, vmax=1,
        annot=True, fmt='.2f',
        linewidths=0.5,
        xticklabels=ratio_labels,
        yticklabels=ratio_labels,
        cbar=False,
    )
    ax.set_title(mname, fontsize=10, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(labelsize=8)

# Hide unused subplots
for idx in range(len(methods_present), len(axes)):
    axes[idx].axis('off')

# Add colorbar
sm = plt.cm.ScalarMappable(cmap='RdYlGn', norm=plt.Normalize(vmin=-0.2, vmax=1))
sm.set_array([])
fig.colorbar(sm, ax=axes[:len(methods_present)], shrink=0.6, label='Mean cosine similarity')

plt.suptitle('Pairwise velocity matrix similarity across reads ratios', fontsize=13, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'velocity_sim_pairwise_all.pdf'),
            dpi=150, bbox_inches='tight')
plt.close()
print("Saved velocity_sim_pairwise_all.pdf")
