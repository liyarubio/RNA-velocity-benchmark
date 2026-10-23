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

ref_idx = ratios.index('1.0')
df_ref = df[df['Ratio_j'] == 1.0].copy()
df_ref = df_ref[df_ref['Ratio_i'] != 1.0]   # exclude self-comparison

pivot = df_ref.pivot(index='Method', columns='Ratio_i', values='CosineSim')
# Reorder rows
present = [m for m in METHOD_ORDER if m in pivot.index]
pivot = pivot.loc[present]

fig, ax = plt.subplots(figsize=(6, 8))
sns.heatmap(
    pivot,
    ax=ax,
    cmap='RdYlGn',
    vmin=0, vmax=1,
    annot=True, fmt='.2f',
    linewidths=0.5,
    cbar_kws={'label': 'Mean cosine similarity to 1.0×'},
)
ax.set_title('Velocity matrix similarity to full reads (1.0×)', fontsize=13, pad=10)
ax.set_xlabel('Reads ratio', fontsize=11)
ax.set_ylabel('')
ax.set_xticklabels([f'{r}×' for r in pivot.columns], rotation=0)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'velocity_sim_to_ref.pdf'), dpi=150)
plt.close()
print("Saved velocity_sim_to_ref.pdf")
