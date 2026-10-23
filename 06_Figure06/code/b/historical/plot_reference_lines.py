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

METHOD_COLORS = {
    'velocyto':      '#D62728FF',
    'scVelo(sto)':   '#cab2d6',
    'scVelo(dyn)':   '#6a3d9a',
    'VeloAE':        '#33a02c',
    'dynamo':        '#FFB547FF',
    'VeloVAE':       '#b2df8a',
    'MultiVelo':     '#843C39FF',
    'UniTVelo':      '#bcbd22',
    'DeepVelo(VAE)': '#1f78b4',
    'cellDancer':    '#e377c2',
    'PhyloVelo':     '#E7CB94FF',
    'veloVI':        '#ff9896',
    'LatentVelo':    '#17BECFFF',
    'DeepVelo(GCN)': '#a6cee3',
    'TFvelo':        '#7B4173FF',
    'STT':           '#ff7f00',
    'biVI':          '#8c564b',
    'SDEvelo':       '#e6ab02',
    'DeepKINET':     '#66c2a5',
}

fig, ax = plt.subplots(figsize=(8, 5))
for mname in methods_present:
    mat = pairwise[mname]
    ref_col = mat[:, ref_idx]  # similarity of each ratio to 1.0x
    xs = [float(r) for r in ratios]
    color = METHOD_COLORS.get(mname, 'grey')
    ax.plot(xs, ref_col, marker='o', markersize=4, linewidth=0.9,
            color=color, label=mname)

ax.axhline(1.0, color='grey', linestyle='--', linewidth=0.7)
ax.set_xlabel('Reads ratio', fontsize=12)
ax.set_ylabel('Mean cosine similarity to 1.0×', fontsize=12)
ax.set_title('Velocity matrix stability across reads downsampling', fontsize=13)
ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_xlim(0.15, 1.05)
ax.legend(fontsize=7.5, ncol=2, loc='lower right', framealpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'velocity_sim_to_ref_line.pdf'), dpi=150)
plt.close()
print("Saved velocity_sim_to_ref_line.pdf")
