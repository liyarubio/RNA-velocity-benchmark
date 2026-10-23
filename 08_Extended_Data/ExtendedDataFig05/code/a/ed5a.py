"""Extended Data 5a topology comparisons. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next((p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file()))
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    pair = read('Extended_Data_Fig06.xlsx', 'Paired cell-level values')
    pair['topology'] = pair.dataset.str.replace('_seed\\d+$', '', regex=True)
    topologies = pair.topology.drop_duplicates().tolist()
    assert len(topologies) == 14 and pair.dataset.nunique() == 42 and len(pair) == 698
    methods = pair.method.drop_duplicates().tolist()
    colors = dict(zip(methods, plt.get_cmap('tab20').colors))
    fig, axes = plt.subplots(3, 5, figsize=(15, 9), constrained_layout=True)
    for ax, top in zip(axes.flat, topologies):
        for m in methods:
            z = pair[(pair.topology == top) & (pair.method == m)][['high', 'low']]
            means = z.mean()
            sd = z.std()
            if len(z):
                ax.errorbar(means.high, means.low, xerr=sd.high if np.isfinite(sd.high) else None, yerr=sd.low if np.isfinite(sd.low) else None, fmt='o', ms=3, lw=0.5, color=colors[m], label=m)
        ax.set(title=top, xlabel='High-dimensional cosine', ylabel='PCA2 cosine')
    axes.flat[-1].axis('off')
    handles, labels = axes.flat[0].get_legend_handles_labels()
    axes.flat[-1].legend(handles, labels, fontsize=6, ncol=2, frameon=False)
    save(fig, 'ED05a', 'Extended_Data_Fig06.xlsx / Paired cell-level values; matched means and sample SD', out=out, log=log)
if __name__ == '__main__':
    run_panel(render)
