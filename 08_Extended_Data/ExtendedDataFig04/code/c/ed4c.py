"""Extended Data 4c full noise series. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    from matplotlib.colors import LinearSegmentedColormap
    plot = read('Extended_Data_Fig05.xlsx', 'c Noise-series plot')
    means = read('Extended_Data_Fig05.xlsx', 'c Noise means')
    methods = plot.method.drop_duplicates().tolist()
    families = plot.family.drop_duplicates().tolist()
    labels = plot.drop_duplicates('method').set_index('method').display_name
    cmap = LinearSegmentedColormap.from_list('light_blue_A', plt.get_cmap('Blues')(np.linspace(0.02, 0.68, 256)))
    cmap.set_bad('#dddddd')
    fig, axes = plt.subplots(2, 3, figsize=(17, 8), constrained_layout=True)
    for r, level in enumerate(['cell', 'gene']):
        for c, family in enumerate(families):
            q = plot[(plot.level == level) & (plot.family == family)].pivot(index='method', columns='noise', values='median').reindex(methods)
            q['Mean'] = means[(means.level == level) & (means.family == family)].set_index('method')['mean']
            ax = axes[r, c]
            im = ax.imshow(q, cmap=cmap, vmin=np.nanmin(q), vmax=np.nanmax(q), aspect='auto')
            ax.set_yticks(range(len(q)), [labels[m] for m in q.index])
            ax.set_xticks(range(len(q.columns)), [str(x) for x in q.columns], rotation=45)
            ax.set_title(f'{family}: {level}')
            for i in range(len(q)):
                for j in range(len(q.columns)):
                    v = q.iloc[i, j]
                    ax.text(j, i, f'{v:.2f}' if np.isfinite(v) else 'x', ha='center', va='center', fontsize=6)
            fig.colorbar(im, ax=ax, shrink=0.7)
    save(fig, 'ED04c', 'Extended_Data_Fig05.xlsx / c Noise-series plot; c Noise means; independent color scale per heatmap', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
