"""Extended Data 5b paired distributions. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    pair = read('Extended_Data_Fig06.xlsx', 'Paired cell-level values')
    methods = pair.method.drop_duplicates().tolist()
    fig, axes = plt.subplots(3, 6, figsize=(15, 8), constrained_layout=True)
    for ax, m in zip(axes.flat, methods):
        v = pair[pair.method == m].sort_values('dataset')[['high', 'low']].to_numpy(float)
        for row in v:
            ax.plot([0, 1], row, c='#bbbbbb', lw=0.35)
        for j, color in [(0, '#08306b'), (1, '#c6dbef')]:
            values = v[:, j]
            values = values[np.isfinite(values)]
            if len(np.unique(values)) > 1:
                violin = ax.violinplot(values, positions=[j], showextrema=False, widths=0.62)
                for body in violin['bodies']:
                    body.set(fc=color, ec='#333', lw=0.6, alpha=1)
            ax.scatter(np.full(len(values), j), values, s=2, c='#333')
        ax.set(title=m, ylim=(-1.05, 1.05), xticks=[0, 1], xticklabels=['High', 'Low'])
    for ax in list(axes.flat)[len(methods):]:
        ax.axis('off')
    save(fig, 'ED05b', 'Extended_Data_Fig06.xlsx / Paired cell-level values', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
