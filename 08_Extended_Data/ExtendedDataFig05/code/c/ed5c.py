"""Extended Data 5c dimension sensitivity. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    d = read('Extended_Data_Fig06.xlsx', 'Dimension accuracy')
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for top, g in d.groupby('topology', sort=False):
        g = g.sort_values('dimension')
        kw = {'color': 'black', 'linewidth': 2} if top == 'total' else {'linewidth': 0.8}
        ax.plot(g.dimension, g.median_velocity_cosine, 'o-', ms=3, label=top, **kw)
    ax.set(xlabel='PCA dimension', ylabel='Median projected velocity cosine', title='Projected velocity accuracy; not score-vector agreement')
    ax.legend(fontsize=6, bbox_to_anchor=(1, 1))
    save(fig, 'ED05c', 'Extended_Data_Fig06.xlsx / Dimension accuracy', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
