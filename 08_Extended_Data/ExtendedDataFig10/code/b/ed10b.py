"""Extended Data 10b paired spatial maps. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    p = read('Extended_Data_Fig16.xlsx', 'Paired spots')
    s = read('Extended_Data_Fig16.xlsx', 'Summary')
    order = s.sort_values('Spatial cosine (median)', ascending=False).Method.tolist()
    cells = pd.read_csv(source_path('data/spatial_completion/spatial_cells.csv'))
    other = cells[~cells.cell_id.isin(p.cell_id)]
    fig, axes = plt.subplots(4, 5, figsize=(14, 10), constrained_layout=True)
    for ax, m in zip(axes.flat, order):
        q = p[p.method == m]
        ax.scatter(other.x, other.y, s=2, c='#d6d6d6', edgecolors='none')
        bad = q[q.spatial_cosine.isna()]
        ax.scatter(bad.x, bad.y, s=8, c='black', marker='x', linewidths=0.5)
        im = ax.scatter(q.x, q.y, s=4, c=q.spatial_cosine, vmin=-1, vmax=1, cmap='RdBu_r', edgecolors='none')
        ax.set_aspect('equal')
        ax.set_title(m + f' ({q.spatial_cosine.median():.3f})')
        ax.axis('off')
    fig.colorbar(im, ax=list(axes.flat), label='Spatial cosine', shrink=0.6)
    save(fig, 'ED10b', 'Extended_Data_Fig16.xlsx / Paired spots; Summary', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
