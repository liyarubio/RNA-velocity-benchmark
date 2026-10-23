"""Figure 5b selected paired spatial maps. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    data = read('Figure05.xlsx', '5b Paired directions')
    order = read('Figure05.xlsx', '5b Selected methods').method.tolist()
    cells = pd.read_csv(source_path('data/spatial_completion/spatial_cells.csv'))
    fig, axes = plt.subplots(1, len(order), figsize=(15, 3.5), constrained_layout=True)
    for ax, m in zip(axes, order):
        q = data[data.method == m]
        other = cells[~cells.cell_id.isin(q.cell_id)]
        ax.scatter(other.x, other.y, s=2, c='#d6d6d6', edgecolors='none')
        missing = q[q.spatial_cosine.isna()]
        ax.scatter(missing.x, missing.y, s=8, c='black', marker='x', linewidths=0.5)
        im = ax.scatter(q.x, q.y, s=4, c=q.spatial_cosine, cmap='RdBu_r', vmin=-1, vmax=1, edgecolors='none')
        ax.set_aspect('equal')
        ax.set_title(m)
        ax.axis('off')
    fig.colorbar(im, ax=list(axes), label='Paired spatial direction cosine', shrink=0.6)
    save(fig, 'Figure05b', 'Figure05.xlsx / 5b Paired directions; 5b Selected methods', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
