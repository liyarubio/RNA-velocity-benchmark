"""Figure 5d KPTracer directions. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    from field_common import field
    v = read('Figure05.xlsx', '5d Plot vectors')
    cells = read('Extended_Data_Fig17.xlsx', 'Cells and coordinates').set_index('cell_id')
    p = read('Extended_Data_Fig17.xlsx', 'Cell-type palette')
    palette = dict(zip(p.cell_type, p.color))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for ax, m in zip(axes, ['PhyloVelo', 'UniTVelo', 'DeepVelo(VAE)']):
        q = v[v.method == m]
        labels = cells.loc[q.cell_id, 'cell_type'].to_numpy()
        field(ax, q[['x', 'y']].to_numpy(), q[['vx', 'vy']].to_numpy(), labels, palette, m + (' (backward)' if m == 'PhyloVelo' else ''))
        ax.set_aspect('equal')
    save(fig, 'Figure05d', 'Figure05.xlsx / 5d Plot vectors; native PhyloVelo backward vectors retained', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
