"""ED08c reference and method fields. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    from field_common import field
    cells = pd.read_csv(source_path('data/reference_restoration/scLT_reference_cells.csv'))
    proj = pd.read_csv(source_path('data/reference_restoration/scLT_projected_directions.csv'))
    cos = pd.read_csv(source_path('data/reference_restoration/scLT_per_cell_cosine.csv')).set_index('cell_id')
    labels = cells.cell_type.to_numpy()
    sys.path.insert(0, str(ROOT / '09_Core/visualization/ed08_celltype'))
    from ed08_celltype import PALETTE
    palette = PALETTE
    order = cos.median().sort_values(ascending=False, kind='stable').index.tolist()
    fig, axes = plt.subplots(4, 5, figsize=(14, 10), constrained_layout=True)
    for ax, m in zip(axes.flat, ['Reference'] + order):
        if m == 'Reference':
            v = cells[['vx', 'vy']].to_numpy()
            title = 'Reference'
        else:
            q = proj[proj.method == m].set_index('cell_id').reindex(cells.cell_id)
            v = q[['vx', 'vy']].to_numpy()
            title = m + f' ({cos[m].median():.3f})'
        field(ax, cells[['x', 'y']].to_numpy(), v, labels, palette, title)
    save(fig, 'ED08c', 'data/reference_restoration/scLT_reference_cells.csv; projected_directions.csv; per_cell_cosine.csv', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
