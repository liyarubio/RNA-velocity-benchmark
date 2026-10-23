"""Extended Data 7a annotation. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    from field_common import field
    cells = pd.read_csv(source_path('data/reference_restoration/scNT_reference_cells.csv'))
    labels = cells.cell_type.to_numpy()
    palette = {k: plt.get_cmap('tab10')(i) for i, k in enumerate(pd.unique(labels))}
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    field(ax, cells[['x', 'y']].to_numpy(), None, labels, palette, 'scNT annotation', annotation=True, point_size=12)
    save(fig, 'ED07a', 'data/reference_restoration/scNT_reference_cells.csv', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
