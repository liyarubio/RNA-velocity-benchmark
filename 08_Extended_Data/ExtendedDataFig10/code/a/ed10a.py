"""Extended Data 10a cell-type medians. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    d = read('Extended_Data_Fig16.xlsx', 'By cell type')
    s = read('Extended_Data_Fig16.xlsx', 'Summary')
    order = s.sort_values('Spatial cosine (median)', ascending=False).Method.tolist()
    grid = d[d.metric == 'spatial'].pivot(index='cell_type', columns='method', values='median').reindex(columns=order)
    heat(grid, 'ED10a', 'Spatial direction cosine by cell type', 'Extended_Data_Fig16.xlsx / By cell type', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
