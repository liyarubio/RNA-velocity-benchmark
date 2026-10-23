"""Extended Data 4b simulation cosine. saved-data diagnostic; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    sheet = 'Cell-level PCA2 cosine'
    label = 'b'
    heat(matrix('Extended_Data_Fig05.xlsx', sheet), 'ED04' + label, 'dyngen: ' + sheet, 'Extended_Data_Fig05.xlsx / ' + sheet, out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
