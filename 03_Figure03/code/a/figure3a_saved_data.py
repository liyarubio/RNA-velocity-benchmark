"""Figure 3a high-dimensional cosine check. saved-data diagnostic; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    heat(matrix('Figure03.xlsx', '3a High dim'), 'Figure03a', 'Simulation: high-dimensional cosine', 'Figure03.xlsx / 3a High dim', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
