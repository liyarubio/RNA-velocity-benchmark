"""Figure 4a GDC scores. score-only diagnostic; original streams and schematic retained as artwork."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    f4 = read('Figure04.xlsx', '4a GDC').set_index('Method')
    bar(f4.GDC, 'Figure04a_GDC', 'ReDeeM: saved GDC', 'Figure04.xlsx / 4a GDC', 'GDC', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
