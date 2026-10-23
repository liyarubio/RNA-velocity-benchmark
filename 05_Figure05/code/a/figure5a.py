"""Figure 5a ST_small GDC scores. score-only diagnostic; original fields and schematic retained as artwork."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    f5 = read('Figure05.xlsx', '5a GDC').set_index('Method')
    bar(f5.GDC, 'Figure05a_GDC', 'ST_small: 608 spots; current native TopicVelo', 'Figure05.xlsx / 5a GDC', 'GDC', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
