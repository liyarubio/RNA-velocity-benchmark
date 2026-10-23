"""ED11 archived GDC score check. archived score diagnostic; separate from current field panels."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    n = 17
    sheet = 'GDC direct edges'
    frame = read(f'Extended_Data_Fig{n}.xlsx', sheet)
    name = 'method' if 'method' in frame else 'Method'
    bar(frame.set_index(name).GDC, f'ED{n}_GDC', 'ST_raw (1,292 spots)' if n == 15 else 'KPTracer (754 cells; denominator 706)', f'Extended_Data_Fig{n}.xlsx / ' + sheet, 'GDC', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
