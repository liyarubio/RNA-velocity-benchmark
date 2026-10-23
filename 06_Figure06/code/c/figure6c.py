"""Figure 6c frozen composite. saved-data diagnostic; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    f6 = read('Figure06.xlsx', '6c Summary').set_index('Method')
    bar(f6.Overall_Score, 'Figure06c_frozen', 'Frozen composite; stability excluded', 'Figure06.xlsx / 6c Summary', '0.6 accuracy rank + 0.4 scalability rank', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
