"""Figure 6c weighted composite. saved-data diagnostic; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    f6 = read('Figure06.xlsx', '6c Summary').set_index('Method')
    bar(f6.Overall_Score, 'Figure06c_weighted', 'Weighted composite: 0.50 / 0.25 / 0.25', 'Figure06.xlsx / 6c Summary', '0.50 accuracy + 0.25 scalability + 0.25 stability rank index', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
