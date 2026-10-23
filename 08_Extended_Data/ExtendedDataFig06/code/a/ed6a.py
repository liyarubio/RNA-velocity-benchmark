"""Extended Data 6a HuBMAP raw GDC. raw-order heatmap diagnostic; does not reproduce formal clustering or organ strip."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    heat(matrix('Extended_Data_Fig07.xlsx', 'HuBMAP raw GDC').T, 'ED07_HuBMAP', 'HuBMAP: 176 datasets; raw source order', 'Extended_Data_Fig07.xlsx / HuBMAP raw GDC', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
