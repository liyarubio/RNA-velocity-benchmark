"""ED07c reference cosine. median-only diagnostic; formal boxplots retained as artwork."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    n = 13
    name = 'scNT'
    frame = read(f'Extended_Data_Fig{n}.xlsx', 'Method summary').set_index('method')
    bar(frame['median'], f'ED{n}_median', name + ': median per-cell reference cosine', f'Extended_Data_Fig{n}.xlsx / Method summary', 'Median cosine; saved pipelines', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
