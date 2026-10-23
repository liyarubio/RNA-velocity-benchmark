"""Extended Data 10 spatial medians. median-only diagnostic; formal a-b are separate."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    paired = pd.read_csv(source_path('data/spatial_context_consistency/paired_spots.csv'))
    bar(paired.groupby('method').spatial_cosine.median(), 'ED16_spatial', 'Same 608 target spots in two saved workflows', 'data/spatial_context_consistency/paired_spots.csv', 'Median spatial direction cosine', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
