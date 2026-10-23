"""Figure 5e KPTracer scores. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    d = read('Figure05.xlsx', '5e KPTracer GDC').sort_values('display_order')
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    colors = ['#b57a4b' if m == 'PhyloVelo' else '#40769c' for m in d.method]
    ax.barh(d.method, d.GDC, color=colors)
    ax.invert_yaxis()
    ax.set(xlabel='Displayed GDC', title='PhyloVelo +1 is an adapted diagnostic, excluded from comparable ranking')
    save(fig, 'Figure05e', 'Figure05.xlsx / 5e KPTracer GDC; comparison_status retained', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
