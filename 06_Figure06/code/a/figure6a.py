"""Figure 6a runtime. saved-data diagnostic; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    rt = read('Figure06.xlsx', '6a Runtime seconds').set_index('Cell')
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for c in rt:
        ax.plot(rt.index, rt[c], label=c, lw=0.8)
    ax.set(xlabel='Cells', ylabel='Seconds', yscale='log', title='Figure 6a: recorded runtime')
    ax.legend(fontsize=5, ncol=3)
    save(fig, 'Figure06a_runtime', 'Figure06.xlsx / 6a Runtime seconds', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
