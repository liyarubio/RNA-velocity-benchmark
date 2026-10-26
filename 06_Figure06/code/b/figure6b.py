"""Figure 6b read-depth consistency. saved-data diagnostic; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    reads = read('Figure06.xlsx', '6b Plotted vs full')
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for m, g in reads.groupby('Method'):
        g = g.sort_values('Ratio_i')
        ax.plot(g.Ratio_i, g.CosineSim, marker='.', lw=0.8, label=m)
    ax.set(xlabel='Read-depth fraction', ylabel='Mean velocity cosine similarity\nto full-depth prediction', title='Figure 6b: consistency; TopicVelo uses a transition proxy')
    ax.legend(fontsize=5, ncol=3)
    save(fig, 'Figure06b_stability', 'Figure06.xlsx / 6b Plotted vs full', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
