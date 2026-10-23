"""Extended Data 8a sampling times. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    d = pd.read_csv(source_path('data/reference_restoration/scLT_timepoints.csv'))
    fig, axes = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)
    for ax, samples in zip(axes, [['d2_5'], ['d5r1', 'd5r2']]):
        for sample in samples:
            q = d[d['sample'] == sample]
            ax.scatter(q.x, q.y, s=2, label=sample, alpha=0.8, edgecolors='none')
        ax.set(xlabel='CCE1', ylabel='CCE2')
        ax.legend()
    save(fig, 'ED08a', 'data/reference_restoration/scLT_timepoints.csv', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
