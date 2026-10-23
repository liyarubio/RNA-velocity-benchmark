"""Extended Data 5 pooled pairing diagnostic. pooled diagnostic only; formal a-c are separate."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    pairs = read('Extended_Data_Fig06.xlsx', 'Paired cell-level values')
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    ax.scatter(pairs.high, pairs.low, s=4, alpha=0.35)
    ax.plot([-1, 1], [-1, 1], c='grey', lw=0.7)
    ax.set(xlabel='High-dimensional cosine', ylabel='PCA2 cosine', title='ED6: 698 paired method–dataset observations')
    save(fig, 'ED06_paired', 'Extended_Data_Fig06.xlsx / Paired cell-level values', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
