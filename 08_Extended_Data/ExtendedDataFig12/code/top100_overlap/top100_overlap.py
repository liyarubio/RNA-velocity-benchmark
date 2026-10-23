"""Extended Data 12 top100 overlap. intersection-matrix diagnostic; original chord artwork retained."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    n = 100
    members = pd.read_csv(source_path(f'data/metascape/top{n}/mapped_membership.csv'))
    counts = pd.read_csv(source_path(f'data/metascape/top{n}/list_counts.csv'))
    binary = members[counts.method].astype(int)
    overlap = binary.T @ binary
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    im = ax.imshow(overlap, cmap='Blues')
    ax.set_xticks(range(9))
    ax.set_xticklabels(overlap.columns, rotation=65, ha='right')
    ax.set_yticks(range(9))
    ax.set_yticklabels(overlap.index)
    for i in range(9):
        for j in range(9):
            ax.text(j, i, str(overlap.iloc[i, j]), ha='center', va='center', fontsize=6, color='white' if overlap.iloc[i, j] > n * 0.6 else 'black')
    ax.set_title(f'ED18: top{n} mapped-gene intersections')
    fig.colorbar(im, ax=ax, label='Shared mapped genes')
    save(fig, f'ED18_top{n}_overlap', f'data/metascape/top{n}/mapped_membership.csv', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
