"""Extended Data 11a cell annotation. annotation-only redraw; literature schematic retained as artwork."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    d = read('Extended_Data_Fig17.xlsx', 'Cells and coordinates')
    p = read('Extended_Data_Fig17.xlsx', 'Cell-type palette')
    colors = dict(zip(p.cell_type, p.color))
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    for cell_type, g in d.groupby('cell_type', sort=False):
        ax.scatter(g.x, g.y, s=8, color=colors[cell_type], label=cell_type, edgecolors='none')
    ax.set(xlabel='UMAP1', ylabel='UMAP2', title='KPTracer: 754 cells; reference schematic retained separately')
    ax.set_aspect('equal')
    ax.legend(fontsize=6)
    save(fig, 'ED11a_annotation', 'Extended_Data_Fig17.xlsx / Cells and coordinates; Cell-type palette', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
