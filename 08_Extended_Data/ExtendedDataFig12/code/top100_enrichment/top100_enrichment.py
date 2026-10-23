"""Extended Data 12 top100 enrichment. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    d = read('Extended_Data_Fig18.xlsx', 'MS100 Heatmap').set_index('Description')
    values = -d.filter(like='_LogP_').apply(pd.to_numeric)
    values.columns = values.columns.str.replace('_LogP_', '', regex=False)
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    cmap = plt.get_cmap('YlOrRd').copy()
    cmap.set_bad('#dddddd')
    im = ax.imshow(values.mask(values == 0), aspect='auto', cmap=cmap, vmin=0)
    ax.set_yticks(range(len(values)), values.index)
    ax.set_xticks(range(len(values.columns)), values.columns, rotation=65, ha='right')
    ax.set_title('Top100: supplied representative enrichment clusters')
    fig.colorbar(im, ax=ax, label='-log10(P)')
    save(fig, 'ED12_top100_enrichment', 'Extended_Data_Fig18.xlsx / MS100 Heatmap; original table order', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
