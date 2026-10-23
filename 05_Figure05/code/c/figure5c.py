"""Figure 5c paired spatial distribution. saved-data redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    import json
    stats = read('Figure05.xlsx', '5c Box statistics')
    items = []
    for row in stats.to_dict('records'):
        q05, q25, med, q75, q95 = json.loads(row['q05_q25_median_q75_q95'])
        items.append(dict(label=row['method'], whislo=q05, q1=q25, med=med, q3=q75, whishi=q95, fliers=[]))
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    ax.bxp(items, vert=False, showfliers=False)
    ax.invert_yaxis()
    ax.set(xlim=(-1.05, 1.05), xlabel='Spatial direction cosine', title='Saved paired-context consistency; 5th-95th percentile whiskers')
    save(fig, 'Figure05c', 'Figure05.xlsx / 5c Box statistics', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
