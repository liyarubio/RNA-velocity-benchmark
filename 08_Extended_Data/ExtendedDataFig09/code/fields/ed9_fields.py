"""Extended Data 9 spatial fields. saved-vector redraw; not the formal layout."""
from pathlib import Path
import sys
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '12_Provenance/current_sources.json').is_file())
sys.path.insert(0, str(ROOT / '09_Core/visualization'))
from panel_common import np, pd, plt, source_path, read, matrix, save, heat, bar, run_panel

def render(out, log):
    from field_common import field
    cells = read('Extended_Data_Fig15.xlsx', 'Cells and coordinates')
    proj = read('Extended_Data_Fig15.xlsx', 'Projected directions')
    topo = pd.read_csv(source_path('data/topovelo_historical_20261021/derived/ST_raw_directions.csv'))
    topo['method'] = 'TopoVelo'
    topic = read('Extended_Data_Fig15.xlsx', 'TopicVelo ST_raw')
    topic['method'] = 'TopicVelo'
    proj = pd.concat([proj[~proj.method.isin(['TopoVelo', 'TopicVelo'])], topo, topic], ignore_index=True)
    order = read('Extended_Data_Fig15.xlsx', 'GDC display order').method.tolist()
    target = {'Inner EPC 1': '#A8DADC', 'Inner EPC 2': '#FFD6A5', 'Outer EPC 1': '#B7E4C7', 'Outer EPC 2': '#FFADAD', 'P-TGC (Sec. P-TGC 1)': '#CDB4DB', 'P-TGC (Sec. P-TGC 2)': '#DDBEA9'}
    palette = {k: target.get(k, '#d6d6d6') for k in cells.cell_type.unique()}
    fig, axes = plt.subplots(4, 5, figsize=(14, 11), constrained_layout=True)
    for ax, m in zip(axes.flat, order):
        q = proj[proj.method == m].set_index('cell_id').reindex(cells.cell_id)
        assert len(q) == 1292 and q[['x', 'y']].notna().all().all(), m
        field(ax, cells[['x', 'y']].to_numpy(), q[['vx', 'vy']].to_numpy(), cells.cell_type.to_numpy(), palette, m)
        ax.set_aspect('equal')
    save(fig, 'ED09_fields', 'Extended_Data_Fig15.xlsx / Projected directions; TopicVelo ST_raw; GDC display order; archived TopoVelo ST_raw_directions.csv', out=out, log=log)

if __name__ == "__main__":
    run_panel(render)
