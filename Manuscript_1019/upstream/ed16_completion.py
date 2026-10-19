"""Saved KPTracer rerun evidence and workbook updates; no inference on import."""
from pathlib import Path
import json
B=Path(__file__).resolve().parents[1]
D=B/'Source_data/kptracer_completion'
METHODS={'tfvelo':'TFvelo','topicvelo':'TopicVelo'}
HISTORY={'GDC direct edges':'Previous GDC direct edges','Missing data':'Previous missing data',
         'Panel provenance':'Previous panel provenance','README':'Previous README'}
def ready():return all((D/(k+'.json')).exists() for k in METHODS)
def results():return {name:json.loads((D/(key+'.json')).read_text()) for key,name in METHODS.items()}
def update_workbook(w):
    if not ready():return
    for name,legacy in HISTORY.items():
        if legacy not in w.sheetnames:
            copy=w.copy_worksheet(w[name]);copy.title=legacy
    data=results()
    for row in w['GDC direct edges']:
        if row[0].value in data:row[1].value=data[row[0].value]['GDC']
    missing=w['Missing data']
    if missing.max_row>1:missing.delete_rows(2,missing.max_row-1)
    provenance=w['Panel provenance'];provenance.delete_rows(1,provenance.max_row)
    provenance.append(['scope','current_artwork','evidence','status'])
    provenance.append(['17 historical fields','Extended_Data_Fig16.pdf','Previous panel provenance; original vector paths aligned to shared cell bounds; enlarged cell markers','preserved'])
    for key,method in METHODS.items():
        provenance.append([method,'Extended_Data_Fig16.pdf','Source_data/kptracer_completion/'+key+'.json','completed and validated'])
    provenance.append(['UMAP','Extended_Data_Fig16.pdf','754 original coordinates and colors','preserved'])
    provenance.append(['Schematic','Extended_Data_Fig16.pdf','Yang et al. Cell 2022 Fig. 5I; GDC inclusion controls line style','redrawn without background regions, following Fig. 15'])
    readme=w['README'];readme.delete_rows(1,readme.max_row);readme.append(['note'])
    for note in ['KPTracer TFvelo refit and native TopicVelo transition completion, 2026-10-10.',
                 'Current scores: GDC direct edges. Current reruns: Rerun summary, Rerun directions and Rerun protocol.',
                 'Previous sheets and Output validation retain historical evidence; current rerun validation is in Source_data/kptracer_completion.',
                 'The 17 historical scores and fields are preserved; TopicVelo score is a native-transition displacement proxy.',
                 'Source inputs use historical normalized continuous values; method gene sets differ.',
                 'Cell-level distributions are not independent donor/culture replicates.']:
        readme.append([note])
    for name in ['Rerun summary','Rerun directions','Rerun protocol','GDC display order']:
        if name in w.sheetnames:del w[name]
    s=w.create_sheet('Rerun summary');s.append(['method','GDC','cells','genes','nonzero_direction_cells','reference_cell_denominator','prediction_type','source','source_sha256'])
    for method,r in data.items():s.append([method]+[r[k] for k in ['GDC','cells','genes','nonzero_direction_cells','reference_cell_denominator','prediction_type','source','source_sha256']])
    s=w.create_sheet('Rerun directions');s.append(['method','cell_id','x','y','vx','vy','cell_type'])
    import csv
    if (D/'GDC_display_order.csv').exists():
        ordered=w.create_sheet('GDC display order')
        with (D/'GDC_display_order.csv').open() as f:
            for row in csv.reader(f):ordered.append(row)
    for key,method in METHODS.items():
        with (D/(key+'_directions.csv')).open() as f:
            for r in csv.DictReader(f):s.append([method,r['cell_id']]+[float(r[k]) for k in ['x','y','vx','vy']]+[r['cell_type']])
    s=w.create_sheet('Rerun protocol');s.append(['field','value'])
    for row in [
        ['Scope','ED16 missing methods only; existing 17 method values and the frozen main Figure6 ranking remain unchanged'],
        ['Reference denominator','706 cells in six states involved in four direct-evidence edges; all 754 cells are displayed'],
        ['TFvelo','New fit; uppercase symbols following bundled official demo; velo_hat / fit_scaling_y; finite fitted genes'],
        ['TopicVelo','Saved fastTopics fit resumed, topic 7 of 8, percentile 35, 490 selected topic cells; native integrated transition matrix'],
        ['TopicVelo GDC proxy','ICCoh uses T @ Ms - Ms; CBDir uses scVelo projection of explicit native T; historical base-scVelo velocity is not used'],
        ['Input boundary','Historical normalized continuous spliced/unspliced inputs; method gene sets differ; this is not a matched raw-count benchmark'],
        ['Control','Historical scVelo(sto) GDC reproduced within 1e-6 using saved neighbor indices and denominator 706']]:s.append(row)

    record=B/'Provenance/ed16_single_page.json'
    if record.exists():
        scheme=json.loads(record.read_text()).get('literature_schematic')
        if scheme:
            for name in ['Literature schematic','Literature connections']:
                if name in w.sheetnames:del w[name]
            sheet=w.create_sheet('Literature schematic');sheet.append(['field','value'])
            for key in ['source_figure','doi','url','hollow_nodes','line_style_semantics']:
                sheet.append([key,scheme[key]])
            sheet.append(['context_only_nodes',', '.join(scheme['context_only_nodes'])])
            sheet.append(['omitted_unconnected_context_nodes',', '.join(scheme.get('omitted_unconnected_context_nodes',[]))])
            sheet.append(['Grouping','Original Gastric/Endoderm progression expanded; Early EMT-1/2 split to identify the evaluated state.'])
            sheet=w.create_sheet('Literature connections')
            sheet.append(['source','target','used_for_GDC','line_style','directed','literature_meaning'])
            for row in scheme['edges']:
                sheet.append([row[k] for k in ['source','target','GDC','line_style','directed','literature_meaning']])
