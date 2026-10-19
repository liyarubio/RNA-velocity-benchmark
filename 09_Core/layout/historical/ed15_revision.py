"""Current ED15 reference schematic and caption; frozen inputs remain historical."""
CAPTION_EN = (
    'a, Two views of all 1,292 supplied ST_raw spots in identical fixed spatial coordinates: all 12 cell types in their original colors, and the six target cell types (608 spots) highlighted in the Fig. 5a colors with the other 684 spots in grey. '
    'The reference schematic is reused directly from Fig. 5a. Its five arrows indicate predefined evaluation relations, not measured spatial trajectories. '
    'b, GDC values displayed in descending order. Saved corrected eval2 GDC is supplemented with DeepVelo(VAE) re-evaluation from its saved velocity and displayed spatial projection (GDC = -0.493). '
    'The corrected eval2 weighting uses the 608 spots in the six reference types as denominator; three archived scores were reproduced within 1e-6 to verify this protocol. '
    'TFvelo retains its separately recorded evaluation with its original denominator and preprocessing; the bars are not a new uniform-protocol ranking. '
    'c, Nineteen available method fields arranged left to right and top to bottom in the same descending-GDC order as panel b, with identical spot coordinates and a shared grey background for non-target types. '
    'TopicVelo uses native stochastic transitions rerun from the saved topic fits after excluding genes without variable U/S signal and topics with fewer than two usable genes. Eleven topics are integrated for ST_raw and five for ST_small, replacing the historical last-topic-only configuration. Its high-dimensional direction is the transition-displacement proxy T @ Ms - Ms, and its spatial field is projected directly from the native transition matrix. GDC is -0.419 for ST_raw and 0.484 for ST_small. One ST_raw spot has an undefined direction and is excluded from directional scoring; the reference denominator remains 608 spots. '
    'Historical base-scVelo output is not used as TopicVelo. Existing directions for other methods are retained. All panels are combined on one page.'
)
CAPTION_ZH = (
    'a, Two spatial views of all 1,292 ST_raw spots use identical coordinates. One retains the original colors of all 12 cell types; the other highlights the 608 spots in six target types using main Fig. 5a colors, with the remaining 684 spots in grey. The reference schematic is reused from main Fig. 5a; its five arrows denote predefined evaluation relationships, not observed spatial trajectories. '
    'b, GDC scores are ordered from highest to lowest. Corrected eval2 GDC is retained, with an additional score for saved DeepVelo(VAE) velocities and their current spatial projection (-0.493). eval2 uses the 608 spots in six reference types as the weighting denominator. Three existing scores were reproduced by the same workflow with error below 1e-6. '
    'TFvelo retains its separate evaluation, original denominator and preprocessing; these bars are not presented as a new ranking under a uniform protocol. '
    'c, Fields for 19 valid methods follow the descending GDC order in b, read left to right and top to bottom, with common coordinates and non-target types in grey. TopicVelo native stochastic transitions were re-inferred from saved topic fits after excluding genes without variable U/S signal and topics with fewer than two valid genes. ST_raw integrates 11 topics and ST_small five topics, replacing the historical last-topic-only configuration. High-dimensional directions use the native transition-displacement proxy T @ Ms - Ms, and spatial fields are projected directly from native transitions; GDC is -0.419 and 0.484, respectively. One ST_raw spot has an undefined direction and is excluded from direction evaluation; the reference denominator remains 608 spots. Other method directions are unchanged. Panels a-c occupy one page.'
)
CAPTION = CAPTION_EN + '\n' + CAPTION_ZH

def update_workbook(w, base):
    import csv, json
    if 'Previous GDC source values' not in w.sheetnames:
        w['GDC source values'].title = 'Previous GDC source values'
    for name in ['GDC source values','Spatial revision 1010','GDC display order','Full cell-type palette']:
        if name in w.sheetnames:del w[name]
    s = w.create_sheet('GDC source values')
    for row in csv.reader((base/'Source_data/spatial_completion/current_scores.csv').open()):
        if row[0]!='method':row[1]=float(row[1]) if row[1] else None
        s.append(row)
    s = w.create_sheet('Spatial revision 1010')
    s.append(['field','value'])
    s.append(['Target spots',608]);s.append(['Other spots shown grey',684])
    s.append(['DeepVelo(VAE) GDC',-0.49306762019139627])
    s.append(['TopicVelo ST_raw and ST_small','Native multi-topic transitions complete after U/S signal QC; see topicvelo_retry/summary.json'])
    s.append(['Reference schematic','Exact Fig5a reference artwork'])
    s.append(['Panel a views','Original full cell-type colors; target types highlighted'])
    s.append(['Panel b order','Descending GDC; all 19 configurations available'])
    for title,name in [('GDC display order','GDC_display_order.csv'),('Full cell-type palette','all_cell_type_palette.csv')]:
        sheet=w.create_sheet(title)
        for row in csv.reader((base/'Source_data/spatial_completion'/name).open()):sheet.append(row)
    for scope in ['ST_raw','ST_small']:
        title='TopicVelo '+scope
        if title in w.sheetnames:del w[title]
        sheet=w.create_sheet(title)
        for row in csv.reader((base/'Source_data/spatial_completion/topicvelo_retry'/(scope+'_directions.csv')).open()):sheet.append(row)

def revise_document(doc):
    from lxml import etree
    section = doc.get_element_by_id('ed15')
    figures = section.xpath('./figure')
    for figure in figures[1:]:section.remove(figure)
    for node in figures[0].iter():
        for attr in ('src','href'):
            if node.get(attr):node.set(attr,node.get(attr).replace('Fig15_overview','Fig15'))
    figures[0].find('h3').text='Extended Data Fig. 15 | Spatial transcriptomics: reference, scores and method fields'
    cap = figures[0].find('figcaption')
    for child in list(cap):
        cap.remove(child)
    cap.text = None
    for lang, cls, content in [('en', 'en', CAPTION_EN), ('zh-CN', 'zh', CAPTION_ZH)]:
        etree.SubElement(cap, 'p', lang=lang, attrib={'class': cls}).text = content
    revise_main_caption(doc)
    maincaps = doc.xpath('//*[@id="ppt-figure-5"]//figcaption')
    for n in (maincaps[0].iter() if maincaps else []):
        if n.text:
            n.text=n.text.replace('with 19 method configurations.', 'with 18 available method fields and one failed-result slot for TopicVelo. Cell-marker radius is increased 1.6-fold for visibility. Native TopicVelo transition inference was rerun from the saved ST_small topic fit but failed (SVD did not converge); the historical base-scVelo field and score are removed.')
    revise_main_caption(doc)
    replacements = {
        'native TopicVelo outputs are missing for both settings': 'spatial TopicVelo is available from documented multi-topic re-inference, while native TopicVelo remains unavailable for KPTracer',
        'Both scenarios lack native TopicVelo outputs': 'Spatial TopicVelo outputs have been completed through documented multi-topic re-inference; native KPTracer TopicVelo outputs remain missing',
        '608-cell ST_small subset': '608-spot ST_small subset',
        'all 1,292 cells in the supplied ST_raw input': 'all 1,292 spots in the supplied ST_raw input',
        'the ST_small subset of 608 cells': 'the ST_small subset of 608 spots',
        'all 1,292 cells in the ST_raw input': 'all 1,292 spots in the ST_raw input',
        'Centroid arrows illustrate the relations used for evaluation, not measured trajectories of individual cells.':
        'The spatial reference relations are shown in a separate labeled schematic; the KPTracer centroid arrows illustrate its evaluation relations. Neither represents measured trajectories of individual cells.',
        'Arrows between centroids represent evaluation relationships, not measured trajectories of individual cells.':
        'Spatial reference relationships are shown in a separate labeled branching schematic; KPTracer centroid arrows represent its evaluation relationships. Neither depicts measured single-cell trajectories.',
    }
    for node in doc.iter():
        for attr in ('text', 'tail'):
            value = getattr(node, attr)
            if value:
                for old, new in replacements.items():
                    value = value.replace(old, new)
                setattr(node, attr, value)


def revise_main_caption(doc):
    old='with 18 available method fields and one failed-result slot for TopicVelo. Cell-marker radius is increased 1.6-fold for visibility. Native TopicVelo transition inference was rerun from the saved ST_small topic fit but failed (SVD did not converge); the historical base-scVelo field and score are removed.'
    new='with 19 available method fields. Cell-marker radius is increased 1.6-fold for visibility. TopicVelo uses five native stochastic topic transitions rerun from saved fits after U/S signal QC, replacing the historical last-topic-only configuration. Its GDC (0.484) uses T @ Ms - Ms as a transition-displacement proxy; the spatial field is projected from native transitions. This is not its historical base-scVelo output.'
    for n in doc.xpath('//*[@id="ppt-figure-5"]//figcaption')[0].iter() if doc.xpath('//*[@id="ppt-figure-5"]//figcaption') else []:
        if n.text:n.text=n.text.replace(old,new).replace('(608 cells)','(608 spots)').replace('TFvelo is a new run with distinct preprocessing; the other method panels are historical.','TFvelo has distinct preprocessing; TopicVelo uses the documented rerun above; the other method panels are historical.')
