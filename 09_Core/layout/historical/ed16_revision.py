"""Current ED16 UMAP axes, separate differentiation schematic and captions."""
CAPTION_EN = (
    'a, Left, all 754 supplied cells in the fixed UMAP, with UMAP1 and UMAP2 axes and colors from the source Cluster-Name annotation. '
    'Right, a separate differentiation schematic uses color-matched circular nodes, cell-state labels and black arrows, following the visual style of Fig. 4a. '
    'The four arrows denote the direct-evidence reference relations used by corrected eval1; schematic positions are not UMAP coordinates or measured cellular trajectories. '
    'No centroid arrows are overlaid on the UMAP. Lung Mixed in historical evaluation labels corresponds to Lung progenitor-like here. '
    'b, Saved corrected eval1 GDC. Neither eval1_raw nor likely-edge eval2 is mixed in. '
    'TopicVelo is omitted because native transition output is absent; TFvelo output was not found. '
    'These are descriptive results from this supplied dataset, not independent biological replication.'
)
CAPTION_ZH = (
    'a, Left: a fixed UMAP of all 754 cells, with UMAP1/UMAP2 labels and colors from the source Cluster-Name annotations. '
    'Right: a separate differentiation schematic following main Fig. 4a, with matching colored nodes, cell-state names and black arrows. '
    'Four arrows denote the direct-evidence reference relationships used in corrected eval1. Schematic positions are neither UMAP coordinates nor measured single-cell trajectories. '
    'Centroid arrows are no longer overlaid on the UMAP. Historical Lung Mixed corresponds to Lung progenitor-like here. '
    'b, Retain corrected eval1 GDC without mixing eval1_raw or eval2 based on inferred relationships. '
    'Native TopicVelo transitions and matching TFvelo results are unavailable and remain missing. Results describe this dataset only and do not constitute independent biological replication.'
)
from ed16_completion import ready
COMPLETED=ready()
FIELD_CAPTION=None
if COMPLETED:
    CAPTION_EN=CAPTION_EN.split('b, Saved corrected eval1 GDC.')[0]+(
        'b, GDC under the four direct-evidence relations, normalized by the 706 cells in the six participating states. '
        'The 17 saved corrected eval1 scores are retained; outlined bars add the KPTracer TFvelo refit and TopicVelo native-transition completion. '
        'TopicVelo is a transition-displacement proxy: within-state coherence uses T @ Ms - Ms, and boundary directions use the native T projected into the fixed UMAP. '
        'It does not use the historical base-scVelo velocity. Gene sets and preprocessing differ by method; scores do not establish absolute kinetic accuracy. '
        'The frozen main Figure 6 composite inputs are unchanged. These are descriptive results from one supplied dataset.'
    )
    CAPTION_ZH=CAPTION_ZH.split('b, Retain corrected eval1')[0]+(
        'b, Evaluation uses four direct-evidence relationships and a denominator of 706 cells across the six participating types. The 17 existing corrected eval1 scores are retained; outlined bars add the TFvelo refit and completed native TopicVelo transitions. '
        'TopicVelo uses a transition-displacement proxy: within-cluster coherence uses T @ Ms - Ms, while boundary directions use native T projected into the fixed UMAP, without historical baseline scVelo velocities. '
        'Gene coverage and preprocessing vary across methods; scores do not measure absolute kinetic accuracy. Frozen composite inputs for main Fig. 6 are unchanged. Results describe this single dataset only.'
    )
    FIELD_CAPTION=(
        'c, Nineteen available method fields on the same 754 cells. Seventeen verified historical vector fields are retained. '
        'TFvelo was refitted after uppercasing gene symbols following its bundled demo; finite fitted genes determine its velocity graph. '
        'TopicVelo resumes the saved fastTopics fit with topic 7 of 8 and the original 35th-percentile threshold (490 topic cells), producing a native integrated transition matrix over all 754 cells. '
        'Its field uses explicit native transitions, with 727 nonzero projected directions; the package fallback for cells outside the selected topic is retained. '
        'New fields are drawn from saved per-cell directions; their interpolation is recorded separately from the preserved historical streams. Deterministic dynamo remains omitted.'
    )
CAPTION_EN=(
    'a, All 754 KPTracer cells in the fixed UMAP, with UMAP1/UMAP2 axes and source cell-type colors. '
    'The schematic uses the plain node-and-arrow style of Fig. 15, without background regions, and is adapted from Yang et al. (Cell 2022, doi:10.1016/j.cell.2022.04.015), Fig. 5I. Solid arrows mark the four relations used for GDC; dashed connections show other literature relations and are excluded from GDC. The dashed horizontal connection denotes an undirected coupling. Line style encodes GDC inclusion, not the original evidence categories. Hollow nodes are literature states absent from this subset. Unconnected mesenchymal states are omitted for clarity. The source Gastric/Endoderm progression is expanded into the evaluated Endoderm relation and a contextual Gastric branch; Early EMT-1/2 is split to identify the evaluated EMT-1 state. Node positions are schematic, not UMAP coordinates or measured trajectories. '
    'Historical Lung Mixed labels correspond to Lung progenitor-like. '
    'b, GDC ranked in descending order, using the 706 cells in the six reference-participating states as denominator. '
    'The 17 historical corrected eval1 scores and the completed TFvelo/TopicVelo scores are unchanged. '
    'TopicVelo denotes a native-transition displacement proxy: within-state coherence uses T @ Ms - Ms; boundary directions use native T projected into the fixed UMAP, without historical base-scVelo velocity. '
    'c, Nineteen method fields in the same order as b, read left to right and then top to bottom. '
    'Saved streamline paths are aligned by per-panel affine transforms to identical cell-coordinate bounds, with uniformly enlarged cell markers. Panel a also uses enlarged cell markers. TFvelo uses 1,950 finite fitted genes. '
    'TopicVelo resumes topic 7 of 8 from the saved fit at the 35th-percentile threshold (490 topic cells; 727 nonzero projected directions), retaining the native package fallback. '
    'Gene sets and preprocessing differ by method. These are descriptive results from one dataset; frozen main Fig. 6 composite inputs remain unchanged.'
)
CAPTION_ZH=(
    'a, Fixed UMAP of 754 KPTracer cells, retaining UMAP1/UMAP2 labels and source cell-type colors. '
    'The schematic uses the background-free node-and-arrow style of Fig. 15 and is adapted from Fig. 5I of Yang et al. (Cell 2022, doi:10.1016/j.cell.2022.04.015). Solid arrows show the four GDC relationships; dashed arrows show other literature relationships excluded from GDC, and horizontal dashed lines denote undirected coupling. Line styles encode inclusion in GDC, not the original evidence grading. Open nodes indicate literature states absent from this subset; the unconnected mesenchymal state is omitted for clarity. Progression of the original Gastric/Endoderm group is expanded into evaluated Endoderm edges and a contextual Gastric branch. Early EMT-1/2 are separated to identify evaluated EMT-1. Node positions are neither UMAP coordinates nor measured trajectories. Historical Lung Mixed corresponds to Lung progenitor-like. '
    'b, GDC is ordered from highest to lowest, using 706 cells in six reference-participating types as the denominator. The 17 historical corrected eval1 scores and completed TFvelo/TopicVelo scores are retained. '
    'TopicVelo denotes a native transition-displacement proxy: within-cluster coherence uses T @ Ms - Ms, and boundary directions use native T projected into the fixed UMAP, without historical baseline scVelo velocities. '
    'c, Fields for 19 methods follow the order in b, read left to right and top to bottom. Existing streams are aligned to common cell-coordinate bounds by per-panel affine transforms, and cell markers are enlarged consistently, including in a. TFvelo uses 1,950 valid fitted genes. '
    'TopicVelo retains saved fitted topics 7/8 and the 35th-percentile threshold, with 490 topic cells, 727 nonzero projected directions and native fallback rules. '
    'Gene coverage and preprocessing differ across methods. Results describe this single dataset only; frozen composite inputs for main Fig. 6 are unchanged.'
)
FIELD_CAPTION=None
CAPTION=CAPTION_EN+'\n'+CAPTION_ZH

def revise_document(doc):
    from lxml import etree
    section=doc.get_element_by_id('ed16');figures=section.xpath('./figure')
    for figure in figures[1:]:section.remove(figure)
    for node in figures[0].iter():
        for attr in ('src','href'):
            if node.get(attr):node.set(attr,node.get(attr).replace('Fig16_overview','Fig16').replace('Fig16_method_fields','Fig16'))
    figures[0].find('h3').text='Extended Data Fig. 16 | KPTracer: reference, scores and method fields'
    cap=figures[0].find('figcaption')
    for child in list(cap):cap.remove(child)
    cap.text=None
    for lang,cls,content in [('en','en',CAPTION_EN),('zh-CN','zh',CAPTION_ZH)]:
        etree.SubElement(cap,'p',lang=lang,attrib={'class':cls}).text=content
    if FIELD_CAPTION:
        field=doc.get_element_by_id('ed16').xpath('.//figcaption')[1]
        for child in list(field):field.remove(child)
        field.text=None
        etree.SubElement(field,'p',lang='en',attrib={'class':'en'}).text=FIELD_CAPTION
    replacements={
        'The spatial reference relations are shown in a separate labeled schematic; the KPTracer centroid arrows illustrate its evaluation relations. Neither represents measured trajectories of individual cells.':
        'Spatial and KPTracer reference relations are shown in separate labeled schematics. The KPTracer cell distribution is shown in a fixed UMAP with UMAP1 and UMAP2 axes. Solid schematic arrows represent GDC evaluation relations; dashed connections add literature context. Neither represents measured trajectories of individual cells.',
        'Spatial reference relationships are shown in a separate labeled branching schematic; KPTracer centroid arrows represent its evaluation relationships. Neither depicts measured single-cell trajectories.':
        'Spatial and KPTracer reference relationships are presented in separate labeled schematics. KPTracer cells use a fixed UMAP with UMAP1/UMAP2 labels. Solid schematic arrows denote GDC evaluation relationships and dashed arrows provide literature context; neither represents measured single-cell trajectories.',
    }
    replacements.update({
        'Schematic arrows represent evaluation relations, not measured trajectories of individual cells.': 'Solid schematic arrows represent GDC evaluation relations; dashed connections add literature context. Neither represents measured trajectories of individual cells.',
        'Schematic arrows denote evaluation relationships, not measured single-cell trajectories.': 'Solid schematic arrows denote GDC evaluation relationships; dashed arrows add literature context. Neither represents measured single-cell trajectories.',
    })
    if COMPLETED:
        replacements.update({
            'Saved score and field provenance is retained; native TopicVelo outputs are missing for both settings, and TFvelo is missing for KPTracer. Neither gap is filled using an output from another method or configuration.':
            'Saved score and field provenance is retained. KPTracer TFvelo has now been refitted and native TopicVelo transitions completed; its transition-displacement GDC proxy is identified separately. Spatial TopicVelo results are documented separately. No base-scVelo output is substituted for native TopicVelo transitions.',
            'Provenance for existing scores and fields is retained. Both scenarios lack native TopicVelo outputs, and KPTracer also lacks TFvelo. Outputs from other methods or configurations are not substituted for missing results.':
            'Provenance for existing scores and fields is retained. KPTracer TFvelo has been refitted and native TopicVelo transitions completed; the transition-displacement GDC proxy is labeled separately. Spatial TopicVelo results are documented independently. Baseline scVelo outputs do not replace native transitions.',
            'while native TopicVelo remains unavailable for KPTracer, and TFvelo is missing for KPTracer. Neither gap is filled using an output from another method or configuration.':
            'and KPTracer now includes a TFvelo refit and completed native TopicVelo transitions. TopicVelo GDC uses the documented transition-displacement proxy; no historical base-scVelo output is substituted.',
            'KPTracer still lacks native TopicVelo outputs and TFvelo. Outputs from other methods or configurations are not substituted for missing results.':
            'Native KPTracer TopicVelo transitions and the TFvelo refit are now complete. TopicVelo GDC uses the documented transition-displacement proxy, not historical baseline scVelo outputs.',
        })
    for node in doc.iter():
        for attr in ('text','tail'):
            value=getattr(node,attr)
            if value:
                for old,new in replacements.items():value=value.replace(old,new)
                setattr(node,attr,value)
