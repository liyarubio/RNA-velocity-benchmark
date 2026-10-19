"""Current reference-figure captions, workbook import and text corrections."""
import csv
COMMON_EN='Numbers below method names in the stream panels are the median per-cell reference-direction cosines, shown to three decimal places in the style of main Fig. 4a. Boxes show the median and interquartile range, with 1.5-IQR whiskers and outliers; methods are ordered by median. All 19 configurations use barcode-aligned cells and the same reference coordinates. Fourteen configurations come from the supplied historical H5AD results; SDEvelo uses the previous preprocessed-input run, and TopicVelo, biVI, DeepKINET and cell2fate use the aligned raw-count extension. These runs do not constitute a uniform-preprocessing ranking. TFvelo projections were restored from saved velo_hat / fit_scaling_y fits using the original graph-and-embedding workflow because the stored projection duplicated scVelo(sto); no model was retrained. Cosines were calculated from the displayed per-cell directions. Main Fig. 6 retains its previously selected run summaries and was not recomputed.'
COMMON_ZH='Numbers below method names in stream panels are median per-cell reference-direction cosines, shown to three decimal places following main Fig. 4a. Boxes show medians and interquartile ranges, with 1.5-IQR whiskers and outliers; methods are ordered by median. All 19 configurations use barcode-aligned cells and common reference coordinates. Fourteen configurations come from supplied historical H5AD results; SDEvelo uses the previous preprocessed-input run, while TopicVelo, biVI, DeepKINET and cell2fate use the aligned raw-count extension. These runs are not a uniform-preprocessing ranking. The stored TFvelo projection duplicated scVelo(sto), so its velocity graph and projection were reconstructed from saved velo_hat / fit_scaling_y fits using the original workflow, without model retraining. Cosines use the displayed per-cell directions. Main Fig. 6 retains previously selected run summaries and was not recomputed. '
CAPTIONS={
13:('a, UMAP annotation of 1,947 scNT-seq cells, restoring the original Extended Data Fig. 7 design and cell-type colors. b, Velocity streams for the dynamo metabolic-labeling reference and all 19 method configurations, arranged in five columns. c, Per-cell cosine similarity to the metabolic-labeling reference. '+COMMON_EN+'\nRestored original Extended Data Fig. 7. a, Cell-type annotation and original colors for a UMAP of 1,947 scNT-seq cells. b, Five-column stream plots for the dynamo metabolic-labeling reference and all 19 method configurations. c, Per-cell direction cosine relative to the metabolic-labeling reference. '+COMMON_ZH),
14:('a, CCE views of Day 2.5 (1,756 cells) and Day 5 (675 and 755 cells in two replicates), restoring the original Extended Data Fig. 8 design. b, Original clone-transition examples for clone IDs 4, 6, 7, 23 and 32 retained as historical artwork; the exact clone-ID membership table is unavailable, so these arrows were not reconstructed from the separately saved expected-target pairs. c, Clone-based reference streams and all 19 method configurations for 1,756 cells, using the original cell-type grouping and colors. d, Per-cell cosine similarity to the clone-based reference. '+COMMON_EN+' TopicVelo has 1,724 valid cosines; its 32 zero-direction cells remain undefined and are excluded from its boxplot. All other configurations have 1,756 valid cosines.\nRestored original Extended Data Fig. 8. a, CCE distributions for 1,756 Day 2.5 cells and two Day 5 replicates with 675 and 755 cells. b, Historical transition examples for clones 4, 6, 7, 23 and 32 are retained because exact clone-ID membership is unavailable; separately saved expected-target pairs are not substituted. c, Clone-reference streams and all 19 configurations for 1,756 reference cells, with original cell-type groups and colors. d, Per-cell direction cosine relative to the clone reference. '+COMMON_ZH+'TopicVelo has 1,724 valid cosines; 32 zero-direction cells remain undefined and are excluded from the boxplot. All other configurations have 1,756 valid cosines.')}
def legacy_sheet_name(name):return ('Previous '+name)[:31]
def revise_reference_document(document):
    replacements={
      'The current comparison includes 1,947 reference cells and seven available configurations: stochastic and dynamical scVelo, SDEvelo, TopicVelo, biVI, DeepKINET and cell2fate. scVelo and SDEvelo use previously preprocessed inputs, whereas the other four methods come from an aligned raw-count extension.':'The restored comparison includes 1,947 reference cells and all 19 configurations: 14 supplied historical method results, the previous SDEvelo run, and the aligned raw-count extensions for TopicVelo, biVI, DeepKINET and cell2fate. Method-specific preprocessing and run provenance are retained in the source data.',
      'The current figure contains 1,947 reference cells and seven available configurations: scVelo stochastic and dynamical, SDEvelo, TopicVelo, biVI, DeepKINET and cell2fate. scVelo and SDEvelo use existing preprocessed inputs; the other four methods use the aligned raw-count extension.':'The restored figure contains 1,947 reference cells and all 19 configurations: 14 newly supplied historical method outputs, the existing SDEvelo run, and aligned raw-count extensions for TopicVelo, biVI, DeepKINET and cell2fate. Method-specific preprocessing and run provenance are recorded in the source data.',
      'per-cell cosine similarities for seven available configurations':'per-cell cosine similarities for all 19 configurations',
      'and per-cell cosine similarities between seven available configurations and the reference':'and per-cell cosine similarities between all 19 configurations and the reference'}
    for section in document.xpath('//*[@id="P044" or @id="P045"]'):
        for node in section.iter():
            for attr in ['text','tail']:
                value=getattr(node,attr)
                if value:
                    for old,new in replacements.items():value=value.replace(old,new)
                    setattr(node,attr,value)
def update_workbook(w,n,B):
    # Preserve every old sheet, including its original values and provenance.
    for s in w:s.title=legacy_sheet_name(s.title)
    ds='scNT' if n==13 else 'scLT';root=B/'Source_data/reference_restoration'
    files={'Reference cells':ds+'_reference_cells.csv','Projected directions':ds+'_projected_directions.csv','Per-cell cosine':ds+'_per_cell_cosine.csv','Method summary':ds+'_method_summary.csv'}
    if n==14:files.update({'Full timepoint coordinates':'scLT_timepoints.csv','Saved expected targets':'scLT_saved_expected_targets.csv'})
    for name,file in files.items():
        s=w.create_sheet(name)
        with (root/file).open() as f:
            rows=csv.reader(f);head=next(rows);s.append(head)
            for row in rows:
                values=[]
                for value in row:
                    try:value=float(value) if value else None
                    except ValueError:pass
                    values.append(value)
                s.append(values)
        s.freeze_panes='B2';s.auto_filter.ref=s.dimensions
    s=w.create_sheet('README');s.append(['Field','Value'])
    for k,v in [('Current figure',n),('Original figure',n-6),('Methods',19),('Caption',CAPTIONS[n]),('Provenance','Provenance/reference_restoration.json; Provenance/tfvelo_projection_restoration.json'),('Reproduction','Code/plots/restore_reference_figures.py'),('Previous sheets','Archived pre-update 7-method tables; current panels use the unprefixed sheets.'),('Undefined cosine','Zero/nonfinite direction is missing, not zero.'),('Composite','Main Fig. 6 keeps its previously selected input runs.')]:s.append([k,v])
    w._sheets.sort(key=lambda sheet:sheet.title.startswith('Previous '))
    w.active=0
def revise_paragraph(id,en,zh):
    if id=='P051':
        en=en.replace('the seven-configuration reference panels show the subset with matched current projection evidence. These two coverage scopes must not be conflated.','the restored reference panels now show all 19 configurations from traceable saved projections. Their run selection is documented separately, and this figure restoration does not update the composite scores.')
        zh=zh.replace('The seven-configuration reference panels show the subset with matching current-projection evidence; the two coverage scopes must not be conflated.','Restored reference panels show traceable projections for all 19 configurations, with run selection recorded separately; restoring the artwork does not update composite scores.')
    if id=='P122':
        en=en.replace('the underlying per-cell vectors and cosine values are unchanged.','19 configurations are aligned by barcode. TFvelo graph projections are restored from saved fits; other saved per-cell vectors are retained, and cosines are recomputed against the common reference without model retraining.')
        zh=zh.replace('Underlying per-cell vectors and cosine scores are unchanged.','The 19 configurations are barcode-aligned. TFvelo graph projections are reconstructed from saved fits; other per-cell vectors are retained and cosines are recomputed against the common reference, without model retraining.')
    if id=='P120':
        en=en.replace('the seven available current configurations replace whole method summaries','the seven configurations available when the composite was assembled replaced whole method summaries')
        en+=' The later 19-configuration restoration of the reference figures is documented separately and does not alter these frozen composite inputs.'
        zh=zh.replace('the seven currently available configurations replace the complete method summaries','the seven configurations available when the composite table was assembled replaced the complete method summaries')
        zh+=' The subsequently restored 19-configuration reference figures are documented separately and do not change these frozen composite-score inputs.'
    return en,zh
