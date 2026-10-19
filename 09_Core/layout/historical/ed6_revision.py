"""Restore original ED6c accuracy summaries, using the current saved method coverage."""
import csv

OLD_EN = ('c, Original-style topology-colored dimension-sensitivity curves: cosine similarity between the vectors of high-dimensional and dimension-specific performance scores, paired by method and dataset. This is agreement between performance profiles, not mean RNA-velocity accuracy. The total curve pools all finite pairs; pair counts and derived values are supplied. Seventeen available configurations are shown; missing pairs are excluded without imputation.')
NEW_EN = ('c, Cosine similarity between predicted and simulated true RNA velocity in PCA spaces of 2–30 dimensions, restoring the original accuracy-curve definition. Each colored point is the median of finite method–dataset cell-level cosine scores within one topology at that dimension; the black total curve is the median across all 42 datasets and available methods. These are projected velocity-accuracy scores, distinct from the high/low-dimensional score-vector consistency in Fig. 3d. Seventeen available configurations are included; missing scores are excluded without imputation. Source scores and valid counts are supplied.')
OLD_ZH = ('c, Dimension-sensitivity curves colored by topology pair high-dimensional and selected-dimensional performance scores for the same method and dataset, then calculate cosine similarity between the score vectors; total combines all finite pairs. This measures agreement between performance distributions, not mean RNA velocity accuracy. There are 17 configurations. Missing pairs are not imputed and source values are unchanged; pairing counts and derived curves are supplied.')
NEW_ZH = ('c, Restored accuracy curves compare predicted and simulated ground-truth velocities in 2-30 PCA dimensions. Each colored point is the median finite cell-level cosine score across methods and datasets for that topology and dimension. The black total curve is the median across all 42 datasets and available method scores. This measures projected velocity accuracy, distinct from agreement between high- and low-dimensional performance-score vectors in main Fig. 3d. The 17 available configurations are used; missing scores are excluded without imputation, and raw scores and valid counts are supplied.')

def revise_text(value):
    if not isinstance(value,str):return value
    return value.replace(OLD_EN,NEW_EN).replace(OLD_ZH,NEW_ZH).replace(
        'High- and low-dimensional velocity consistency','High- and low-dimensional velocity evaluation')

def update_workbook(workbook,base):
    # Replace the superseded derived statistic; raw score sheets remain unchanged.
    if 'Dimension consistency' in workbook:del workbook['Dimension consistency']
    if 'Dimension accuracy' in workbook:del workbook['Dimension accuracy']
    with (base/'Source_data/derived/ExtFig06c_dimension_accuracy.csv').open() as handle:
        reader=csv.DictReader(handle)
        ws=workbook.create_sheet('Dimension accuracy');ws.append(reader.fieldnames)
        for r in reader:
            ws.append([r['topology'],int(r['dimension']),float(r['median_velocity_cosine']),
                       int(r['valid_method_dataset_scores']),int(r['datasets']),int(r['methods'])])
