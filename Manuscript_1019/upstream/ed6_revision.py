"""Restore original ED6c accuracy summaries, using the current saved method coverage."""
import csv

OLD_EN = ('c, Original-style topology-colored dimension-sensitivity curves: cosine similarity between the vectors of high-dimensional and dimension-specific performance scores, paired by method and dataset. This is agreement between performance profiles, not mean RNA-velocity accuracy. The total curve pools all finite pairs; pair counts and derived values are supplied. Seventeen available configurations are shown; missing pairs are excluded without imputation.')
NEW_EN = ('c, Cosine similarity between predicted and simulated true RNA velocity in PCA spaces of 2–30 dimensions, restoring the original accuracy-curve definition. Each colored point is the median of finite method–dataset cell-level cosine scores within one topology at that dimension; the black total curve is the median across all 42 datasets and available methods. These are projected velocity-accuracy scores, distinct from the high/low-dimensional score-vector consistency in Fig. 3d. Seventeen available configurations are included; missing scores are excluded without imputation. Source scores and valid counts are supplied.')
OLD_ZH = ('c，按拓扑着色的维度敏感性曲线：将同一方法、同一数据集的高维与指定维度性能分数配对，计算两个分数向量的余弦相似度；total 汇总全部有限配对。这是性能分布的一致性，不是平均 RNA 速度准确率。当前有 17 种配置，缺失配对不插补，来源表值不变；配对数量及派生曲线另附。')
NEW_ZH = ('c，恢复原版准确性曲线：在 2–30 维 PCA 空间比较预测速度与模拟真实速度。每个彩色点为该拓扑、该维度内有限的方法—数据集细胞级余弦得分的中位数；黑色 total 曲线为全部 42 个数据集及可用方法得分的中位数。该指标为投影后的速度准确性，与主图 3d 的高低维性能分数向量一致性不同。采用当前 17 种可用配置，缺失得分排除且不插补；附原始得分及有效数量。')

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
