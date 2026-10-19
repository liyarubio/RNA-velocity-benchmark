"""Current reference-figure captions, workbook import and text corrections."""
import csv
COMMON_EN='Numbers below method names in the stream panels are the median per-cell reference-direction cosines, shown to three decimal places in the style of main Fig. 4a. Boxes show the median and interquartile range, with 1.5-IQR whiskers and outliers; methods are ordered by median. All 19 configurations use barcode-aligned cells and the same reference coordinates. Fourteen configurations come from the supplied historical H5AD results; SDEvelo uses the previous preprocessed-input run, and TopicVelo, biVI, DeepKINET and cell2fate use the aligned raw-count extension. These runs do not constitute a uniform-preprocessing ranking. TFvelo projections were restored from saved velo_hat / fit_scaling_y fits using the original graph-and-embedding workflow because the stored projection duplicated scVelo(sto); no model was retrained. Cosines were calculated from the displayed per-cell directions. Main Fig. 6 retains its previously selected run summaries and was not recomputed.'
COMMON_ZH='流线图的方法名下方标注逐细胞参考方向余弦的中位数，保留三位小数，形式参考主图 4a。箱体表示中位数及四分位距，须线为 1.5 IQR，显示离群点，按方法中位数排序。全部 19 个配置按细胞条形码对齐，采用相同参考坐标。其中 14 个配置来自新提供的历史 H5AD，SDEvelo 来自既有预处理运行，TopicVelo、biVI、DeepKINET、cell2fate 来自对齐原始计数扩展；并非统一预处理排名。TFvelo 原存储投影与 scVelo(sto) 重复，现按原流程从已保存的 velo_hat / fit_scaling_y 拟合重建速度图及投影，未重训模型。余弦由图示逐细胞方向计算。主图 6 保留此前选定运行的汇总，未随本次重算。'
CAPTIONS={
13:('a, UMAP annotation of 1,947 scNT-seq cells, restoring the original Extended Data Fig. 7 design and cell-type colors. b, Velocity streams for the dynamo metabolic-labeling reference and all 19 method configurations, arranged in five columns. c, Per-cell cosine similarity to the metabolic-labeling reference. '+COMMON_EN+'\n恢复原版扩展图 7。a，1947 个 scNT-seq 细胞的 UMAP 类型注释及原版配色。b，基于代谢标记的 dynamo 参考与全部 19 个方法配置的五列流线图。c，相对代谢标记参考的逐细胞方向余弦。'+COMMON_ZH),
14:('a, CCE views of Day 2.5 (1,756 cells) and Day 5 (675 and 755 cells in two replicates), restoring the original Extended Data Fig. 8 design. b, Original clone-transition examples for clone IDs 4, 6, 7, 23 and 32 retained as historical artwork; the exact clone-ID membership table is unavailable, so these arrows were not reconstructed from the separately saved expected-target pairs. c, Clone-based reference streams and all 19 method configurations for 1,756 cells, using the original cell-type grouping and colors. d, Per-cell cosine similarity to the clone-based reference. '+COMMON_EN+' TopicVelo has 1,724 valid cosines; its 32 zero-direction cells remain undefined and are excluded from its boxplot. All other configurations have 1,756 valid cosines.\n恢复原版扩展图 8。a，Day 2.5 的 1756 个细胞与 Day 5 两次重复的 675、755 个细胞的 CCE 分布。b，原版克隆 4、6、7、23、32 的转移示例，因缺少精确克隆编号成员表而保留历史图稿，未用另存的期望目标配对替代。c，1756 个参考细胞的克隆参考及全部 19 个配置，采用原版细胞类型分组与配色。d，相对克隆参考的逐细胞方向余弦。'+COMMON_ZH+'TopicVelo 有 1724 个有效余弦，32 个零方向细胞保持未定义并从箱线图排除；其他配置均为 1756 个有效余弦。')}
def legacy_sheet_name(name):return ('Previous '+name)[:31]
def revise_reference_document(document):
    replacements={
      'The current comparison includes 1,947 reference cells and seven available configurations: stochastic and dynamical scVelo, SDEvelo, TopicVelo, biVI, DeepKINET and cell2fate. scVelo and SDEvelo use previously preprocessed inputs, whereas the other four methods come from an aligned raw-count extension.':'The restored comparison includes 1,947 reference cells and all 19 configurations: 14 supplied historical method results, the previous SDEvelo run, and the aligned raw-count extensions for TopicVelo, biVI, DeepKINET and cell2fate. Method-specific preprocessing and run provenance are retained in the source data.',
      '当前图包含 1947 个参考细胞及七种可用配置：scVelo 的 stochastic 和 dynamical 模式、SDEvelo、TopicVelo、biVI、DeepKINET 与 cell2fate。其中 scVelo 和 SDEvelo 使用既有预处理输入，另外四种方法来自对齐原始计数的扩展分析。':'恢复后的图包含 1947 个参考细胞及全部 19 个配置：新提供的 14 个历史方法结果、既有 SDEvelo 运行，以及 TopicVelo、biVI、DeepKINET、cell2fate 的对齐原始计数扩展。各方法的预处理和运行来源分别记录在源数据中。',
      'per-cell cosine similarities for seven available configurations':'per-cell cosine similarities for all 19 configurations',
      '及七种可用配置与参考之间的逐细胞余弦相似度':'及全部 19 个配置与参考之间的逐细胞余弦相似度'}
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
        zh=zh.replace('七配置参考面板展示的是具有匹配当前投影证据的子集，两种覆盖范围不能混同。','恢复后的参考面板展示全部 19 个配置的可追溯投影，运行选择单独记录；本次恢复图稿不更新综合分数。')
    if id=='P122':
        en=en.replace('the underlying per-cell vectors and cosine values are unchanged.','19 configurations are aligned by barcode. TFvelo graph projections are restored from saved fits; other saved per-cell vectors are retained, and cosines are recomputed against the common reference without model retraining.')
        zh=zh.replace('底层逐细胞向量与余弦分数不变。','19 个配置按条形码对齐；TFvelo 从已保存拟合重建速度图投影，其他逐细胞向量保留，并对共同参考重新计算余弦，未重训模型。')
    if id=='P120':
        en=en.replace('the seven available current configurations replace whole method summaries','the seven configurations available when the composite was assembled replaced whole method summaries')
        en+=' The later 19-configuration restoration of the reference figures is documented separately and does not alter these frozen composite inputs.'
        zh=zh.replace('七种当前可用配置替换整项方法汇总','综合表组装时可用的七种配置替换了整项方法汇总')
        zh+=' 后续恢复的 19 配置参考图单独记录，未改变这些冻结的综合评分输入。'
    return en,zh
