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
    'a，左侧为全部 754 个细胞的固定 UMAP，标注 UMAP1 和 UMAP2，按来源 Cluster-Name 注释着色。'
    '右侧参照主图 Fig. 4a，以同色圆点、细胞状态名称和黑色箭头单独展示分化关系示意图。'
    '四条箭头表示 corrected eval1 使用的直接证据参考关系；示意图位置不是 UMAP 坐标，也不是单个细胞的实测轨迹。'
    'UMAP 上不再叠加质心箭头。历史评价中的 Lung Mixed 对应此处的 Lung progenitor-like。'
    'b，沿用 corrected eval1 的 GDC，未混入 eval1_raw 或使用推测关系的 eval2。'
    'TopicVelo 缺少原生转移输出、TFvelo 未找到匹配结果，均保留缺失。结果仅描述该数据集，不构成独立生物学重复。'
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
    CAPTION_ZH=CAPTION_ZH.split('b，沿用 corrected eval1')[0]+(
        'b，按四条直接证据关系评价，以涉及的六类细胞共 706 个细胞为分母。保留已有 17 项 corrected eval1 分数；描边柱新增 TFvelo 重拟合和 TopicVelo 原生转移补全结果。'
        'TopicVelo 为转移位移代理：群内一致性使用 T @ Ms - Ms，边界方向来自原生 T 在固定 UMAP 的投影，未使用历史基础 scVelo 速度。'
        '方法间基因范围和预处理不同，分数不代表绝对动力学准确性。主图 6 的冻结综合评分输入保持原样。结果仅描述此单一数据集。'
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
    'a，754 个 KPTracer 细胞的固定 UMAP，保留 UMAP1/UMAP2 轴标及来源细胞类型配色。'
    '示意图采用 Fig. 15 的无背景圆点和箭头样式，改绘自 Yang 等（Cell 2022，doi:10.1016/j.cell.2022.04.015）Fig. 5I。实线箭头为 GDC 使用的四条关系，虚线为未纳入 GDC 的其他文献关系；水平虚线表示无向耦合。线型表示 GDC 纳入情况，不沿用原文的证据等级。空心节点为本数据子集中未出现的文献状态；为简化展示，省略无连线的间充质状态。原文 Gastric/Endoderm 群体的进展关系展开为参与评价的 Endoderm 连线和仅作背景的 Gastric 分支；Early EMT-1/2 拆开展示以标明参与评价的 EMT-1。节点位置不表示 UMAP 坐标或实测轨迹。历史 Lung Mixed 对应 Lung progenitor-like。'
    'b，GDC 按降序排列，以六类参考关系参与细胞共 706 个细胞为分母；保留 17 项历史 corrected eval1 分数及已补跑的 TFvelo、TopicVelo 分数。'
    'TopicVelo 表示原生转移位移代理：群内一致性使用 T @ Ms - Ms，边界方向来自原生 T 的固定 UMAP 投影，不使用历史基础 scVelo 速度。'
    'c，19 个方法流场与 b 同序，按从左到右、从上到下阅读；原有流线按每图的仿射变换统一到相同细胞坐标范围，细胞点统一增大；a 的细胞点也已增大。TFvelo 使用 1950 个有效拟合基因。'
    'TopicVelo 沿用已保存拟合的第 7/8 号主题和 35 百分位阈值，包含 490 个主题细胞、727 个非零投影方向，并保留原生回退规则。'
    '方法间基因范围和预处理不同，结果仅描述此单一数据集；主图 6 的冻结综合评分输入不变。'
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
        '空间数据的参考关系以独立且带名称的分支示意图展示；KPTracer 的质心箭头表示其评价关系。两者均不表示单个细胞的实测轨迹。':
        '空间数据与 KPTracer 的参考关系均以独立且带名称的示意图展示。KPTracer 的细胞分布采用固定 UMAP，并标注 UMAP1 和 UMAP2。示意图实线箭头表示 GDC 评价关系，虚线补充文献关系，均不表示单个细胞的实测轨迹。',
    }
    replacements.update({
        'Schematic arrows represent evaluation relations, not measured trajectories of individual cells.': 'Solid schematic arrows represent GDC evaluation relations; dashed connections add literature context. Neither represents measured trajectories of individual cells.',
        '示意图箭头表示评价关系，不是单个细胞的实测轨迹。': '示意图实线箭头表示 GDC 评价关系，虚线补充文献关系，均不表示单个细胞的实测轨迹。',
    })
    if COMPLETED:
        replacements.update({
            'Saved score and field provenance is retained; native TopicVelo outputs are missing for both settings, and TFvelo is missing for KPTracer. Neither gap is filled using an output from another method or configuration.':
            'Saved score and field provenance is retained. KPTracer TFvelo has now been refitted and native TopicVelo transitions completed; its transition-displacement GDC proxy is identified separately. Spatial TopicVelo results are documented separately. No base-scVelo output is substituted for native TopicVelo transitions.',
            '保留已有分数与流场的来源记录；两个场景均缺少 TopicVelo 原生输出，KPTracer 还缺少 TFvelo。缺口不使用其他方法或配置的输出替代。':
            '保留已有分数与流场的来源记录。KPTracer 的 TFvelo 已重新拟合，TopicVelo 原生转移已补全，其转移位移 GDC 代理单独标明。空间数据 TopicVelo 结果另有记录；不使用基础 scVelo 输出代替原生转移。',
            'while native TopicVelo remains unavailable for KPTracer, and TFvelo is missing for KPTracer. Neither gap is filled using an output from another method or configuration.':
            'and KPTracer now includes a TFvelo refit and completed native TopicVelo transitions. TopicVelo GDC uses the documented transition-displacement proxy; no historical base-scVelo output is substituted.',
            'KPTracer 的 TopicVelo 原生输出仍缺失，KPTracer 还缺少 TFvelo。缺口不使用其他方法或配置的输出替代。':
            'KPTracer 的 TopicVelo 原生转移和 TFvelo 重拟合现已补全。TopicVelo GDC 使用有记录的转移位移代理，不以历史基础 scVelo 输出替代。',
        })
    for node in doc.iter():
        for attr in ('text','tail'):
            value=getattr(node,attr)
            if value:
                for old,new in replacements.items():value=value.replace(old,new)
                setattr(node,attr,value)
