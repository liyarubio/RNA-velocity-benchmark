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
    'a，全部 1292 个 ST_raw spot的两种空间着色视图使用相同坐标：一张保留原始全部 12 类细胞配色，另一张将目标六类的 608 个 spot按主图 5a 配色高亮，其他 684 个 spot统一为灰色。直接复用主图 5a 的参考示意图，五条箭头表示预定义评价关系，不是实测空间轨迹。'
    'b，GDC 从高到低排列。沿用修正 eval2 GDC，并补算 DeepVelo(VAE) 已保存速度及当前空间投影的评分（-0.493）。eval2 以六类参考类型的 608 个 spot为权重分母；用同一流程复现三个已有分数，误差小于 1e-6。'
    'TFvelo 保留独立评价及其原有分母、预处理，不将这些柱形描述为统一协议新排名。'
    'c，19 个有效方法流场按与 b 相同的 GDC 降序从左到右、从上到下排列，使用相同坐标，非目标类型统一置灰。TopicVelo 从保存主题拟合补跑原生随机转移推断，剔除无可变 U/S 信号的基因和不足两个有效基因的主题，ST_raw 整合 11 个主题，ST_small 整合 5 个主题，替代历史仅用末主题的配置。高维方向采用原生转移位移代理 T @ Ms - Ms，空间流场直接由原生转移矩阵投影；GDC 分别为 -0.419 和 0.484。ST_raw 有 1 个 spot 方向未定义，方向评价中排除该 spot，参考分母仍为 608 spots。其他方法方向保持原值；a–c 合并为单页。'
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
        '两个场景均缺少 TopicVelo 原生输出': '空间 TopicVelo 已通过有记录的多主题重推断补全，KPTracer 的 TopicVelo 原生输出仍缺失',
        '608-cell ST_small subset': '608-spot ST_small subset',
        'all 1,292 cells in the supplied ST_raw input': 'all 1,292 spots in the supplied ST_raw input',
        '608 个细胞的 ST_small 子集': '608 个 spot 的 ST_small 子集',
        'ST_raw 输入中的全部 1292 个细胞': 'ST_raw 输入中的全部 1292 个 spot',
        'Centroid arrows illustrate the relations used for evaluation, not measured trajectories of individual cells.':
        'The spatial reference relations are shown in a separate labeled schematic; the KPTracer centroid arrows illustrate its evaluation relations. Neither represents measured trajectories of individual cells.',
        '质心之间的箭头表示评价所依据的关系，不是单个细胞的实测轨迹。':
        '空间数据的参考关系以独立且带名称的分支示意图展示；KPTracer 的质心箭头表示其评价关系。两者均不表示单个细胞的实测轨迹。',
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
