"""Active cellDancer scope: the original three scenarios, 30 datasets.

Frozen V2 inputs and model outputs remain historical provenance. This module
applies the requested exclusion when building current readers and workbooks.
"""
from ed6_revision import revise_text
FAMILIES = ('Mono-kinetic', 'Transcriptional_boost', 'Multi-forward_branching')
TITLE = 'cellDancer conditions across noise levels'
CAPTION_EN = ('Original-style blue bubble encoding of the three original cellDancer scenarios. '
    'a-c, Median cell-level cosine; d-f, median gene-level cosine. Each scenario contains ten '
    'noise levels (0.2-2.0), for 30 datasets in total. Lighter blue indicates better within-condition '
    'rank, while circle area uses a common min-minus-0.01/max normalization of the retained median '
    'cosine values for display only. Shared finite gene counts appear above columns. Grey crosses '
    'show missing results, including nine failed forward-branching STT outputs. BiVI (dagger) '
    'remains exploratory on continuous inputs. Each condition is one saved dataset, not an '
    'independent-replicate average; preprocessing and retained gene sets differ across datasets. '
    'Retained raw scores and the main composite ranking are unchanged.')
CAPTION_ZH = ('将原版蓝色气泡编码用于 cellDancer 的三类原始场景。a–c 为细胞级中位数，'
    'd–f 为基因级中位数。每类场景包含噪声 0.2–2.0 的十个数据集，共 30 个数据集。'
    '颜色越浅表示条件内排名越好，圆面积按保留数据的统一最小值减 0.01 至最大值归一化，仅用于显示。'
    '列上方数字为共有有限值基因数。灰叉表示缺失，包括九个前向分支 STT 失败输出。'
    'biVI† 为连续输入上的探索性拟合。每个条件只有一个已保存数据集，不能当作独立重复均值；'
    '方法预处理及保留基因范围存在差异。保留的原始得分和主文综合排名不变。')

REPLACEMENTS = {
    '; ten additional backward-branching datasets are analyzed separately.': '.',
    '，另将十个反向分支数据集单列分析。': '。',
    ' Ten saved backward-branching datasets form a separate extension. They are not added to the original 72-dataset design or to the main composite ranking.': '',
    '另有十个已保存的反向分支数据集，单列为扩展，不并入原始 72 数据集设计或主文综合排名。': '',
    'the noise series and separate backward-branching datasets': 'the complete noise series for these three scenarios',
    '完整噪声序列和单列的反向分支数据': '三类场景的完整噪声序列',
    'The present supplementary re-evaluation retains these original 30 datasets and separately includes ten saved backward-branching datasets. Their noise mapping was not verified, and they are shown by series index. These additional datasets are not included in the original 72-dataset total (42 dyngen plus 30 cellDancer datasets).': 'The present supplementary re-evaluation includes the original 30 cellDancer datasets. Together with the 42 dyngen datasets, the simulation design contains 72 datasets.',
    '本次附图重新评价保留上述原始 30 个数据集，另行展示已有的十个反向分支数据集；后者噪声对应关系尚未验证，因此按序列编号展示。这十个附加数据集不计入原稿 72 个模拟数据集的口径（42 个 dyngen 加 30 个 cellDancer）。': '本次附图重新评价包含上述原始 30 个 cellDancer 数据集；加上 42 个 dyngen 数据集，模拟设计共包含 72 个数据集。',
    'all 40 saved datasets': 'all 30 retained datasets',
    'all 40 datasets': 'all 30 retained datasets',
    '全部 40 个已有数据集': '保留的全部 30 个数据集',
    '全部 40 个数据集': '保留的全部 30 个数据集',
    '31 of 40 datasets': '21 of 30 datasets',
    '31/40 个数据集': '21/30 个数据集',
    'The original 30 cellDancer datasets and ten additional backward-branching datasets': 'The original 30 cellDancer datasets',
    'The full noise series and separate backward-branching extension remain': 'The full noise series for the three original scenarios is shown',
    '完整噪声序列与反向分支扩展保留在': '三类原始场景的完整噪声序列展示在',
    'cellDancer conditions and backward-branching extension': TITLE,
    'Extended Data Fig. 5c shows the three original scenarios at noise 0.2, and Extended Data Fig. 7 extends the comparison to the complete noise series for these three scenarios.': 'Extended Data Fig. 5c shows the three original scenarios at noise 0.2; the complete noise-series results are retained in the source-data summary.',
    '扩展数据图 5c 展示噪声 0.2 的三类原始场景，扩展数据图 7 将比较扩展到三类场景的完整噪声序列。': '扩展数据图 5c 展示噪声 0.2 的三类原始场景；完整噪声序列结果保留在源数据汇总中。',
    'The full noise series for the three original scenarios is shown in Extended Data Fig. 7.': 'The full noise-series results are retained in the source-data summary.',
    '三类原始场景的完整噪声序列展示在扩展数据图 7。': '完整噪声序列结果保留在源数据汇总中。',
    'STT is unavailable for forward branching; missing results remain missing. BiVI (dagger) is exploratory because it was fitted to continuous simulated U/S values. ': 'Methods without valid results are omitted from the corresponding subplot. ',
    '细胞和基因不作为独立模拟重复；前向分支的 STT 缺失保留为空。biVI† 在连续 U/S 上的拟合仍属探索性结果；': '细胞和基因不作为独立模拟重复；无有效结果的方法不在对应子图中展示。',
}

def scope_text(value):
    if not isinstance(value, str):
        return value
    if value.startswith('Original-style blue bubble encoding applied to the additional cellDancer condition summary.'):
        return CAPTION_EN + ('\n' + CAPTION_ZH if '将原版蓝色气泡编码' in value else '')
    if value.startswith('将原版蓝色气泡编码用于 cellDancer 条件扩展。'):
        return CAPTION_ZH
    for old, new in REPLACEMENTS.items():
        value = value.replace(old, new)
    return revise_text(value)

def current_rows(rows):
    """Filter excluded observation/provenance rows; retain surviving values exactly."""
    for row in rows:
        if any(isinstance(v, str) and
               (v.startswith('Multi-backward_branching') or v == 'backward branching')
               for v in row):
            continue
        yield tuple(scope_text(v) for v in row)

def scope_workbook(workbook):
    for sheet in workbook:
        rows = list(sheet.values)
        filtered = list(current_rows(rows))
        if rows != filtered:
            sheet.delete_rows(1, sheet.max_row)
            for row in filtered:
                sheet.append(row)

def scope_document(document):
    for node in document.xpath('//*[@id="simulation-extension"]'):
        node.getparent().remove(node)
    for node in document.iter():
        if not isinstance(node.tag, str):
            continue
        node.text = scope_text(node.text)
        node.tail = scope_text(node.tail)
        href = node.get('href', '')
        for filename in ('dataset_inventory.csv', 'dataset_summary.csv', 'method_coverage.csv', 'evaluation_protocol.json'):
            if href == 'Source_data/cellDancer_simulations/' + filename:
                node.set('href', 'Source_data/cellDancer_current_30/' + filename)
