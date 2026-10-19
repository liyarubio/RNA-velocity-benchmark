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
CAPTION_ZH = ('The original blue-bubble encoding is applied to the three original cellDancer scenarios. Panels a-c show cell-level medians, '
    'and d-f show gene-level medians. Each scenario includes ten datasets with noise from 0.2 to 2.0, totaling 30 datasets. '
    'Lighter colors indicate better within-condition ranks. Bubble areas use the common retained-data minimum minus 0.01 to maximum range for display only. '
    'Numbers above columns indicate shared genes with finite values. Grey crosses mark missing results, including nine failed STT outputs in the forward-branching scenario. '
    'biVI dagger-marked results are exploratory fits to continuous inputs. Each condition has one saved dataset and cannot be treated as a mean over independent replicates; '
    'method preprocessing and retained gene sets differ. Retained raw scores and the main-text overall ranking are unchanged.')

REPLACEMENTS = {
    '; ten additional backward-branching datasets are analyzed separately.': '.',
    ', with ten backward-branching datasets analyzed separately.': '.',
    ' Ten saved backward-branching datasets form a separate extension. They are not added to the original 72-dataset design or to the main composite ranking.': '',
    'Ten additional saved backward-branching datasets are shown separately as an extension and are not included in the original 72-dataset design or main-text composite ranking.': '',
    'the noise series and separate backward-branching datasets': 'the complete noise series for these three scenarios',
    'the full noise series and separately analyzed backward-branching datasets': 'the full noise series for three scenarios',
    'The present supplementary re-evaluation retains these original 30 datasets and separately includes ten saved backward-branching datasets. Their noise mapping was not verified, and they are shown by series index. These additional datasets are not included in the original 72-dataset total (42 dyngen plus 30 cellDancer datasets).': 'The present supplementary re-evaluation includes the original 30 cellDancer datasets. Together with the 42 dyngen datasets, the simulation design contains 72 datasets.',
    'This Extended Data reevaluation retains the original 30 datasets and separately shows ten existing backward-branching datasets. Their noise correspondence remains unverified, so they are labeled by sequence index. These additional datasets are not counted in the original 72 simulated datasets (42 dyngen and 30 cellDancer).': 'This Extended Data reevaluation includes the original 30 cellDancer datasets; together with 42 dyngen datasets, the simulation design contains 72 datasets.',
    'all 40 saved datasets': 'all 30 retained datasets',
    'all 40 datasets': 'all 30 retained datasets',
    'all 40 existing datasets': 'all 30 retained datasets',
    'all 40 datasets': 'all 30 retained datasets',
    '31 of 40 datasets': '21 of 30 datasets',
    '31/40 datasets': '21/30 datasets',
    'The original 30 cellDancer datasets and ten additional backward-branching datasets': 'The original 30 cellDancer datasets',
    'The full noise series and separate backward-branching extension remain': 'The full noise series for the three original scenarios is shown',
    'The full noise series and backward-branching extension are retained in': 'The full noise series for the three original scenarios is shown in',
    'cellDancer conditions and backward-branching extension': TITLE,
    'Extended Data Fig. 5c shows the three original scenarios at noise 0.2, and Extended Data Fig. 7 extends the comparison to the complete noise series for these three scenarios.': 'Extended Data Fig. 5c shows the three original scenarios at noise 0.2; the complete noise-series results are retained in the source-data summary.',
    'Extended Data Fig. 5c shows the three original scenarios at noise 0.2; Extended Data Fig. 7 extends the comparison to the full noise series for these scenarios.': 'Extended Data Fig. 5c shows the three original scenarios at noise 0.2; results for the full noise series are retained in the source-data summaries.',
    'The full noise series for the three original scenarios is shown in Extended Data Fig. 7.': 'The full noise-series results are retained in the source-data summary.',
    'The full noise series for the three original scenarios is shown in Extended Data Fig. 7.': 'Results for the full noise series are retained in the source-data summaries.',
    'STT is unavailable for forward branching; missing results remain missing. BiVI (dagger) is exploratory because it was fitted to continuous simulated U/S values. ': 'Methods without valid results are omitted from the corresponding subplot. ',
    'Cells and genes are not independent simulation replicates. Missing STT results for forward branching remain unavailable. biVI dagger-marked fits to continuous U/S remain exploratory; ': 'Cells and genes are not independent simulation replicates; methods without valid results are omitted from the corresponding panels.',
}

def scope_text(value):
    if not isinstance(value, str):
        return value
    if value.startswith('Original-style blue bubble encoding applied to the additional cellDancer condition summary.'):
        return CAPTION_EN + ('\n' + CAPTION_ZH if 'The original blue-bubble encoding' in value else '')
    if value.startswith('The original blue-bubble encoding is applied to the cellDancer condition extension.'):
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
