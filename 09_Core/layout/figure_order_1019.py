"""Current display numbering; legacy source IDs and numerical assets stay stable."""
import re

LEGACY_ORDER = [1, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18]
NUMBER = {old: new for new, old in enumerate(LEGACY_ORDER, 1)}
ITEM = r'\d+[a-z]?(?:[––-][a-z])?(?:,\s*[a-z](?![a-z]))*'
ED_PATTERN = re.compile(
    '(Extended Data (?:Figs?\\.|Figures?)\\s*)'
    r'(' + ITEM + '(?:(?:\\s+and\\s+|\\s*(?:and|,)\\s*|,\\s*(?:and\\s+)?|[–-])' + ITEM + r')*)'
)

def renumber_text(text):
    # These are provenance statements about older layouts, not current callouts.
    replacements = {
        'The original Extended Data Fig. 4 layout is retained': 'The original simulation-comparison layout is retained',
        'restoring the original Extended Data Fig. 7 design': 'restoring the original metabolic-labeling design',
        'restoring the original Extended Data Fig. 8 design': 'restoring the original lineage-tracing design',
        'plain node-and-arrow style of Fig. 15': 'plain node-and-arrow style of Extended Data Fig. 15',
    }
    for a, b in replacements.items():
        text = text.replace(a, b)
    def replace(m):
        nums = re.sub(r'\d+', lambda n: str(NUMBER[int(n[0])]), m[2])
        return m[1] + nums
    return ED_PATTERN.sub(replace, text)

def apply_figure_order(doc):
    # The original overview belongs at the first discussion of those datasets.
    node = doc.find(id='P035')
    if node:
        for lang, old, new in [
            ('en', '(Supplementary Table 3)', '(Extended Data Fig. 4; Supplementary Table 3)'),
            ('zh-CN', '(Supplementary Table 3)', '(Extended Data Fig. 4; Supplementary Table 3)'),
        ]:
            p = node.find(lang=lang)
            if p:
                for t in list(p.find_all(string=True)):
                    if old in str(t):
                        t.replace_with(str(t).replace(old, new, 1))
                        break
    for t in list(doc.find_all(string=True)):
        if t.parent.name in ['script', 'style']:
            continue
        new = renumber_text(str(t))
        if new != str(t):
            t.replace_with(new)
    for tag in doc.find_all(True):
        oldid = tag.get('id', '')
        if re.fullmatch(r'ed\d+', oldid):
            tag['data-legacy-figure'] = oldid
            tag['id'] = 'ed' + str(NUMBER[int(oldid[2:])])
        for attr in ['href', 'src']:
            if not tag.has_attr(attr):
                continue
            value = tag[attr]
            value = re.sub(r'#ed(\d+)\b', lambda m: '#ed' + str(NUMBER[int(m[1])]), value)
            value = re.sub(r'figures/extended/Extended_Data_Fig(\d+)\.(png|svg|pdf)',
                lambda m: f'figures/extended_numbered/Extended_Data_Fig{NUMBER[int(m[1])]:02d}.{m[2]}', value)
            tag[attr] = value
        for attr in ['alt', 'title']:
            if tag.has_attr(attr):
                tag[attr] = renumber_text(tag[attr])
        if tag.name == 'a' and tag.get('href', '').startswith('source_data/Extended_Data_Fig'):
            old = int(re.search(r'Fig(\d+)', tag['href'])[1])
            tag.string = f'Source workbook (current ED{NUMBER[old]:02d}; legacy source ED{old:02d})'
    gallery = doc.find(id='extended-data-gallery')
    if gallery:
        if gallery.h2:
            gallery.h2.string = 'Extended Data Figures 1-12'
        for num in range(1, len(NUMBER) + 1):
            section = gallery.find(id=f'ed{num}')
            assert section is not None, num
            gallery.append(section.extract())

def figure_occurrences(doc, lang='en'):
    """Callouts in narrative/methods, excluding legends and references."""
    occurrences = []
    for block in doc.select('section.pair'):
        if block.has_attr('data-reference-number'):
            continue
        p = block.find(lang=lang, recursive=False)
        if not p:
            continue
        text = p.get_text(' ', strip=True)
        for match in ED_PATTERN.finditer(text):
            for number in re.findall(r'\d+', match[2]):
                occurrences.append({'paragraph': block['id'], 'number': int(number), 'callout': match[0]})
    return occurrences
