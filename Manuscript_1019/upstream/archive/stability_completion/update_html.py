"""Localized canonical HTML update preserving manuscript baseline markers."""
from pathlib import Path
import sys,base64
from lxml import html
R=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(R))
from revised_figure_assets import TEXT,CAPTIONS
for name in ['Manuscript_0528_minimal_bilingual.html','Manuscript_0528_minimal_bilingual_standalone.html']:
    p=R/name;doc=html.document_fromstring(p.read_bytes(),parser=html.HTMLParser(encoding='utf-8',huge_tree=True))
    for pid in ['P051','P120']:
        for lang,text in TEXT[pid].items():
            parent=doc.get_element_by_id(pid).xpath('./p[@lang="'+lang+'"]')[0]
            spans=parent.xpath('./span[@data-figure-data-update]')
            assert len(spans)==1
            span=spans[0]
            for child in list(span):span.remove(child)
            span.text=text
    for n in [5,6]:
        fig=doc.get_element_by_id('ppt-figure-'+str(n))
        fig.xpath('.//figcaption')[0].text=CAPTIONS[n]
        img=fig.xpath('.//img')[0];image=R/f'Figures/Figure{n}/Figure{n}_revised.png'
        img.set('src','data:image/png;base64,'+base64.b64encode(image.read_bytes()).decode() if 'standalone' in name else str(image.relative_to(R)))
    doc.set('data-stability-completion','20261008')
    note=doc.get_element_by_id('figure5-scope-note')
    note.text='5a：参考流向；5b：细胞数与运行时间；5c：reads 降采样稳定性。已补入 cell2fate 和 TopicVelo*；后者使用原生联合转移矩阵在固定全深度基因坐标中的投影方向，不使用基础 scVelo velocity。Figure 6 取 0.2–0.8× 相对 1.0× 的四档均值，稳定性不计入综合分。'
    p.write_bytes(html.tostring(doc,encoding='utf-8',doctype='<!doctype html>'))
print('Updated canonical bilingual HTML and embedded Figure5/6 previews')
