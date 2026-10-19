"""Verify numerical artifacts, selected prose changes and all current readers."""
from pathlib import Path
import base64,csv,hashlib,json,zipfile,sys
from urllib.parse import urlsplit,unquote
from lxml import html,etree
HERE=Path(__file__).resolve().parent;R=HERE.parent;M=R.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
report={'status':'pass','readers':{},'figures':{},'allowed_prose_changes':['P051','P120']}
paths=[R/'Manuscript_0528_minimal_bilingual.html',R/'Manuscript_0528_minimal_bilingual_standalone.html',
       M/'audit_manuscript_20260930/manuscript_nature_methods_full.html',M/'audit_manuscript_20260930/manuscript_nature_methods_standalone.html']
sys.path.insert(0,str(R))
from figure1_citations import baseline_text
expected={r['id']:r['revised'] for r in json.loads((R/'paragraphs.json').read_text()) if r['id']!='P183'}
canonical=None
for p in paths:
    d=html.document_fromstring(p.read_bytes(),parser=html.HTMLParser(encoding='utf-8',huge_tree=True))
    assert len(d.xpath('//*[@id="all-figures-gallery"]/figure'))==7
    assert d.get('data-stability-completion')=='20261008'
    ids=d.xpath('//@id');assert len(ids)==len(set(ids))
    english=d.xpath('//*[contains(concat(" ",@class," ")," en ")]')
    assert {n.getparent().get('id'):baseline_text(n) for n in english}==expected
    # A concurrent algorithm/citation revision changed the earlier snapshot.
    # Validate against the current canonical baseline, preserving that revision.
    prose={pid:{lang:[n.text_content() for n in d.get_element_by_id(pid).xpath('./p[@lang="'+lang+'"]')]
                 for lang in ['en','zh-CN']} for pid in expected}
    if canonical is None:canonical=prose
    else:assert prose==canonical,'Readers disagree with current canonical prose'
    for n in range(1,8):
        f=d.get_element_by_id('ppt-figure-'+str(n));img=f.xpath('.//img')[0]
        target=R/f'Figures_original_style/Figure{n}/Figure{n}_original_style.png'
        if img.get('src').startswith('data:'):assert base64.b64decode(img.get('src').split(',',1)[1])==target.read_bytes()
        else:assert (p.parent/unquote(img.get('src'))).resolve()==target.resolve()
    for value in d.xpath('//@href|//@src'):
        u=urlsplit(value)
        if u.scheme or u.netloc:continue
        if u.path:assert (p.parent/unquote(u.path)).exists(),(p,value)
        elif u.fragment:assert unquote(u.fragment) in ids,(p,value)
    relevant=' '.join(d.get_element_by_id(pid).text_content() for pid in ['P051','P120','ppt-figure-5','ppt-figure-6'])
    assert '17 methods' not in relevant and 'stability data are unavailable' not in relevant
    assert 'transition' in relevant and '472' in relevant
    report['readers'][str(p)]={'sha256':sha(p),'seven_figure_gallery':True,'images_current':True,'links_valid':True,'current_canonical_baseline_matches':True}
for folder in ['Figures','Figures_original_style']:
    for n in [5,6]:
        name=f'Figure{n}_'+('revised' if folder=='Figures' else 'original_style')
        p=R/f'{folder}/Figure{n}/{name}.svg'
        d=etree.parse(str(p),parser=etree.XMLParser(huge_tree=True));ids=d.xpath('//@id')
        assert len(ids)==len(set(ids))
        images=len(d.findall('.//{http://www.w3.org/2000/svg}image'))
        if folder=='Figures':assert images==0
        assert p.with_suffix('.png').is_file() and p.with_suffix('.pdf').is_file()
        report['figures'][str(p.relative_to(R))]={'sha256':sha(p),'embedded_images':images}
for p in list((R/'Figures/Figure5').glob('*.layout.json'))+list((R/'Figures/Figure6').glob('*.layout.json'))+list((R/'Figures_original_style/Figure5').glob('*.layout.json')):
    assert not json.loads(p.read_text())['out_of_bounds_text'],p
with zipfile.ZipFile(R/'Source_data_and_Tables_20261008.zip') as z:
    assert len(z.namelist())==9
    for name in z.namelist():assert z.read(name)==(R/name).read_bytes()
report['xlsx_downloads']=9
with (M/'Figure5_sources/data/Figure5c_reads_velocity_cosine_similarity.csv').open() as f:rows=list(csv.DictReader(f))
assert len(rows)==525 and len({r['Method'] for r in rows})==21
report['figure5_pairwise_rows']=525
unchanged=[]
for row in json.loads((HERE/'before_manifest.json').read_text()):
    p=Path(row['path'])
    if '/Figures/Figure3/' in str(p) or '/Figures/Figure4/' in str(p) or p.name in ['Fig3_Accuracy.xlsx','Fig4_ReDeeM.xlsx']:
        assert sha(p)==row['sha256'],p
        unchanged.append(str(p))
report['unrelated_artifacts_preserved']=unchanged
report['merge']=json.loads((HERE/'merge_validation.json').read_text())
report['compatibility']=json.loads((HERE/'compatibility_validation.json').read_text())
(HERE/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print('PASS: 525 rows, four synchronized readers, current gallery images and canonical baseline, unchanged unrelated figures, nine XLSX downloads')
