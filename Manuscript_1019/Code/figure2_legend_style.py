"""Figure 2 caption typography; preserve all caption text and scientific artwork."""
import re
from copy import deepcopy
from xml.etree import ElementTree as ET

TITLE = 'Fig. 2 | Schematic overview of the benchmarking workflow.'

def format_html(doc):
    figure = doc.find(id='ppt-figure-2')
    caption = figure.find('figcaption') or figure.find('p',recursive=False)
    before = ' '.join(caption.get_text(' ',strip=True).split())
    if not caption.find('strong',attrs={'data-legend-title':'figure2'}):
        for text in list(caption.find_all(string=True)):
            if TITLE in str(text):
                first,rest=str(text).split(TITLE,1)
                strong=doc.new_tag('strong',attrs={'data-legend-title':'figure2'})
                strong.string=TITLE
                text.replace_with(first,strong,rest)
                break
        else: raise AssertionError('Figure 2 caption title missing')
    caption['style']='text-align:justify;text-align-last:left;white-space:normal'
    assert ' '.join(caption.get_text(' ',strip=True).split()) == before

def format_docx(root):
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    q=lambda name:'{'+ns['w']+'}'+name
    for p in root.findall('.//w:p',ns):
        text=''.join(n.text or '' for n in p.findall('.//w:t',ns))
        if not text.startswith(TITLE): continue
        props=p.find('w:pPr',ns)
        if props is None: props=ET.Element(q('pPr'));p.insert(0,props)
        for old in props.findall('w:jc',ns): props.remove(old)
        ET.SubElement(props,q('jc'),{q('val'):'both'})
        remaining=len(TITLE)
        for run in list(p.findall('w:r',ns)):
            if remaining<=0: break
            nodes=run.findall('w:t',ns)
            count=sum(len(n.text or '') for n in nodes)
            if not count: continue
            if count>remaining:
                # Pandoc stores this caption's title in one plain text run.
                assert len(nodes)==1 and len([n for n in run if n.tag!=q('rPr')])==1
                tail=deepcopy(run);tail.find('w:t',ns).text=nodes[0].text[remaining:]
                nodes[0].text=nodes[0].text[:remaining]
                tail.find('w:t',ns).set('{http://www.w3.org/XML/1998/namespace}space','preserve')
                p.insert(list(p).index(run)+1,tail)
            prop=run.find('w:rPr',ns)
            if prop is None:prop=ET.Element(q('rPr'));run.insert(0,prop)
            if prop.find('w:b',ns) is None:ET.SubElement(prop,q('b'))
            remaining-=min(count,remaining)
        assert remaining==0
        assert ''.join(n.text or '' for n in p.findall('.//w:t',ns))==text
        return
    raise AssertionError('Figure 2 DOCX caption missing')

def draw_caption(page,box,caption,fitz):
    """Typeset mixed-weight Arial with justified full lines and a ragged last line."""
    caption=' '.join(caption.split());assert caption.startswith(TITLE)
    fonts={False:fitz.Font(fontfile='/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'),
           True:fitz.Font(fontfile='/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf')}
    for bold in [False,True]:page.insert_font(fontname='LegendBold' if bold else 'LegendRegular',fontbuffer=fonts[bold].buffer)
    tokens=[(word,True) for word in TITLE.split()]
    tokens += [(word,bool(re.fullmatch(r'[a-z],',word))) for word in caption[len(TITLE):].split()]
    size=10;space=fonts[False].text_length(' ',fontsize=size)
    lines=[];line=[];width=0
    for word,bold in tokens:
        length=fonts[bold].text_length(word,fontsize=size)
        if line and width+space+length>box.width:lines.append(line);line=[];width=0
        width+=length+(space if line else 0);line.append((word,bold,length))
    if line:lines.append(line)
    y=box.y0+fonts[False].ascender*size
    assert y+(len(lines)-1)*13+3<page.rect.height-32
    for i,line in enumerate(lines):
        gap=(box.width-sum(w[2] for w in line))/(len(line)-1) if i<len(lines)-1 and len(line)>1 else space
        x=box.x0
        for word,bold,length in line:
            page.insert_text((x,y),word,fontname='LegendBold' if bold else 'LegendRegular',fontsize=size)
            x+=length+gap
        y+=13
    return {'lines':len(lines),'font_size_pt':10,'line_height_pt':13,'title_bold':True,'alignment':'justified; last line left',
        'panel_markers_present':re.findall(r'(?<!\w)[a-z],(?=\s)',caption[len(TITLE):])}
