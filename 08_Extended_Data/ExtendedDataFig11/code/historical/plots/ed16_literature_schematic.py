"""Redraw Yang et al. Cell 2022 Fig. 5I; line style encodes GDC inclusion.

Fig. 5I groups Gastric/Endoderm and Early EMT-1/2. The former group-level
progression is expanded to the saved Endoderm reference and a contextual
Gastric branch; the latter is split to distinguish the evaluated EMT-1 state.
Unconnected mesenchymal context nodes are omitted for visual clarity.
"""
import pymupdf as fitz

DOI = '10.1016/j.cell.2022.04.015'
SOURCE_URL = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9452598/figure/F5/'

def draw(page, text, palette, reference_edges):
    nodes = {
        'AT2-like': (264, 88), 'Early gastric': (205, 111),
        'AT1-like': (264, 119), 'Gastric-like': (218, 151),
        'Endoderm-like': (194, 180), 'High plasticity': (264, 151),
        'Lung progenitor-like': (321, 151), 'Pre-EMT': (264, 205),
        'Early EMT-1': (244, 237), 'Early EMT-2': (287, 237),
    }
    labels = {
        'AT2-like': ['AT2-like'], 'Early gastric': ['Early', 'gastric'],
        'AT1-like': ['AT1-like'], 'Gastric-like': ['Gastric-', 'like'],
        'Endoderm-like': ['Endoderm-', 'like'],
        'High plasticity': ['High', 'plasticity'],
        'Lung progenitor-like': ['Lung progenitor-', 'like'],
        'Pre-EMT': ['Pre-EMT'], 'Early EMT-1': ['Early EMT-1'],
        'Early EMT-2': ['Early EMT-2'],
    }
    # White background and simple nodes/edges, matching Fig. 15.

    reference_edges = [tuple(e) for e in reference_edges]
    context = [
        ('Gastric-like','Pre-EMT','group-level progression'),
        ('High plasticity','Pre-EMT','likely intermediate states'),
        ('Lung progenitor-like','Pre-EMT','rare examples'),
        ('Pre-EMT','Early EMT-2','split of Early EMT-1/2'),
        ('High plasticity','Lung progenitor-like','undirected coupling'),
    ]
    records = []
    for src,dst,meaning in [(a,b,'GDC reference') for a,b in reference_edges] + context:
        solid = (src,dst) in reference_edges
        directed = meaning != 'undirected coupling'
        a,b = fitz.Point(nodes[src]),fitz.Point(nodes[dst])
        # Start below text labels; all other edges run between node boundaries.
        if src == 'AT2-like':
            a += fitz.Point(0,5)
            bend = fitz.Point(nodes[dst][0],112)
            end = b-fitz.Point(0,6)
            points = [a,bend,end]
        elif src == 'Gastric-like' and dst == 'Pre-EMT':
            points=[fitz.Point(223,151),fitz.Point(235,151),fitz.Point(242,182),fitz.Point(260,201)]
        elif src == 'Lung progenitor-like' and dst == 'Pre-EMT':
            points=[fitz.Point(318,155),fitz.Point(290,158),fitz.Point(279,181),fitz.Point(268,201)]
        else:
            v=(b-a)/abs(b-a);a+=v*5;b-=v*5
            points=[a,b]
        color=(.13,.13,.13) if solid else (.42,.42,.42)
        width=.95 if solid else .7
        dash=None if solid else '[2.6 2] 0'
        page.draw_polyline(points,color=color,width=width,dashes=dash)
        if directed:
            end=points[-1];v=(end-points[-2])/abs(end-points[-2]);q=end-v*3.3;n=fitz.Point(-v.y,v.x)*1.35
            page.draw_polyline([q+n,end,q-n],color=color,width=width)
        records.append({'source':src,'target':dst,'GDC':solid,
                        'directed':directed,'literature_meaning':meaning,
                        'path':[list(p) for p in points],
                        'line_style':'solid' if solid else 'dashed'})
    for name,(x,y) in nodes.items():
        color=palette.get(name)
        fill=tuple(int(color[i:i+2],16)/255 for i in (1,3,5)) if color else (1,1,1)
        page.draw_circle((x,y),3.5,color=None if color else (.5,.5,.5),fill=fill,width=.65)
        lines=labels[name]
        # Put branching-node labels above or beside the outgoing edges.
        if name == 'Pre-EMT':
            text(x+11,y+2,'Pre-EMT',5.9)
            continue
        if name in ['AT2-like','Early gastric','AT1-like','High plasticity']:
            first=y-8-(len(lines)-1)*6.6
        elif name in ['Early EMT-1','Early EMT-2']:
            first=y+10
        elif name == 'Lung progenitor-like':
            first=y+24
        else:
            first=y+10
        for i,line in enumerate(lines):
            text(x,first+i*6.6,line,5.9,center=True)
    # These styles deliberately encode evaluation membership, not the source
    # article's distinct solid/dotted/gray levels of biological evidence.
    page.draw_line((187,267),(202,267),color=(.13,.13,.13),width=.95)
    text(207,269,'Used for GDC',5.8)
    page.draw_line((273,267),(288,267),color=(.42,.42,.42),width=.7,dashes='[2.6 2] 0')
    text(292,269,'Other relations',5.8)
    return nodes, {'source_figure':'Yang et al., Cell (2022), Fig. 5I',
                   'doi':DOI,'url':SOURCE_URL,'edges':records,
                   'context_only_nodes':[n for n in nodes if n not in palette],
                   'hollow_nodes':'not present in the supplied 754-cell subset',
                   'line_style_semantics':'GDC inclusion, not original evidence strength',
                   'omitted_unconnected_context_nodes':['Mesenchymal-1','Mesenchymal-2'],
                   'line_legend_y_pt':267,'background_regions':False,'style_reference':'Fig. 15'}
