"""Align saved vector fields by their 754 cell coordinates, enlarge markers."""
import math
import numpy as np
import pymupdf as fitz

POINT_RADIUS = 1.30
CELL_BOUNDS = (6., 8., 99., 85.)

def rgb(s):
    return tuple(int(s[i:i+2],16)/255 for i in (1,3,5))

def replay(shape, drawing, matrix):
    for item in drawing['items']:
        kind=item[0]
        if kind=='l':shape.draw_line(item[1]*matrix,item[2]*matrix)
        elif kind=='c':shape.draw_bezier(*(v*matrix for v in item[1:]))
        elif kind=='re':
            r=item[1];pts=[r.tl,r.tr,r.br,r.bl]
            if item[2]<0:pts=pts[::-1]
            shape.draw_polyline([p*matrix for p in pts+[pts[0]]])
        elif kind=='qu':shape.draw_quad(item[1]*matrix)
        else:raise ValueError(kind)
    width=(drawing['width'] or 0)*math.sqrt(abs(matrix.a*matrix.d))
    # PDF width 0 is a hairline. Shape.finish(width=0) instead removes stroke
    # color, so emit the PDF operator explicitly while retaining its color.
    if drawing['color'] is not None and width==0:
        shape.draw_cont='0 w\n'+shape.draw_cont
        width=1
    shape.finish(color=drawing['color'],fill=drawing['fill'],
                 width=width,
                 lineCap=int(max(drawing['lineCap'] or [0])),lineJoin=int(drawing['lineJoin'] or 0),
                 dashes=drawing['dashes'],closePath=bool(drawing['closePath']),
                 even_odd=bool(drawing['even_odd']),
                 fill_opacity=drawing['fill_opacity'] if drawing['fill_opacity'] is not None else 1,
                 stroke_opacity=drawing['stroke_opacity'] if drawing['stroke_opacity'] is not None else 1)

def prepare(source_page,cells,palette):
    drawings=source_page.get_drawings()
    centers=np.array([[(d['rect'].x0+d['rect'].x1)/2,(d['rect'].y0+d['rect'].y1)/2] for d in drawings])
    colors={name:rgb(color) for name,color in palette.items()}
    is_cell=np.array([bool(d['fill']) and d['rect'].width<4 and d['rect'].height<4 and
                     any(np.max(np.abs(np.array(d['fill'])-c))<.005 for c in colors.values()) for d in drawings])
    xy=np.array([[r[1],r[2]] for r in cells]);lo=xy.min(0);span=np.ptp(xy,axis=0)
    l,t,r,b=CELL_BOUNDS
    positions=np.column_stack([l+(xy[:,0]-lo[0])/span[0]*(r-l),b-(xy[:,1]-lo[1])/span[1]*(b-t)])
    return drawings,centers,is_cell,positions,colors

def render(clip,cells,prepared):
    drawings,centers,is_cell,positions,colors=prepared
    selected=(centers[:,0]>=clip.x0)&(centers[:,0]<=clip.x1)&(centers[:,1]>=clip.y0)&(centers[:,1]<=clip.y1)
    point_ids=np.flatnonzero(selected&is_cell)
    assert len(point_ids)==754
    points=centers[point_ids];lo=points.min(0);hi=points.max(0)
    l,t,r,b=CELL_BOUNDS;sx=(r-l)/(hi[0]-lo[0]);sy=(b-t)/(hi[1]-lo[1])
    matrix=fitz.Matrix(sx,0,0,sy,l-lo[0]*sx,t-lo[1]*sy)
    mapped=np.array([list(fitz.Point(p)*matrix) for p in points])
    # The saved scatter must really be the same cells, separately for each type.
    errors=[]
    for name,color in colors.items():
        old=np.array([mapped[i] for i,j in enumerate(point_ids) if np.max(np.abs(np.array(drawings[j]['fill'])-color))<.005])
        new=positions[[row[3]==name for row in cells]]
        assert len(old)==len(new),(name,len(old),len(new))
        dist=np.sqrt(((old[:,None,:]-new[None,:,:])**2).sum(2).min(1))
        assert dist.max()<.003,(name,dist.max())
        # Tiny source-coordinate roundoff may make exact duplicates ambiguous.
        errors.append(float(dist.max()))
    doc=fitz.open();page=doc.new_page(width=105,height=98)
    shape=page.new_shape()
    for name,color in colors.items():
        for p in positions[[row[3]==name for row in cells]]:
            shape.draw_circle(p,POINT_RADIUS)
            shape.finish(color=None,fill=color)
    shape.commit()
    field_ids=[int(i) for i in np.flatnonzero(selected&~is_cell)
               if drawings[i]['fill']!=(1.,1.,1.)]
    # All selected non-background paths are copied with the same affine map.
    shape=page.new_shape()
    for i in field_ids:replay(shape,drawings[i],matrix)
    shape.commit()
    return doc,{'source_cell_bounds':[float(lo[0]),float(lo[1]),float(hi[0]),float(hi[1])],
                'target_cell_bounds':list(CELL_BOUNDS),'source_to_panel_matrix':list(matrix),
                'cells':754,'point_radius_pt':POINT_RADIUS,
                'source_coordinate_max_error_pt':max(errors),'vector_path_count':len(field_ids),
                'source_vector_indices':field_ids}
