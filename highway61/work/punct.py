"""Punctuation made from the letters of each font (lines and circles only, so exact):
  period       the i's dot, set on the baseline
  comma        the y's tail below the baseline, raised so its top is level with the period's top
  quoteright / quotesingle   the comma raised so its top is on the cap line
  quoteleft    the comma turned through 180 degrees, its top on the cap line
  exclam       the i turned upside down, its dot on the baseline
  quotedblright / quotedbl, quotedblleft   two of those, a stroke apart
  hyphen       the H's crossbar between the stems, at the crossbar's height
make(glyph_polygon) -> {glyph name: strokeop segment contours}; glyph_polygon(name) gives a shapely geometry."""
import numpy as np
from shapely.geometry import Polygon, box, LineString
from shapely import affinity
import lcoutline, lcgeom
CAP=690.0
def _segs(geom,prims=()):
    return [[tuple(s) for s in segs] for segs in lcoutline.contours(geom,list(prims))]
def _dot(gi):
    """centre and radius of the i's dot: its highest piece"""
    parts=sorted(getattr(gi,'geoms',[gi]),key=lambda p:-p.bounds[3])
    x0,y0,x1,y1=parts[0].bounds; r=(x1-x0)/2; return r
def _crossbar(gh):
    x0,y0,x1,y1=gh.bounds; rows=[]
    for y in np.arange(y0+2,y1-2,1.0):
        s=gh.intersection(LineString([(x0-10,y),(x1+10,y)]))
        n=len(getattr(s,'geoms',[s])); rows.append((y,n,s))
    bar=[y for y,n,s in rows if n==1]
    yb,yt=min(bar),max(bar)
    s=[s for y,n,s in rows if n==2 and y>yt][0]
    a,b=sorted(getattr(s,'geoms',[s]),key=lambda q:q.bounds[0])
    # exact bar edges: refine to the outline's horizontal edges
    ys=sorted({round(q[1],4) for q in np.array(gh.exterior.coords) if yb-2<q[1]<yt+2})
    yb=min(ys,key=lambda v:abs(v-yb)); yt=min(ys,key=lambda v:abs(v-yt))
    return a.bounds[2],b.bounds[0],yb,yt
def make(poly):
    out={}
    r=_dot(poly('i'))
    lcgeom.REG.clear(); P=lcgeom.ellipse(r,r,r,r)
    out['period']=_segs(P,[(r,r,r,r)])
    gy=poly('y')
    tail=gy&box(-1000,-1000,3000,0.0)
    tail=max(getattr(tail,'geoms',[tail]),key=lambda p:p.area)
    comma=affinity.translate(tail,0,2*r)                 # top level with the period's top
    out['comma']=_segs(comma)
    apo=affinity.translate(comma,0,CAP-comma.bounds[3]); out['quoteright']=_segs(apo); out['quotesingle']=out['quoteright']
    ql=affinity.rotate(comma,180,origin=comma.centroid)
    ql=affinity.translate(ql,0,CAP-ql.bounds[3]); out['quoteleft']=_segs(ql)
    # exclamation point: the i turned upside down, its dot on the baseline
    gi=poly('i'); ex=affinity.rotate(gi,180,origin=gi.centroid); ex=affinity.translate(ex,0,-ex.bounds[1])
    dot=min(getattr(ex,'geoms',[ex]),key=lambda p:p.bounds[1]); x0,y0,x1,y1=dot.bounds
    out['exclam']=_segs(ex,[((x0+x1)/2,(y0+y1)/2,(x1-x0)/2,(y1-y0)/2)])
    w=apo.bounds[2]-apo.bounds[0]; xl,xr,yb,yt=_crossbar(poly('H')); T=yt-yb
    gap=w+0.5*T                                            # side by side, half a stroke apart
    out['quotedblright']=_segs(apo|affinity.translate(apo,gap,0)); out['quotedbl']=out['quotedblright']
    out['quotedblleft']=_segs(ql|affinity.translate(ql,gap,0))
    out['hyphen']=_segs(box(xl,yb,xr,yt))
    return out
