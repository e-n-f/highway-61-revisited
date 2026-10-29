"""Punctuation made from the letters of each font (lines and circles only, so exact):
  period       the i's dot, set on the baseline
  comma        the y's tail below the baseline, raised so its top is level with the period's top
  quoteright / quotesingle   the comma raised so its top is on the cap line
  quoteleft    the comma turned through 180 degrees, its top on the cap line
  exclam       the i turned upside down, its dot on the baseline
  quotedblright / quotedbl, quotedblleft   two of those, a stroke apart
  hyphen       the H's crossbar between the stems, at the crossbar's height
  slash        (slash(), not from a letter) as drawn in the fractions of the 1958 Interstate signing manual
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
def _strip(tail):
    """the tail as the band between its two long edges, flat at its top and bottom: cut from the y at the baseline,
    it carries a spur of the left arm's junction (and, with the lower case's own offsets, a slight bend), which on
    the comma and quotation marks read as serifs"""
    Q=np.array(tail.exterior.coords)[:-1]; x0,y0,x1,y1=tail.bounds; cx=(x0+x1)/2
    E=[(Q[i],Q[(i+1)%len(Q)]) for i in range(len(Q))]
    longest=lambda side: max([e for e in E if side*((e[0][0]+e[1][0])/2-cx)>0],key=lambda e: np.hypot(*(e[1]-e[0])))
    def x_at(e,y): (a,b)=e; return a[0]+(b[0]-a[0])*(y-a[1])/(b[1]-a[1])
    L,R=longest(-1),longest(1)
    return Polygon([(x_at(L,y0),y0),(x_at(R,y0),y0),(x_at(R,y1),y1),(x_at(L,y1),y1)])
def make(poly):
    out={}
    r=_dot(poly('i'))
    lcgeom.REG.clear(); P=lcgeom.ellipse(r,r,r,r)
    out['period']=_segs(P,[(r,r,r,r)])
    gy=poly('y')
    tail=gy&box(-1000,-1000,3000,0.0)
    tail=max(getattr(tail,'geoms',[tail]),key=lambda p:p.area)
    tail=_strip(tail)
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
def slash(T,f,gap):
    """the slash of the 1958 Manual for Signing and Pavement Marking of the National System of Interstate and Defense
    Highways (figures 3 and 5, "1/2"; there in Series E(M) and E): a parallelogram with flat ends, from 0.225 cap
    heights below the baseline to 1.275 above (1.5 cap heights, centred on the capitals; measured -0.21..1.27 and
    -0.24..1.29), its sides straight and 30 degrees from vertical (measured 29.6 and 30.5), as thick as the capitals'
    stroke T across them (measured 1.06 and 1.14 of the I's stem, with the print's spread).  In other series its
    horizontal run scales by the capitals' width factor f, as their other diagonals do.  Side bearings: gap (the
    drawings' nearest gaps to the figures are 1.1-1.8 strokes) between it and a straight-sided neighbour, nearest at
    the baseline on the left and at the cap line on the right.
    Returns (corner points, advance width)."""
    yb,yt=-0.225*CAP,1.275*CAP; H=yt-yb
    run=np.tan(np.radians(30))*H*f; w=T*np.hypot(run,H)/H          # horizontal width for thickness T across
    x0=gap-run*(0-yb)/H                                            # left edge at the baseline is gap from the origin
    pts=[(x0,yb),(x0+w,yb),(x0+run+w,yt),(x0+run,yt)]
    adv=x0+w+run*(CAP-yb)/H+gap                                    # right edge at the cap line is gap from the advance
    return pts,int(round(adv))
