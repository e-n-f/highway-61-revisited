"""Series E(M) lower case, 1966 sheets (pp. 33-35), as constructions of straight lines and quarter-ellipses.

Units: font units, baseline 0, x-height X = 6 squares, ascender 8 squares, descender 3 squares; one square = 86.25
units (1 inch of an 8-inch letter).  Each letter's construction takes a dict of dimensions in squares; DIM holds the
adopted values (snapped to 1/32 square where the scan allows).  lcfitall.py fits them to the scans."""
import numpy as np
from shapely.geometry import box, Polygon
from lcgeom import Quad, ellipse, half, slope_dir
SQ=690/8; X=6*SQ
def U(p): return {k:(v*SQ if isinstance(v,(int,float)) else v) for k,v in p.items()}
def quad(p,pre):
    g=lambda k: p[pre+k]
    return Quad(g('L'),g('yl'),g('R'),g('yr'),g('xt'),g('T'),g('xb'),g('B'))
# ------------------------------------------------------------------------------------------------------------------
def letter_o(p):
    """outer and inner ellipses (quad ovals), concentric"""
    return quad(p,'O').poly().difference(quad(p,'I').poly())
def letter_n(p):
    """stems w wide at 0 and x2; shoulder: right quarter of an ellipse tangent to the right stem's outer edge, left
    quarter of a flatter one, sharing the top point (xt, T); counter: the upper half of an ellipse (ccx, ccy, ca, cb) wider than the
    counter, cut by the stems' inner edges"""
    w,x2=p['w'],p['x2']; xr=x2+w
    O=Quad(p['oL'],p['oyl'],xr,p['oyr'],p['oxt'],p['T'],p['oxt'],p['oyl']-1).poly()&box(-50,min(p['oyl'],p['oyr']),xr,800)
    stems=box(0,0,w,X)|box(x2,0,xr,p['oyr'])
    cnt=(box(w,-10,x2,p['ccy'])|ellipse(p['ccx'],p['ccy'],p['ca'],p['cb']))&box(w,-10,x2,800)
    body=(stems|(O&box(w,-10,xr,800))).difference(cnt)
    return body
def side(pt,th,ref):
    """half-plane bounded by the line through pt at angle th (degrees), on the side of the point ref"""
    th=np.radians(th); u=np.array([np.cos(th),np.sin(th)]); n=np.array([-u[1],u[0]])
    k=1 if (np.asarray(ref)-np.asarray(pt))@n>0 else -1
    return half(pt,u,k)
def bowl(p,pre='B'):
    """bowl joined to a stem on its right: the union of two quad ovals sharing the left extreme and the top and
    bottom points, so that the upper-right and lower-right quarters (which run into the stem) are independent"""
    g=lambda k: p[pre+k]
    return (Quad(g('L'),g('yl'),g('R'),g('yr'),g('xt'),g('T'),g('xb'),g('B')).poly()|
            Quad(g('L'),g('yl'),g('R2'),g('yr2'),g('xt'),g('T'),g('xb'),g('B')).poly())
def term(p,ref,pre='t'):
    """terminal cut: the line through (tx, ty) at angle tth (degrees), keeping the side of ref"""
    return side((p[pre+'x'],p[pre+'y']),p[pre+'th'],ref)
def letter_a(p):
    """stem xs..xr; arch: quarter-ellipses sharing the top point (axt, aT), the right one tangent to the stem's outer
    edge at height ayr, the left one reaching aL at height ayl, cut by the terminal line; underside: circle (ux, uy, ur);
    bowl: ellipse B; counter: ellipse C; gap and counter end at the stem's inner edge"""
    xs,xr=p['xs'],p['xr']
    A=Quad(p['aL'],p['ayl'],xr,p['ayr'],p['axt'],p['aT'],p['axt'],min(p['ayl'],p['ayr'])-1).poly()
    A=A&box(-100,min(p['ayl'],p['ayr']),xr,800)
    Uc=ellipse(p['ux'],p['uy'],p['ur'],p['ur'])
    B=bowl(p); C=ellipse(p['cx'],p['cy'],p['ca'],p['cb'])
    arch=A.difference(Uc)&term(p,(xr,X))
    stem=box(xs,0,xr,p['ayr']+1)|(A&box(xs,0,xr,800))
    left=box(-100,-100,xs,800)
    return (arch|stem|B).difference(Uc.difference(B)&left).difference(C&left)
def letter_e(p):
    """outer quad oval O cut by the vertical xc (the bar's end); bar from yb to yt; upper counter: ellipse UC above
    the bar; lower counter: ellipse LC with, to its right, everything above the terminal line"""
    O=quad(p,'O').poly()&box(-100,-100,p['xc'],800)
    UC=ellipse(p['ux'],p['uy'],p['ua'],p['ub'])&box(-100,p['yt'],800,800)
    H=term(p,(p['xc']+200,p['yb']))
    LC=(ellipse(p['lx'],p['ly'],p['la'],p['lb'])|H)&box(-100,-100,800,p['yb'])
    return O.difference(UC).difference(LC)
def letter_g(p):
    """stem xs..xr up to the x-height; bowl: ellipse B, counter: ellipse C cut by the stem's inner edge; tail: outer
    ellipse TO (right extreme on the stem's outer edge) and inner ellipse TI, below their centres, cut by the terminal
    line; the bowl is a quad oval"""
    xs,xr=p['xs'],p['xr']
    TO=ellipse(xr-p['toa'],p['toy'],p['toa'],p['tob']); TI=ellipse(p['tix'],p['tiy'],p['tia'],p['tib'])
    outer=box(-100,p['toy'],xr,X)|(TO&box(-100,-400,xr,p['toy']))
    inner=box(-100,p['tiy'],xs,800)|(TI&box(-100,-400,800,p['tiy']))
    tail=outer.difference(inner)&term(p,(xr,0))
    B=bowl(p); C=ellipse(p['cx'],p['cy'],p['ca'],p['cb'])
    return (tail|B).difference(C&box(-100,-100,xs,800))
LETTERS={'o':letter_o,'n':letter_n,'a':letter_a,'e':letter_e,'g':letter_g}
def dims():
    """adopted dimensions (squares): the scan fits in lc1966/fit.json snapped to 1/32 square, angles to 0.5 degree"""
    import json, os
    fn=os.path.join(os.path.dirname(os.path.abspath(__file__)),'lc1966','fit.json')
    out={}
    for ch,r in json.load(open(fn)).items():
        out[ch]={k:(round(v*2)/2 if k.endswith('th') else round(v*32)/32) for k,v in r.items()}
    return out
def series():
    """{letter: (geometry, quarter-ellipses used)} in font units"""
    import lcgeom
    res={}
    for ch,d in dims().items():
        lcgeom.REG.clear()
        g=LETTERS[ch]({k:(v if k.endswith('th') else v*SQ) for k,v in d.items()})
        res[ch]=(g,list(lcgeom.REG))
    return res
