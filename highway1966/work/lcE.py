"""Series E(M) lower case, 1966 sheets (pp. 33-35), as constructions of straight lines and quarter-ellipses.

Units: font units, baseline 0, x-height X = 6 squares, ascender 8 squares, descender 3 squares; one square = 86.25
units (1 inch of an 8-inch letter).  Each letter's construction takes a dict of dimensions in squares; DIM holds the
adopted values (snapped to 1/32 square where the scan allows).  lcfitall.py fits them to the scans."""
import numpy as np
from shapely.geometry import box, Polygon
from lcgeom import Quad, ellipse, half, slope_dir
SQ=690/8; X=6*SQ
NOCUT=False    # True: build letters without their terminal cuts and slanted stem ends (for lccondense)
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
    quarter of a flatter one, sharing the top point (xt, T); counter: quarter-ellipses tangent to the stems' inner edges at cyl, cyr,
    sharing the top point (cxt, cT)"""
    w,x2=p['w'],p['x2']; xr=x2+w
    O=Quad(p['oL'],p['oyl'],xr,p['oyr'],p['oxt'],p['T'],p['oxt'],p['oyl']-1).poly()&box(-50,min(p['oyl'],p['oyr']),xr,800)
    stems=box(0,0,w,X)|box(x2,0,xr,p['oyr'])|box(w,0,x2,max(p['oyl'],p['oyr'])+0.01)
    cnt=counter_arch(w,x2,p['cyl'],p['cyr'],p['cxt'],p['cT'])
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
    if NOCUT: return box(-3000,-3000,3000,3000)
    return side((p[pre+'x'],p[pre+'y']),p[pre+'th'],ref)
def letter_a(p):
    """stem xs..xr; arch: quarter-ellipses sharing the top point (axt, aT), the right one tangent to the stem's outer
    edge at height ayr, the left one reaching aL at height ayl, cut by the terminal line; underside: quad oval U (top point (uxt, uT), left extreme uL at uyl,
    tangent to the stem's inner edge at uyr);
    bowl: ellipse B; counter: ellipse C; gap and counter end at the stem's inner edge"""
    xs,xr=p['xs'],p['xr']
    A=Quad(p['aL'],p['ayl'],xr,p['ayr'],p['axt'],p['aT'],p['axt'],min(p['ayl'],p['ayr'])-1).poly()
    A=A&box(-100,min(p['ayl'],p['ayr']),xr,800)
    Uc=Quad(p['uL'],p['uyl'],xs,p['uyr'],p['uxt'],p['uT'],p['uxt'],2*p['uyl']-p['uT']).poly()
    B=bowl(p); C=ellipse(p['cx'],p['cy'],p['ca'],p['cb'])
    arch=A.difference(Uc)&term(p,(xr,X))
    stem=box(xs,0,xr,p['ayr']+1)|(A&box(xs,0,xr,800))
    left=box(-100,-100,xs,800)
    return (arch|stem|B).difference(Uc.difference(B)&left).difference(C&left)
def letter_e(p):
    """outer quad oval O, whose right side curves down to the bar's lower edge (no cut); bar from yb to yt; upper counter: ellipse UC above
    the bar; lower counter: ellipse LC with, to its right, everything above the terminal line"""
    O=quad(p,'O').poly()
    UC=ellipse(p['ux'],p['uy'],p['ua'],p['ub'])&box(-100,p['yt'],800,800)
    H=Polygon() if NOCUT else term(p,(p['xc']+200,p['yb']))
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
def slanted(x0,x1,y0,ytop,s):
    """stem from x0 to x1, bottom y0, top cut by the line through (x1, ytop) of slope s"""
    if NOCUT: return box(x0,y0,x1,max(ytop,ytop-s*(x1-x0)))
    return Polygon([(x0,y0),(x1,y0),(x1,ytop),(x0,ytop-s*(x1-x0))])
def slanted_bottom(x0,x1,ytop,ybot,s):
    """stem from x0 to x1, top ytop, bottom cut by the line through (x0, ybot) of slope s"""
    if NOCUT: return box(x0,min(ybot,ybot+s*(x1-x0)),x1,ytop)
    return Polygon([(x0,ybot),(x1,ybot+s*(x1-x0)),(x1,ytop),(x0,ytop)])
def xform(g,start,mx=None,rot=None,my=None):
    """mirror about x = mx, or rotate 180 degrees about rot, both the geometry and the ellipses recorded since start"""
    import lcgeom
    from shapely import affinity
    R=lcgeom.REG
    for i in range(start,len(R)):
        cx,cy,a,b=R[i]
        R[i]=(2*mx-cx,cy,a,b) if mx is not None else (cx,2*my-cy,a,b) if my is not None else (2*rot[0]-cx,2*rot[1]-cy,a,b)
    if mx is not None: return affinity.scale(g,-1,1,origin=(mx,0))
    if my is not None: return affinity.scale(g,1,-1,origin=(0,my))
    return affinity.rotate(g,180,origin=rot)
def arch(p,xl,w,x2,pre='',clip=None):
    """an n-arch between a stem ending at xl+w and the stem x2..x2+w: shoulder quarter-ellipses sharing the top point
    (oxt, T), the right one tangent to the right stem's outer edge at oyr, the left reaching oL at oyl; counter: an
    arch of quarter-ellipses tangent to the stems' inner edges at heights cyl and cyr, sharing the top point (cxt, cT).
    Returns (ink above the stems, counter)"""
    g=lambda k: p[pre+k]; xr=x2+w
    O=Quad(xl+g('oL'),g('oyl'),xr,g('oyr'),xl+g('oxt'),g('T'),xl+g('oxt'),min(g('oyl'),g('oyr'))-1).poly()
    O=O&box(xl+w if clip is None else clip,min(g('oyl'),g('oyr')),xr,800)
    cnt=counter_arch(xl+w,x2,g('cyl'),g('cyr'),xl+g('cxt'),g('cT'))
    return fuse(O|box(x2,0,xr,g('oyr'))|box(xl+w,0,x2,max(g('oyl'),g('oyr'))+0.01)),cnt
def counter_arch(x0,x2,yl,yr,xt,T):
    """a counter between vertical edges x0 and x2, closed at the top by quarter-ellipses tangent to the edges at
    heights yl and yr and sharing the top point (xt, T)"""
    lo=min(yl,yr)
    Q=Quad(x0,yl,x2,yr,xt,T,xt,lo-1).poly()&box(x0,lo,x2,800)
    return fuse(box(x0,-10,x2,lo)|Q|box(x0,lo-1,xt,yl)|box(xt,lo-1,x2,yr))
def letter_l(p):
    return slanted(0,p['w'],0,p['ytop'],p['s'])
def letter_i(p):
    return box(p['x0'],0,p['x0']+p['w'],X)|ellipse(p['dx'],p['dy'],p['dr'],p['dr'])
def letter_h(p):
    ink,cnt=arch(p,0,p['w'],p['x2'])
    return (slanted(0,p['w'],0,p['ytop'],p['s'])|ink).difference(cnt)
def letter_m(p):
    w=p['w']; ink1,c1=arch(p,0,w,p['x2'],'')
    ink2,c2=arch(p,p['x2'],w,p['x3'],'r',clip=p['x2'])
    return fuse(box(0,0,w,X)|ink1|ink2).difference(c1).difference(c2)
def letter_u(p):
    """the n turned through 180 degrees"""
    import lcgeom
    k=len(lcgeom.REG); g=letter_n(p)
    return xform(g,k,rot=((p['x2']+p['w'])/2,X/2))
def letter_r(p):
    """stem and an n shoulder, cut by the terminal line"""
    ink,cnt=arch(p,0,p['w'],p['x2'])
    return (box(0,0,p['w'],X)|((ink.difference(cnt))&term(p,(0,X)))).difference(box(-100,-100,800,0))
def bowl_left(p,pre='B'):
    """bowl joined to a stem on its left: two quad ovals sharing the right extreme and the top and bottom points"""
    g=lambda k: p[pre+k]
    return (Quad(g('L'),g('yl'),g('R'),g('yr'),g('xt'),g('T'),g('xb'),g('B')).poly()|
            Quad(g('L2'),g('yl2'),g('R'),g('yr'),g('xt'),g('T'),g('xb'),g('B')).poly())
def letter_d(p):
    C=ellipse(p['cx'],p['cy'],p['ca'],p['cb'])&box(-100,-100,p['xs'],800)
    return (slanted(p['xs'],p['xr'],0,p['ytop'],p['s'])|bowl(p)).difference(C)
def letter_q(p):
    C=ellipse(p['cx'],p['cy'],p['ca'],p['cb'])&box(-100,-100,p['xs'],800)
    return (slanted_bottom(p['xs'],p['xr'],X,p['ybot'],p['sb'])|bowl(p)).difference(C)
def letter_b(p):
    C=ellipse(p['cx'],p['cy'],p['ca'],p['cb'])&box(p['w'],-100,800,800)
    return (slanted(0,p['w'],0,p['ytop'],p['s'])|bowl_left(p)).difference(C)
def letter_p(p):
    C=ellipse(p['cx'],p['cy'],p['ca'],p['cb'])&box(p['w'],-100,800,800)
    return (slanted_bottom(0,p['w'],X,p['ybot'],p['sb'])|bowl_left(p)).difference(C)
def letter_c(p):
    """outer quad oval O, counter ellipse I; the opening is the wedge to the right of both terminal lines (upper
    through (ux, uy) at uth, lower through (lx, ly) at lth)"""
    O=quad(p,'O').poly(); I=ellipse(p['ix'],p['iy'],p['ia'],p['ib'])
    far=(p['OR']+300,p['Oyr'])
    wedge=side((p['ux'],p['uy']),p['uth'],far)&side((p['lx'],p['ly']),p['lth'],far)
    return O.difference(I) if NOCUT else O.difference(I).difference(wedge)
def hook(p,top):
    """a stem x0..x1 from `top` down into a hook turning right: quarter-ellipses sharing the bottom point (hx, hB), the
    left one rounding the stem's outer edge from height hyl, the right one rising to hR at hyr; inside, an ellipse
    (ia, ib) centred at height icy, tangent to the stem's inner edge, below the open side; cut by the terminal line"""
    x0,x1=p['x0'],p['x1']
    BLq=ellipse(p['hx'],p['hyl'],p['hx']-x0,p['hyl']-p['hB'])&box(x0,-1000,p['hx'],p['hyl'])
    if abs(p['tth']-90)<30:                              # a steep terminal cut: the outer curve runs out past it,
        p=dict(p,hR=max(p['hR'],line_x(p['tx'],p['ty'],p['tth'],p['hyr'])+0.75*SQ))   # so the cut ends in a corner
    BRq=ellipse(p['hx'],p['hyr'],p['hR']-p['hx'],p['hyr']-p['hB'])&box(p['hx'],-1000,p['hR'],p['hyr'])
    ink=fuse(box(x0,p['hyl'],x1,top)|BLq|BRq)
    # between the stem's inner edge and the bottom point, down to where the quarters begin, is solid
    ink=fuse(ink|box(min(x1,p['hx']),min(p['hyl'],p['hyr']),max(x1,p['hx']),top))
    # below the right quarter's vertical tangent the outer edge runs straight on, down to the inner ellipse's centre
    if p['icy']>p['hyr']: ink=fuse(ink|box(p['hx'],p['hyr'],p['hR'],p['icy']))
    # and so is the strip beside the inner ellipse, from the hook body up to the ellipse's centre
    if p['icy']>min(p['hyl'],p['hyr']): ink=fuse(ink|box(x1,min(p['hyl'],p['hyr']),x1+p['ia'],p['icy']))
    inner=ellipse(x1+p['ia'],p['icy'],p['ia'],p['ib'])|box(x1,p['icy'],3000,3000)
    ink=ink.difference(inner)
    # a shallow cut (t) applies only beyond the inner ellipse's centre, so it cannot reach the counter's stem side;
    # a steep one (f, j) applies everywhere right of the stem
    right=box(x1+(0 if abs(p['tth']-90)<30 else p['ia']),-3000,3000,3000)
    return fuse(ink.difference(right)|(ink&right&term(p,(x1,p['hB']+1))))
def letter_t(p):
    return fuse(hook(p,p['by'])|slanted(p['x0'],p['x1'],p['hyl'],p['ytop'],p['s'])|box(p['bx0'],p['by'],p['bx1'],X))
def letter_j(p):
    """the hook turned to the left (mirrored about the stem's axis), stem up to the x-height, and the dot"""
    import lcgeom
    k=len(lcgeom.REG); g=hook(p,X); g=xform(g,k,mx=(p['x0']+p['x1'])/2)
    return g|ellipse(p['dx'],p['dy'],p['dr'],p['dr'])
def letter_f(p):
    """the hook turned upside down (about the line y = X/2... mirrored about y = 3 squares), with the bar"""
    import lcgeom
    k=len(lcgeom.REG); my=X/2
    g=hook(p,2*my-0.0)
    g=xform(g,k,my=my)
    return fuse(g|box(p['bx0'],p['by'],p['bx1'],X))
def line_x(x0,y0,th,y):
    t=np.radians(th); return x0+(y-y0)*np.cos(t)/np.sin(t)
def rest(q,m,side):
    """vertical shift for quad oval m so that it rests on quad oval q from above (side=1; m's lower quarters touch q)
    or hangs under it (side=-1; m's upper quarters touch q), touching from outside.  Returns (shift, x of the contact)"""
    Q=m.quarters()
    E=np.vstack([quarter(*Q[k],k) for k in (('BL','BR') if side>0 else ('TL','TR'))])
    if side>0: lo,hi=min(q.yl,q.yr)-m.B,q.T-m.B+10
    else: lo,hi=q.B-m.T-10,max(q.yl,q.yr)-m.T
    for _ in range(60):
        c=(lo+hi)/2; clear=q.dist(E+[0,c]).min()>0
        if (side>0)==clear: hi=c
        else: lo=c
    c=hi if side>0 else lo
    d=q.dist(E+[0,c]); return c,E[int(np.argmin(np.abs(d))),0]
def under(q,x0,side):
    """side=1: everything below q's upper edge to the right of x0 (out to q's right extreme); side=-1: everything above
    q's lower edge to the left of x0 (out to q's left extreme)"""
    Q=q.quarters()
    if side>0:
        E=np.vstack([quarter(*Q['TL'],'TL')[::-1],quarter(*Q['TR'],'TR')[::-1]])   # left extreme -> top -> right extreme
        E=E[E[:,0]>=x0]
        return Polygon(np.vstack([[x0,-3000],E,[q.R,-3000]])).buffer(0)
    E=np.vstack([quarter(*Q['BL'],'BL'),quarter(*Q['BR'],'BR')])               # left extreme -> bottom -> right extreme
    E=E[E[:,0]<=x0]
    return Polygon(np.vstack([[q.L,3000],E,[x0,3000]])).buffer(0)
def shifted(g,pre,dy):
    return Quad(g[pre+'L'],g[pre+'yl']+dy,g[pre+'R'],g[pre+'yr']+dy,g[pre+'xt'],g[pre+'T']+dy,g[pre+'xb'],g[pre+'B']+dy)
def letter_s(p):
    """two outer quad ovals: the upper (U) gives the top arc, the left side and the spine's lower edge; the lower (L)
    the bottom arc, the right side and the spine's upper edge.  The counters are quad ovals UI (resting on L's upper
    edge) and LI (hanging under U's lower edge), each touching the spine edge tangentially; past the touching point
    the counter is bounded by the spine edge, and it opens past its terminal line"""
    from lcgeom import quarter
    globals()['quarter']=quarter
    p=dict(p)
    # the outer ovals run out past the terminal lines, so each cut ends in a sharp corner
    p['LL']=min(p['LL'],line_x(p['ltx'],p['lty'],p['ltth'],p['Lyl'])-2.0)
    p['UR']=max(p['UR'],line_x(p['utx'],p['uty'],p['utth'],p['Uyr'])+2.0)
    UO=quad(p,'U'); LO=quad(p,'L')
    du,xu=rest(LO,quad(p,'ui'),1); dl,xl=rest(UO,quad(p,'li'),-1)
    for _ in range(2):                                   # the counters run out past their terminal lines
        p['uiR']=max(p['uiR'],line_x(p['utx'],p['uty'],p['utth'],p['uiyr']+du)+2.0)
        p['liL']=min(p['liL'],line_x(p['ltx'],p['lty'],p['ltth'],p['liyl']+dl)-2.0)
        du,xu=rest(LO,quad(p,'ui'),1); dl,xl=rest(UO,quad(p,'li'),-1)
    UI=shifted(p,'ui',du).poly(); LI=shifted(p,'li',dl).poly()
    yu=p['uiyl']+du; yl=p['liyr']+dl
    UO=UO.poly(); LO=LO.poly()
    Hu=side((p['utx'],p['uty']),p['utth'],(p['UR']+300,p['UT'])); Hl=side((p['ltx'],p['lty']),p['ltth'],(p['LL']-300,p['LB']))
    if NOCUT: Hu=Hl=Polygon()
    Wu=(UI|Hu|under(shifted(p,'ui',du),xu,1)).difference(LO)      # the upper counter, open to the right
    Wl=(LI|Hl|under(shifted(p,'li',dl),xl,-1)).difference(UO)      # the lower counter, open to the left
    return fuse(UO|LO).difference(Wu).difference(Wl)
LETTERS={'o':letter_o,'n':letter_n,'a':letter_a,'e':letter_e,'g':letter_g,'l':letter_l,'i':letter_i,'h':letter_h,
         'm':letter_m,'u':letter_u,'r':letter_r,
         'b':letter_b,'d':letter_d,'p':letter_p,'q':letter_q,'c':letter_c,'t':letter_t,'j':letter_j,'f':letter_f,'s':letter_s}
def fuse(g,e=1e-3):
    """weld pieces that only touch along an edge (float slivers between unions)"""
    return g.buffer(e,join_style=2).buffer(-e,join_style=2)
def clean(g,e=0.01):
    """weld touching pieces and drop slivers (holes or pieces under one square unit)"""
    g=fuse(g,e)
    ps=[]
    for p in getattr(g,'geoms',[g]):
        if p.area<1: continue
        ps.append(Polygon(p.exterior,[i for i in p.interiors if Polygon(i).area>=1]))
    from shapely.geometry import MultiPolygon
    return ps[0] if len(ps)==1 else MultiPolygon(ps)
def polys():
    import json, os
    fn=os.path.join(os.path.dirname(os.path.abspath(__file__)),'lc1966','poly.json')
    return json.load(open(fn)) if os.path.exists(fn) else {}
def dims():
    """adopted dimensions (squares): the scan fits in lc1966/fit.json snapped to 1/32 square, angles to 0.5 degree"""
    import json, os
    fn=os.path.join(os.path.dirname(os.path.abspath(__file__)),'lc1966','fit.json')
    out={}
    for ch,r in json.load(open(fn)).items():
        out[ch]={k:(round(v*2)/2 if k.endswith('th') else round(v*1000)/1000 if k in ('s','sb') else round(v*32)/32) for k,v in r.items()}
    return out
def build(ch,nocut=False):
    """(geometry, quarter-ellipses) of one constructed letter, optionally without its cuts"""
    global NOCUT
    import lcgeom
    d=dims()[ch]; lcgeom.REG.clear(); NOCUT=nocut
    try: g=clean(LETTERS[ch]({k:(v if (k.endswith('th') or k in ('s','sb')) else v*SQ) for k,v in d.items()}))
    finally: NOCUT=False
    return g,list(lcgeom.REG)
def series():
    """{letter: (geometry, quarter-ellipses used)} in font units"""
    import lcgeom
    res={}
    for ch,d in dims().items():
        lcgeom.REG.clear()
        g=clean(LETTERS[ch]({k:(v if (k.endswith('th') or k in ('s','sb')) else v*SQ) for k,v in d.items()}))
        res[ch]=(g,list(lcgeom.REG))
    for ch,V in polys().items():
        res[ch]=(Polygon(np.array(V)*SQ),[])
    return res
