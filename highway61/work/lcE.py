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
def letter_s_quads(p):
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
def _ey(cx,cy,a,b,x,sg):
    """y of the ellipse's upper (sg=1) or lower (sg=-1) arc at x"""
    t=min(abs(x-cx)/a,1.0); return cy+sg*b*np.sqrt(1-t*t)
def _cut_meets(p,pre,e,hi=True):
    """where the terminal line crosses ellipse e: the upper crossing (hi) or the lower"""
    cx,cy,a,b=e; th=np.radians(p[pre+'th']); u=np.array([np.cos(th),np.sin(th)]); o=np.array([p[pre+'x'],p[pre+'y']])
    d=o-np.array([cx,cy]); A=u[0]**2/a**2+u[1]**2/b**2; B=2*(d[0]*u[0]/a**2+d[1]*u[1]/b**2); C=d[0]**2/a**2+d[1]**2/b**2-1
    disc=B*B-4*A*C
    if disc<0: return o
    pts=[o+u*(-B-np.sqrt(disc))/(2*A),o+u*(-B+np.sqrt(disc))/(2*A)]
    return max(pts,key=lambda q:q[1]) if hi else min(pts,key=lambda q:q[1])
def _F(e,P):
    cx,cy,a,b=e; return ((P[:,0]-cx)/a)**2+((P[:,1]-cy)/b)**2-1
def _touching(cx,a,Fin,Fout,y0,b0):
    """an ellipse centred at x = cx, half-width a, inside the region Fin <= 0 and outside the region Fout < 0,
    touching both (exactly: the extreme over the ellipse is refined continuously): its centre y, half-height and the
    two touching points"""
    from scipy.optimize import fsolve, minimize_scalar
    t=np.linspace(0,2*np.pi,721)[:-1]; h=t[1]-t[0]
    def pt(v,th): return np.array([[cx+a*np.cos(th),v[0]+abs(v[1])*np.sin(th)]])
    def ext(v,F,sg):                               # sg=1: max of F over the ellipse, -1: min; and where
        P=np.c_[cx+a*np.cos(t),v[0]+abs(v[1])*np.sin(t)]; k=np.argmax(sg*F(P))
        r=minimize_scalar(lambda th:-sg*F(pt(v,th))[0],bounds=(t[k]-h,t[k]+h),method='bounded',options={'xatol':1e-12})
        return sg*-r.fun,r.x
    def eqs(v): return [ext(v,Fin,1)[0],ext(v,Fout,-1)[0]]
    v=fsolve(eqs,[y0,b0],xtol=1e-12)
    return v[0],abs(v[1]),pt(v,ext(v,Fin,1)[1])[0],pt(v,ext(v,Fout,-1)[1])[0]
def _tangent(e,x,sg):
    """point on the upper (sg=1) / lower (sg=-1) arc of e at x, and the unit tangent there (pointing to +x)"""
    cx,cy,a,b=e; y=_ey(cx,cy,a,b,x,sg); dydx=-(b*b*(x-cx))/(a*a*(y-cy)) if abs(y-cy)>1e-9 else 0.0
    u=np.array([1.0,dydx]); return np.array([x,y]),u/np.linalg.norm(u)
def s_ellipses(p):
    """the six ellipses of the s, (cx, cy, a, b): U1 U2 UM through the top point, L1 L2 LM through the bottom point"""
    xt,T,xb,B=p['xt'],p['T'],p['xb'],p['B']
    up=lambda k: (xt,T-p[k+'b'],p[k+'a'],p[k+'b']); lo=lambda k: (xb,B+p[k+'b'],p[k+'a'],p[k+'b'])
    return up('U1'),up('U2'),up('UM'),lo('L1'),lo('L2'),lo('LM')
def letter_s_ell(p):
    """the s as ellipses (the user's reading of the 1966 drawing).  Upper family, all through the top point (xt, T):
    the big ellipse U1 (the top stroke's upper edge, to the upper terminal line), U2 (the upper left side and, with
    its bottom arc, the spine's lower edge) and UM (its bottom arc the bottom stroke's upper edge).  Lower family,
    all through the bottom point (xb, B): L1 (the bottom stroke's lower edge, from the lower terminal line), L2 (the
    lower right side and, with its top arc, the spine's upper edge) and LM (its top arc the top stroke's underside).
    The spine's middle is straight: its upper edge leaves L2's top arc at x = s1 along the tangent, its lower edge
    leaves U2's bottom arc at x = s2 along the tangent.  The counters are the bands inside LM above the spine's upper
    edge and inside UM below its lower edge, each ended by a small ellipse (C1, C2: centre x and half-width given)
    touching both of the band's edges."""
    if 'ell' in p: u1,u2,um,l1,l2,lm=p['ell']           # ellipses given directly (a condensed s: lccondense)
    else: u1,u2,um,l1,l2,lm=s_ellipses(p)
    xt,xb=u1[0],l1[0]
    U1,U2,UM,L1,L2,LM=[ellipse(*e) for e in (u1,u2,um,l1,l2,lm)]
    q1,v1=_tangent(l2,p['s1'],1); q2,v2=_tangent(u2,p['s2'],-1)
    above1=lambda P: (P[:,1]-q1[1])*v1[0]-(P[:,0]-q1[0])*v1[1]          # > 0 above the spine's upper edge line
    below2=lambda P: -((P[:,1]-q2[1])*v2[0]-(P[:,0]-q2[0])*v2[1])       # > 0 below the spine's lower edge line
    y1,b1,ta,tb=_touching(p['C1x'],p['C1a'],lambda P:_F(lm,P),above1,_ey(*lm,p['C1x'],1)-0.6*p['C1a'],0.4*p['C1a'])
    y2,b2,sa,sb=_touching(p['C2x'],p['C2a'],lambda P:_F(um,P),below2,_ey(*um,p['C2x'],-1)+0.6*p['C2a'],0.4*p['C2a'])
    C1=ellipse(p['C1x'],y1,p['C1a'],b1); C2=ellipse(p['C2x'],y2,p['C2a'],b2)
    def beyond(q,u,keep):                          # the half-plane past the line through q along u, on keep's side
        return side(q,np.degrees(np.arctan2(u[1],u[0])),keep)
    S1=beyond(q1,v1,(q1[0],q1[1]-100))&box(-3000,-3000,q1[0],3000)     # below the straight part of the upper edge
    S2=beyond(q2,v2,(q2[0],q2[1]+100))&box(q2[0],-3000,3000,3000)      # above the straight part of the lower edge
    Wu=C1|(LM.difference(L2).difference(S1)&beyond(ta,tb-ta,(p['C1x']+1000,y1)))
    Wl=C2|(UM.difference(U2).difference(S2)&beyond(sa,sb-sa,(p['C2x']-1000,y2)))
    yu=_cut_meets(p,'ut',lm,True)[1] if not NOCUT else u1[1]; yl=_cut_meets(p,'lt',um,False)[1] if not NOCUT else l1[1]
    Tu=U1&box(xt,max(u1[1],yu),3000,3000)&term(p,p.get('utref',(xt,u1[1])),'ut')
    Tl=L1&box(-3000,-3000,xb,min(l1[1],yl))&term(p,p.get('ltref',(xb,l1[1])),'lt')
    if not NOCUT:                                   # U2 and L2 also stop at the terminal lines (outside LM, UM)
        U2=U2.difference(box(-3000,-3000,3000,3000).difference(term(p,p.get('utref',(xt,u1[1])),'ut')).difference(LM))
        L2=L2.difference(box(-3000,-3000,3000,3000).difference(term(p,p.get('ltref',(xb,l1[1])),'lt')).difference(UM))
    spine=(S1&box(-3000,-3000,3000,3000)).intersection(S2)&box(min(q1[0],q2[0]),-3000,max(q1[0],q2[0]),3000)
    return fuse((U2|Tu|spine).difference(Wu)|(L2|Tl|spine).difference(Wl))
LETTERS={'o':letter_o,'n':letter_n,'a':letter_a,'e':letter_e,'g':letter_g,'l':letter_l,'i':letter_i,'h':letter_h,
         'm':letter_m,'u':letter_u,'r':letter_r,
         'b':letter_b,'d':letter_d,'p':letter_p,'q':letter_q,'c':letter_c,'t':letter_t,'j':letter_j,'f':letter_f,'s':letter_s_ell}
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
    for ch in 'bdpq':
        lighten(out[ch],ch in 'dq')
    if COND>0: condensed_ae(out,COND)
    if COMP>0 or FHOOK: outer_bands(out,COMP)
    return out
COMP=0.0       # extra vertical erosion (units) of the series being built in the ascender and descender bands; lccondense
def outer_bands(d,delta):
    """the horizontal strokes above the x-height and below the baseline (the f's top, the g's and j's tails) made
    2 delta thicker, by shrinking their counter ellipses' vertical semi-axes: there the condensation only shifts the
    letter, and its erosion (which, with the stretch, makes the x-height band's strokes the series weight) would
    leave them 2 delta too light"""
    t=2*delta/SQ
    d['f']['ib']-=t+FHOOK/SQ; d['j']['ib']-=t; d['g']['tib']-=t
FHOOK=0.0      # further thickening (units) of the f's hook, beyond outer_bands' (lccondense.F_HOOK sets it)
COND=0.0       # condensation (1 - lower-case width factor) of the series being built; lccondense sets it
# The Palo Alto signs (Oregon Expressway, early 1960s, Series C lower case at our C's widths) have the e's bar and the
# a's bowl top as heavy as the top and bottom strokes, grown into the e's lower opening and the a's gap below the
# arch, with the bar and bowl top a little further down and up again (measured from stroke centres, which photographic
# bloom does not move).  Applied in proportion to the condensation, full at Series C's, the weight at most the stroke.
C_REF=1-0.82                 # Series C's condensation (lccondense.F_LC['C'])
E_DROP=0.0256                # extra downward move of the e's bar at C_REF (squares; from 0.494 of the centre span)
A_RISE=0.012                 # extra upward move of the a's bowl top at C_REF (squares; from 0.399)
def condensed_ae(d,c):
    """the e's bar and the a's bowl top for condensation c (squares, E(M) construction)"""
    s=c/C_REF; w=min(1.0,s)
    e=d['e']; top=e['OT']-(e['uy']+e['ub'])                   # top stroke at the centre
    t=e['yt']-e['yb']; t2=t+w*(top-t)                          # bar: towards the top stroke's weight
    e['yt']-=s*E_DROP; e['yb']=e['yt']-t2; e['Oyr']=e['yb']   # keeps the eye's bottom (less the move); right side tied
    a=d['a']; top=a['aT']-a['uT']
    cb0=a['cy']-a['cb']; ct=a['cy']+a['cb']; t=a['BT']-ct; t2=t+w*(top-t)
    ct+=s*A_RISE; a['BT']=ct+t2                                # bowl top over the counter (raised by the move)
    a['cy'],a['cb']=(ct+cb0)/2,(ct-cb0)/2
BOWL_T=(1.53125,1.5)   # side and top/bottom thickness of the b d p q bowls (squares), as the o's; the fitted counters
                       # left them 0.03-0.13 square heavier, which read as heavy next to o and g
def lighten(d,stem_right):
    """enlarge the bowl's counter ellipse (centre x kept) so the bowl is BOWL_T thick at its free side, top and bottom"""
    ts,tt=BOWL_T
    d['ca']=round((d['cx']-(d['BL']+ts) if stem_right else (d['BR']-ts)-d['cx'])*32)/32
    d['cy']=round((d['BT']+d['BB'])/2*32)/32
    d['cb']=round(((d['BT']-d['BB'])/2-tt)*32)/32
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
