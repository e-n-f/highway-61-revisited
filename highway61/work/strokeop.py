"""Change the stroke weight of an exact outline while keeping its vertical metrics.

The outline is a list of contours of segments, ('L', p0, p1) or ('E', p0, p1, (cx, cy, a, b), t0, sweep) — an arc
of an axis-aligned ellipse (a circle when a = b) from parameter t0 through `sweep` radians.  Contours run with the
ink on their right (outline.contours / lcoutline.contours orientation).

offset(): every segment moves sideways by d (d > 0 grows the ink): a line along its normal, an arc to the parallel
arc with the same centre and both semi-axes changed by d (exact for circles; for ellipses the draftsman's parallel
ellipse).  Neighbouring segments are re-joined: at a tangent joint the joint moves along the common normal, at a
corner the two moved curves are intersected.

vmap(): a piecewise-linear map of y (each band affine), which keeps lines straight and maps axis-aligned ellipses
to axis-aligned ellipses; segments are split where they cross a band boundary.

For Series E -> E(M) capitals: offset by (138 - A)/2 then squeeze y back to the baseline and cap line.  For the
lower case the opposite: offset inwards, then stretch each band (descender, x-height, ascender) back to its lines."""
import numpy as np
import outline, lcoutline
# ---------------------------------------------------------------------------------------------------- segments
def from_caps(geom):
    """1961 capitals: lines and circular arcs (outline.contours) as segment lists"""
    out=[]
    for segs in outline.contours(geom):
        cs=[]
        for s in segs:
            if s[0]=='L': cs.append(('L',np.array(s[1],float),np.array(s[2],float)))
            else:
                _,p0,p1,c,r,sw=s; t0=np.arctan2(p0[1]-c[1],p0[0]-c[0])
                cs.append(('E',np.array(p0,float),np.array(p1,float),(c[0],c[1],r,r),t0,sw))
        out.append(cs)
    return out
def from_lc(geom,prims):
    """lower case: lines and elliptical arcs (lcoutline.contours)"""
    return [[(s[0],np.array(s[1],float),np.array(s[2],float))+tuple(s[3:]) if s[0]=='E' else ('L',np.array(s[1],float),np.array(s[2],float))
             for s in segs] for segs in lcoutline.contours(geom,prims)]
def _ept(e,t): cx,cy,a,b=e; return np.array([cx+a*np.cos(t),cy+b*np.sin(t)])
def _etan(e,t): cx,cy,a,b=e; return np.array([-a*np.sin(t),b*np.cos(t)])
def _par(e,p):
    cx,cy,a,b=e; return np.arctan2((p[1]-cy)/b,(p[0]-cx)/a)
def tangent(s,end):
    """unit direction of travel at the start (end=0) or end (end=1) of a segment"""
    if s[0]=='L':
        v=s[2]-s[1]
    else:
        t=s[4]+(s[5] if end else 0.0); v=_etan(s[3],t)*np.sign(s[5])
    return v/np.linalg.norm(v)
# ---------------------------------------------------------------------------------------------------- offset
def _support(n,d):
    """the offset vector for outward normal n: d is a number (a circle) or (dx, dy) (an axis-aligned ellipse: the
    Minkowski sum with an ellipse moves a line by its support distance along the normal)"""
    if callable(d): return n*d(n)                            # an amount for each direction
    if np.isscalar(d): return n*d
    dx,dy=d; sg=1.0 if (dx if dx!=0 else dy)>=0 else -1.0; h=np.hypot(n[0]*dx,n[1]*dy)
    if h<1e-12: return np.zeros(2)
    return sg*np.array([dx*dx*n[0],dy*dy*n[1]])/max(h,1e-12)
def _moved(s,d):
    """the moved carrier of a segment: ('L', point, direction) or ('E', ellipse)"""
    if s[0]=='L':
        u=tangent(s,0); n=np.array([-u[1],u[0]])            # left of travel = away from the ink
        return ('L',s[1]+_support(n,d),u)
    cx,cy,a,b=s[3]; k=np.sign(s[5])                         # travelling anticlockwise: the centre is on the left
    if callable(d): dx,dy=d(np.array([1.0,0.0])),d(np.array([0.0,1.0]))   # semi-axes: the amounts at the axis points
    else: dx,dy=(d,d) if np.isscalar(d) else d
    return ('E',(cx,cy,a-k*dx,b-k*dy))
def _on(c,t):
    return c[1]+c[2]*t if c[0]=='L' else _ept(c[1],t)
def _don(c,t):
    return c[2] if c[0]=='L' else _etan(c[1],t)
def _param(c,p):
    return (p-c[1])@c[2] if c[0]=='L' else _par(c[1],p)
def _intersect(c1,c2,guess):
    """intersection of two carriers near `guess` (Newton on the two parameters)"""
    t1,t2=_param(c1,guess),_param(c2,guess)
    for _ in range(60):
        F=_on(c1,t1)-_on(c2,t2)
        if np.hypot(*F)<1e-9: break
        J=np.c_[_don(c1,t1),-_don(c2,t2)]
        try: dt=np.linalg.solve(J,-F)
        except np.linalg.LinAlgError: return None
        lim=0.3 if c1[0]=='E' or c2[0]=='E' else 1e9
        if c1[0]=='E': dt[0]=np.clip(dt[0],-lim,lim)
        if c2[0]=='E': dt[1]=np.clip(dt[1],-lim,lim)
        t1+=dt[0]; t2+=dt[1]
    p=_on(c1,t1)
    return p if np.hypot(*(_on(c1,t1)-_on(c2,t2)))<1e-6 else None
def offset_contour(segs,d,tags=None):
    """offset one contour; segments that the offset makes vanish or run backwards are dropped and their neighbours
    re-joined, until none are left.  With tags (one per segment), returns (segments, surviving tags)"""
    segs=list(segs); tags=list(tags) if tags is not None else None
    for _ in range(len(segs)):
        out=_offset_once(segs,d)
        if len(segs)<=2: break
        idx=_offset_once.idx
        bad=[i for k,i in enumerate(idx) if i>=0 and _reversed(segs[i],out[k])]
        if not bad: break
        segs=[s for i,s in enumerate(segs) if i!=bad[0]]
        if isinstance(d,list): d=[x for i,x in enumerate(d) if i!=bad[0]]
        if tags is not None: tags=[t for i,t in enumerate(tags) if i!=bad[0]]
    if tags is None: return out
    idx=_offset_once.idx; tg=[]                             # one tag per output segment (a bevel takes the next one's)
    for k_,i in enumerate(idx):
        if i>=0: tg.append(tags[i])
        else:
            nxt=next((j for j in idx[k_+1:] if j>=0),idx[0] if idx[0]>=0 else 0); tg.append(tags[nxt])
    return out,tg
def _reversed(s,o):
    if s[0]=='L':
        return (o[2]-o[1])@(s[2]-s[1])<=0
    if abs(o[5])<1e-9: return True
    if np.sign(o[5])!=np.sign(s[5]): return True
    # an arc whose chord now points against the original chord has turned inside out
    c0=s[2]-s[1]; c1=o[2]-o[1]
    return np.hypot(*c0)>1e-6 and abs(s[5])<np.pi and c0@c1<0
MITRE=None      # when set, corners whose mitre would reach further than MITRE times the offset are bevelled
def _offset_once(segs,d):
    n=len(segs)
    D=d if isinstance(d,list) else [d]*n                    # one amount per segment, or the same for all
    C=[_moved(s,D[i]) for i,s in enumerate(segs)]
    J=[]; BEV={}
    for i in range(n):
        a,b=segs[i-1],segs[i]; P=a[2]                         # joint between segment i-1 and i
        ua,ub=tangent(a,1),tangent(b,0)
        na=np.array([-ua[1],ua[0]]); nb=np.array([-ub[1],ub[0]])
        smooth=np.degrees(np.arccos(np.clip(ua@ub,-1,1)))<0.5
        if smooth:
            nm=(na+nb)/np.linalg.norm(na+nb)
            J.append(P+(_support(nm,D[i-1])+_support(nm,D[i]))/2); continue
        # a corner: intersect the moved carriers near the mitre point
        def amt(dv):
            if callable(dv): mm=(na+nb)/max(np.linalg.norm(na+nb),1e-9); return dv(mm)
            return dv if np.isscalar(dv) else max(abs(dv[0]),abs(dv[1]))*(np.sign(dv[0]) if dv[0]!=0 else np.sign(dv[1]))
        dd=(amt(D[i-1])+amt(D[i]))/2
        m=na+nb; g=P+m/max(np.linalg.norm(m),1e-9)*dd/max(np.cos(np.arccos(np.clip(na@nb,-1,1))/2),0.2)
        q=_intersect(C[i-1],C[i],g)
        q=q if q is not None else g
        def own(dv,nn): return abs(dv(nn)) if callable(dv) else abs(amt(dv))
        big=max(abs(dd),own(D[i-1],na),own(D[i],nb))  # (amounts that differ by side: each along its own normal)
        # a long mitre is bevelled only where the corner is acute enough to spike: a long one at a blunter corner
        # comes from its sides moving differently (the a's counter meeting its stem), which only slides the corner
        # along them, or from a crotch (the m's, which a bevel left flat)
        acute=MITRE is not None and np.sqrt(max((1+na@nb)/2,0.0))*MITRE<1
        if acute and np.hypot(*(q-P))>MITRE*big:
            pa=P+_support(na,D[i-1]); pb=P+_support(nb,D[i]); BEV[i]=(pa,pb); J.append(None); continue
        J.append(q)
    out=[]; IDX=[]
    for i in range(n):
        if i in BEV: out.append(('L',BEV[i][0],BEV[i][1])); IDX.append(-1)
        IDX.append(i)
        j=(i+1)%n
        p0=BEV[i][1] if i in BEV else J[i]; p1=BEV[j][0] if j in BEV else J[j]; s=segs[i]
        if s[0]=='L': out.append(('L',p0,p1)); continue
        e=C[i][1]; t0=_par(e,p0); t1=_par(e,p1); k=np.sign(s[5])
        if n==1: out.append(('E',p0,p1,e,t0,s[5])); continue      # a whole ellipse (a dot)
        sw=(t1-t0)%(2*np.pi) if k>0 else -((t0-t1)%(2*np.pi))
        if abs(sw)>abs(s[5])+np.pi: sw=0.0                    # the arc vanished (crossed over itself)
        out.append(('E',p0,p1,e,t0,sw))
    _offset_once.idx=IDX
    return out
def clean(segs,minlen=0.05,tags=None):
    """drop segments that collapsed to (nearly) nothing; with tags, return (segments, tags)"""
    keep=[i for i,s in enumerate(segs) if np.hypot(*(s[2]-s[1]))>minlen or (s[0]=='E' and abs(s[5])>1e-3)]
    out=[segs[i] for i in keep]
    return (out,[tags[i] for i in keep]) if tags is not None else out
def offset(contours,d):
    return [clean(offset_contour(c,d)) for c in contours]
# ---------------------------------------------------------------------------------------------------- vertical map
def vmap(contours,bands):
    """bands: list of (y0, y1, k, c): y in [y0, y1) maps to c + k*y; the bands must be contiguous (y0 of each = y1 of
    the previous); the first and last extend to infinity"""
    def band(y):
        for i,(y0,y1,k,c) in enumerate(bands):
            if (y<y1 or i==len(bands)-1) and (y>=y0 or i==0): return i
    def f(i,y): y0,y1,k,c=bands[i]; return c+k*y
    cuts=[b[1] for b in bands[:-1]]
    out=[]
    for segs in contours:
        cs=[]
        for s in segs:
            pieces=[s]
            for yc in cuts: pieces=[q for p in pieces for q in _split(p,yc)]
            for p in pieces:
                i=band((p[1][1]+p[2][1])/2 if p[0]=='L' else _ept(p[3],p[4]+p[5]/2)[1])
                M=lambda q: np.array([q[0],f(i,q[1])])
                if p[0]=='L': cs.append(('L',M(p[1]),M(p[2])))
                else:
                    cx,cy,a,b=p[3]; k=bands[i][2]
                    cs.append(('E',M(p[1]),M(p[2]),(cx,f(i,cy),a,b*k),p[4],p[5]))
        out.append(cs)
    return out
def _split(s,yc):
    if s[0]=='L':
        p0,p1=s[1],s[2]
        if (p0[1]-yc)*(p1[1]-yc)<0:
            t=(yc-p0[1])/(p1[1]-p0[1]); m=p0+(p1-p0)*t
            return [('L',p0,m),('L',m,p1)]
        return [s]
    cx,cy,a,b=s[3]; t0,sw=s[4],s[5]
    if abs(yc-cy)>=b: return [s]
    base=np.arcsin((yc-cy)/b); roots=[]
    for r in (base,np.pi-base):
        for kk in range(-3,4):
            t=r+2*np.pi*kk; u=(t-t0)/sw if sw else -1
            if 1e-9<u<1-1e-9: roots.append(u)
    if not roots: return [s]
    roots=sorted(set(np.round(roots,12))); out=[]; prev=0.0; pp=s[1]
    for u in roots+[1.0]:
        ta,tb=t0+sw*prev,t0+sw*u; q=s[2] if u==1.0 else _ept(s[3],tb)
        out.append(('E',pp,q,s[3],ta,tb-ta)); prev=u; pp=q
    return out
# ---------------------------------------------------------------------------------------------------- drawing
def draw(contours,pen,dx=0.0,dy=0.0):
    T=lambda p:(float(p[0]+dx),float(p[1]+dy))
    for segs in contours:
        if not segs: continue
        pen.moveTo(T(segs[0][1]))
        for k,s in enumerate(segs):
            if s[0]=='L':
                if k<len(segs)-1: pen.lineTo(T(s[2]))
            elif s[0]=='B': pen.curveTo(T(s[3]),T(s[4]),T(s[2]))
            else:
                for bz in lcoutline.beziers(s): pen.curveTo(T(bz[1]),T(bz[2]),T(bz[3]))
        pen.closePath()
def to_polygons(contours,n=16,step=None):
    """flattened rings (for checks and overlays); step: largest parameter step in radians (overrides n)"""
    rings=[]
    for segs in contours:
        pts=[]
        for s in segs:
            if s[0]=='L': pts.append(s[1])
            elif s[0]=='B':
                for t in np.linspace(0,1,n,endpoint=False):
                    m_=1-t; pts.append(m_**3*s[1]+3*m_*m_*t*s[3]+3*m_*t*t*s[4]+t**3*s[2])
            else:
                m=n if step is None else max(2,int(np.ceil(abs(s[5])/step)))
                for t in np.linspace(0,1,m,endpoint=False): pts.append(_ept(s[3],s[4]+s[5]*t))
        rings.append(np.array(pts))
    return rings
def clip_y(contours,yc,below=True):
    """cut off whatever crosses y = yc (below it, or above it when below=False): each excursion is replaced by a
    straight line along y = yc between its crossing points"""
    out=[]
    for segs in contours:
        pieces=[q for s in segs for q in _split(s,yc)]
        mid=lambda p: (p[1][1]+p[2][1])/2 if p[0]=='L' else _ept(p[3],p[4]+p[5]/2)[1]
        bad=[(mid(p)<yc-1e-9) if below else (mid(p)>yc+1e-9) for p in pieces]
        if not any(bad): out.append(segs); continue
        n=len(pieces); k=next(i for i in range(n) if not bad[i] and bad[i-1])   # start just after an excursion
        P=[pieces[(k+i)%n] for i in range(n)]; B=[bad[(k+i)%n] for i in range(n)]
        cs=[]; i=0
        while i<n:
            if not B[i]: cs.append(P[i]); i+=1; continue
            j=i
            while j<n and B[j]: j+=1
            cs.append(('L',P[i][1],P[j-1][2])); i=j
        out.append(cs)
    return out

# ---------------------------------------------------------------------------------------------------- condensing
def scale_x(contours,f):
    out=[]
    for segs in contours:
        cs=[]
        for s in segs:
            M=lambda p: np.array([p[0]*f,p[1]])
            if s[0]=='L': cs.append(('L',M(s[1]),M(s[2])))
            else:
                cx,cy,a,b=s[3]; cs.append(('E',M(s[1]),M(s[2]),(cx*f,cy,a*f,b),s[4],s[5]))
        out.append(cs)
    return out
def _carrier(s):
    return ('L',s[1],tangent(s,0)) if s[0]=='L' else ('E',s[3])
def _rebuild(s,p0,p1):
    if s[0]=='L': return ('L',p0,p1)
    e=s[3]; t0=_par(e,p0); t1=_par(e,p1); k=np.sign(s[5])
    sw=(t1-t0)%(2*np.pi) if k>0 else -((t0-t1)%(2*np.pi))
    if abs(s[5])>=2*np.pi-1e-6: sw=s[5]
    return ('E',p0,p1,e,t0,sw)
def reangle(segs,tags,pivot_y=258.75):
    """turn every line tagged with an angle (degrees, direction of travel) back to that angle about its midpoint, and
    re-join it to its neighbours"""
    segs=list(segs); n=len(segs)
    for i,t in enumerate(tags):
        if t is None or segs[i][0]!='L': continue
        s=segs[i]; u=np.array([np.cos(np.radians(t)),np.sin(np.radians(t))])
        # pivot on the end farther from mid x-height: a stem's top corner stays on the ascender line, a descender's
        # bottom corner on its line
        m=s[1] if abs(s[1][1]-pivot_y)>abs(s[2][1]-pivot_y) else s[2]
        if u@(s[2]-s[1])<0: u=-u
        L=('L',m,u)
        a,b=segs[i-1],segs[(i+1)%n]
        p0=_intersect(_carrier(a),L,s[1]); p1=_intersect(L,_carrier(b),s[2])
        if p0 is None or p1 is None: continue
        segs[i-1]=_rebuild(a,a[1],p0); segs[i]=('L',p0,p1); segs[(i+1)%n]=_rebuild(b,p1,b[2])
    return segs
def terminal_tags(segs,maxlen=160.0,turn=15.0,par=0.85):
    """the direction (degrees) of each terminal cut, None elsewhere: a short slanted line meeting both neighbours at
    corners, the neighbours (the stroke's two sides) running roughly parallel"""
    n=len(segs); tags=[None]*n
    for i,s in enumerate(segs):
        if s[0]!='L': continue
        v=s[2]-s[1]; L=np.hypot(*v); ang=np.degrees(np.arctan2(v[1],v[0]))
        if L>maxlen or min(ang%180,abs(ang%180-90),180-ang%180)<1.0: continue
        a=tangent(segs[i-1],1); b=tangent(s,0); c=tangent(s,1); d=tangent(segs[(i+1)%n],0)
        if np.degrees(np.arccos(np.clip(a@b,-1,1)))<turn or np.degrees(np.arccos(np.clip(c@d,-1,1)))<turn: continue
        if abs(a@d)<par: continue
        tags[i]=ang
    return tags
def condense(contours,f,t,tx,ty,X,metric=None,reangle_terms=True):
    """narrow a lower-case letter: offset inwards by t (a thin core), scale x by f, turn the terminal cuts back to their
    angles, map y so the baseline, x-height, ascender and descender come back to their lines once the strokes are
    thickened, then offset outwards by (tx, ty) (horizontal and vertical half-strokes)"""
    out=[]; dots=[]
    for segs in contours:
        tags=terminal_tags(segs) if reangle_terms else [None]*len(segs)
        if len(segs)==1 and segs[0][0]=='E':          # a dot: kept round, its radius changed by the mean half-stroke
            cx,cy,a,b=segs[0][3]; r=a-t+(tx+ty)/2
            yc=cy+(t-ty if cy>X-t else ty-t if cy<t else 0)
            p=np.array([cx*f+r,yc]); dots.append([('E',p,p,(cx*f,yc,r,r),0.0,2*np.pi*np.sign(segs[0][5]))]); continue
        core,tg=offset_contour(segs,-t,tags)
        core,tg=clean(core,tags=tg)
        core=scale_x([core],f)[0]
        core=reangle(core,tg)
        out.append(clean(core))
    k=(X-2*ty)/(X-2*t)
    out=vmap(out,[(-1e9,t,1.0,ty-t),(t,X-t,k,ty-t*k),(X-t,1e9,1.0,t-ty)])
    return [clean(offset_contour(clean(c),(tx,ty))) for c in out]+dots

DIRN=True       # condense2: offsets by edge direction (see there)
PROPORTIONAL=True   # and by the local stroke width (_proportional)
def _seg_mid(sg):
    """(midpoint, outward unit normal) of a segment (contours run with the ink on the right)"""
    if sg[0]=='L':
        p=(sg[1]+sg[2])/2; u=sg[2]-sg[1]
    else:
        t=sg[4]+sg[5]/2; p=_ept(sg[3],t); u=_etan(sg[3],t)*np.sign(sg[5])
    u=u/np.linalg.norm(u); return p,np.array([-u[1],u[0]])
def _even(r):
    """a stroke's width as a fraction of the series stroke, from its source width r (fraction of the source stroke):
    strokes from 0.9 to 1.15 of the source stroke all come out the series stroke, the source's own small variations
    (E(M)'s l 0.98, n 1.02) evened out; thinner ones (where a curve meets a stem) keep their proportion below 0.8,
    with a ramp between.  (Heavier ones, where strokes merge, are left to _proportional: they lose what one of 1.15
    does, so the junctions are not eaten into.)"""
    if r<0.8: return r
    if r<0.9: return 0.8+2*(r-0.8)
    return 1.0
def _proportional(contours,f,k,V,T,ey,X,plain=()):
    """per-segment offset amounts (negative: erosion) that thin each stroke by its width, so that every stroke ends
    as the series weight T (or, where thinner, in proportion: _even): a stroke whose width here (after the narrowing
    and stretching) is w, across an edge of normal n, loses w - T _even(w / (V s)), s = sqrt((f nx)^2 + (k ny)^2).
    The edges on the alignment
    lines (tops in the x-height or ascender zone, bottoms at the baseline or below) move by the full-stroke amount
    (V s - T)/2 that the vertical map assumed, keeping the letter on its lines, and the edge across the stroke from one
    takes the rest; other strokes lose half on each side.  Ends of strokes (w over 1.6 strokes) keep the full-stroke
    amount.  A line is probed at its middle; an arc at its side point and its top or bottom point (or the ends of the
    arc nearest them), the two amounts going to its two semi-axes."""
    from shapely.geometry import Polygon, LineString
    g=Polygon()
    for r in to_polygons(contours,n=24):
        if len(r)>2: g=g.symmetric_difference(Polygon(r).buffer(0))
    bnd=g.boundary
    sk=lambda n: np.sqrt((f*n[0])**2+(k*n[1])**2)
    dV=lambda n: (V*sk(n)-T)/2
    top=X+ey-0.25*V*k; bot=-ey+0.25*V*k
    def is_align(p,no): return (no[1]>0.5 and p[1]>=top) or (no[1]<-0.5 and p[1]<=bot)
    def ray(p,no):
        hit=LineString([p-no*0.01,p-no*(4*V)]).intersection(bnd); w=None; q=None
        for h in getattr(hit,'geoms',[hit]):
            if h.is_empty: continue
            for xy in (h.coords if h.geom_type!='Point' else [h.coords[0]]):
                dd=np.linalg.norm(np.array(xy)-p)
                if dd>0.05 and (w is None or dd<w): w=dd; q=np.array(xy)
        return w,q
    def amount(p,no):
        """erosion (positive) at the edge point p with outward normal no"""
        if is_align(p,no): return dV(no)
        if any(np.hypot(*(p-np.asarray(z)))<r_ for z,r_ in plain): return dV(no)   # near a cut: as the cut is placed
        w,q=ray(p,no)
        if w is None or w>1.6*V*sk(no): return dV(no)
        total=min(w,1.15*V*sk(no))-T*_even(w/(V*sk(no)))
        # the far side: an alignment edge there (its outward normal is about -no) takes the full-stroke amount
        if q is not None and is_align(q,-no): return total-dV(-no)
        return total/2
    out=[]
    for c in contours:
        data=[]
        for sg in c:
            if sg[0]=='L':
                p,no=_seg_mid(sg); data.append([(None,amount(p,no))]+([('pin',dV(no))] if is_align(p,no) else []))
                continue
            e=sg[3]; t0,sw=sg[4],sg[5]; sgn=np.sign(sw)
            ts=[t0+sw*u for u in (0.1,0.5,0.9)+(0.03,0.97)*3]   # near the ends (closest counted thrice: the
                                                     # junctions are there), the middle, and the axis points it passes
            lo,hi=min(t0,t0+sw),max(t0,t0+sw)
            ts+=[jj*np.pi/2 for jj in range(-8,9) if lo<jj*np.pi/2<hi]
            row=[]
            for t in ts:
                p=_ept(e,t); u=_etan(e,t)*sgn; u=u/np.linalg.norm(u); row.append((t,amount(p,np.array([-u[1],u[0]]))))
            data.append(row)
        out.append(_reconcile(c,data,(dV(np.array([1.0,0.0])),dV(np.array([0.0,1.0]))),dV))
    return out
def _arc_coef(a,b,c2,s2):
    """how far a point of an ellipse (a, b) moves along its normal per unit change of each semi-axis, at the point
    with cos^2 t = c2, sin^2 t = s2"""
    r=np.sqrt(b*b*c2+a*a*s2); return np.array([b*c2/r,a*s2/r])
def _reconcile(c,data,prior,dV):
    """per-segment amounts from the probed ones (data: for each segment, (parameter, amount) pairs; parameter None for
    a line), fitted by least squares: a line has one amount, an arc one for each semi-axis, its amount at a point being
    how far that moves the point along its normal; segments that meet smoothly must move their common end alike
    (otherwise the joint is placed between the two moved carriers, re-angling the lines and stepping the arcs: the
    notches in Series A's stems and arches).  prior: the full-stroke amounts (x, y), weakly held where the probes leave
    an arc's other semi-axis free (a short piece near its top or side).  A line's probe counts five times an arc's:
    where they meet, the arc gives way; a line on an alignment line (a stem's foot, the x-height) is held there.  At
    corners up to 60 degrees (and, less firmly, up to 120) the two sides must depart alike from the
    full-stroke amounts dV of their own normals: otherwise their moved carriers cross far from the corner (the
    shift grows as 1 / sin of the angle), leaving slivers where Series B's counters meet the stems and steps in the
    k, and the a's counter bulging past its stem.  Returns the amount functions of the outward normal (negative: erosion)"""
    n=len(c); col=[]; nv=0
    for sg in c: k_=1 if sg[0]=='L' else 2; col.append(list(range(nv,nv+k_))); nv+=k_
    def coef(i,t):
        row=np.zeros(nv); sg=c[i]
        if sg[0]=='L': row[col[i][0]]=1.0; return row
        a,b=sg[3][2],sg[3][3]; row[col[i]]=_arc_coef(a,b,np.cos(t)**2,np.sin(t)**2); return row
    def end(i,e): sg=c[i]; return None if sg[0]=='L' else sg[4]+(sg[5] if e else 0.0)
    A=[]; y=[]
    for i,rows in enumerate(data):
        wt=5.0 if c[i][0]=='L' else 1.0                 # a straight stroke keeps its weight; the arc beside it gives
        for t,v in rows:
            if t=='pin': A.append(1e5*coef(i,None)); y.append(1e5*v); continue      # on an alignment line: stays there
            A.append(wt*coef(i,t)); y.append(wt*v)
        for j,pv in zip(col[i],prior if c[i][0]=='E' else prior[:0]):
            r=np.zeros(nv); r[j]=0.03; A.append(r); y.append(0.03*pv)
    def nrm(u): return np.array([-u[1],u[0]])
    for i in range(n):
        if n<2: break
        ua,ub=tangent(c[i-1],1),tangent(c[i],0)
        if ua@ub<np.cos(np.radians(120.0)): continue          # (150 tangles the m's crotch in A and B)
        wj=1e3 if ua@ub>=np.cos(np.radians(60.0)) else 30.0      # steeper corners: held less firmly
        A.append(wj*(coef(i-1,end(i-1,1))-coef(i,end(i,0)))); y.append(wj*(dV(nrm(ua))-dV(nrm(ub))))
    x=np.linalg.lstsq(np.array(A),np.array(y),rcond=None)[0]
    res=[]
    for i,sg in enumerate(c):
        if sg[0]=='L': res.append(lambda m,d=x[col[i][0]]: -d); continue
        a,b=sg[3][2],sg[3][3]; dx,dy=x[col[i]]
        def fa(m,a=a,b=b,dx=dx,dy=dy):
            cc=(m[0]/b)**2; ss=(m[1]/a)**2; tot=max(cc+ss,1e-24)
            return -float(_arc_coef(a,b,cc/tot,ss/tot)@[dx,dy])
        res.append(fa)
    return res
def _jogs(c,short=10.0):
    """a contour with its jogs taken out: a run of lines shorter than `short` between two longer lines (left by the
    corner offsets where the k's arm meets its leg and stem) becomes the longer lines' intersection, when that is
    near; a speck (a contour no longer than a jog) is dropped"""
    c=list(c)
    if sum(np.hypot(*(s[2]-s[1])) for s in c)<4*short: return None
    L=lambda s: s[0]=='L' and np.hypot(*(s[2]-s[1]))<short
    changed=True
    while changed and len(c)>3:
        changed=False; n=len(c)
        for i in range(n):
            if not L(c[i]) or L(c[i-1]) or c[i-1][0]!='L': continue
            j=i
            while L(c[j%n]) and j-i<n-2: j+=1
            b=c[j%n]; a=c[i-1]
            if b[0]!='L': continue
            u=a[2]-a[1]; v=b[2]-b[1]; den=u[0]*v[1]-u[1]*v[0]
            if abs(den)<1e-9*np.hypot(*u)*np.hypot(*v): continue
            t=((b[1]-a[1])[0]*v[1]-(b[1]-a[1])[1]*v[0])/den; q=a[1]+u*t
            if np.hypot(*(q-a[2]))>3*short or np.hypot(*(q-b[1]))>3*short: continue
            run=[k%n for k in range(i,j)]
            c[i-1]=('L',a[1],q); c[j%n]=('L',q,b[2])
            c=[s for k,s in enumerate(c) if k not in run]; changed=True; break
    return c
def trim_jogs(c,short=15.0,reach=30.0):
    """short lines left between an arc and a long line (where a cut was remade and its ears removed: the t's tail,
    the a's terminal) taken out: the long line is extended to the arc's ellipse, if that is within reach; between
    two lines, the two meet"""
    c=list(c); changed=True
    ln_=lambda sg: np.hypot(*(np.array(sg[2],float)-np.array(sg[1],float)))
    def meet(e,p,d,end):
        """the point of the line p + d s on ellipse e nearest `end`, or None"""
        cx,cy,a,b=e[:4]; q0=np.array([(p[0]-cx)/a,(p[1]-cy)/b]); q1=np.array([d[0]/a,d[1]/b])
        A2,B2,C2=q1@q1,2*q0@q1,q0@q0-1; disc=B2*B2-4*A2*C2
        if disc<0: return None
        return min([p+d*s_ for s_ in ((-B2-np.sqrt(disc))/(2*A2),(-B2+np.sqrt(disc))/(2*A2))],key=lambda r: np.hypot(*(r-end)))
    # first: a run of lines (together under 2 short) between an arc and a long line (a remnant of the edge beside a
    # remade cut, with its jogs: Series A's e) goes, the long line extended to the arc
    while changed and len(c)>4:
        changed=False; n=len(c)
        for i in range(n):
            if c[i][0]!='E': continue
            for dr in (1,-1):
                run=[]; j=(i+dr)%n
                while c[j][0]=='L' and ln_(c[j])<1.5*short and len(run)<4: run.append(j); j=(j+dr)%n
                if len(run)<2 or sum(ln_(c[r]) for r in run)>=2*short or c[j][0]!='L' or ln_(c[j])<3*short: continue
                A_,L_=c[i],c[j]; p=np.array(L_[1],float); d=np.array(L_[2],float)-p
                end=np.array(A_[2] if dr==1 else A_[1],float); q=meet(A_[3],p,d,end)
                if q is None or np.hypot(*(q-end))>reach: continue
                e=A_[3]; tq=_par(e[:4],q); sw=A_[5]
                if dr==1:
                    t0=A_[4]; nsw=(tq-t0)%(2*np.pi) if sw>0 else -((t0-tq)%(2*np.pi))
                    if abs(nsw-sw)>0.5: continue
                    c[i]=('E',A_[1],q,e,t0,nsw); c[j]=('L',q,np.array(L_[2],float))
                else:
                    t1=A_[4]+sw; nsw=(t1-tq)%(2*np.pi) if sw>0 else -((tq-t1)%(2*np.pi))
                    if abs(nsw-sw)>0.5: continue
                    c[i]=('E',q,A_[2],e,tq,nsw); c[j]=('L',np.array(L_[1],float),q)
                c=[sg for k,sg in enumerate(c) if k not in run]; changed=True; break
            if changed: break
    changed=True
    while changed and len(c)>3:
        changed=False; n=len(c)
        for i in range(n):
            sg=c[i]
            if sg[0]!='L' or np.hypot(*(np.array(sg[2])-np.array(sg[1])))>=short: continue
            P_,N_=c[(i-1)%n],c[(i+1)%n]
            if P_[0]=='L' and N_[0]=='L':                   # line, jog, line: the two lines meet
                p=np.array(P_[1],float); u=np.array(P_[2],float)-p; q_=np.array(N_[1],float); v=np.array(N_[2],float)-q_
                den=u[0]*v[1]-u[1]*v[0]
                if abs(den)>1e-9*np.hypot(*u)*np.hypot(*v):
                    t_=((q_-p)[0]*v[1]-(q_-p)[1]*v[0])/den; X_=p+u*t_
                    if np.hypot(*(X_-np.array(P_[2],float)))<reach and np.hypot(*(X_-q_))<reach and 0<t_:
                        c[(i-1)%n]=('L',p,X_); c[(i+1)%n]=('L',X_,np.array(N_[2],float)); del c[i]; changed=True; break
                continue
            for arc_i,line_i in (((i-1)%n,(i+1)%n),((i+1)%n,(i-1)%n)):
                A_,L_=c[arc_i],c[line_i]
                if A_[0]!="E" or L_[0]!="L" or np.hypot(*(np.array(L_[2])-np.array(L_[1])))<1.5*short: continue
                e=A_[3]; p=np.array(L_[1],float); d=np.array(L_[2],float)-p
                cx,cy,a,b=e[:4]
                # points of the line p + d s on the ellipse: ((x-cx)/a)^2 + ((y-cy)/b)^2 = 1
                q0=np.array([(p[0]-cx)/a,(p[1]-cy)/b]); q1=np.array([d[0]/a,d[1]/b])
                A2,B2,C2=q1@q1,2*q0@q1,q0@q0-1; disc=B2*B2-4*A2*C2
                end=np.array(A_[2] if arc_i==(i-1)%n else A_[1],float)
                if disc>=0:
                    cand=[p+d*s_ for s_ in ((-B2-np.sqrt(disc))/(2*A2),(-B2+np.sqrt(disc))/(2*A2))]
                    q=min(cand,key=lambda r: np.hypot(*(r-end)))
                    newline=None
                else:                                   # the ellipse falls short of the line (the r's arm end): the
                    nl=np.array([-d[1],d[0]])/np.linalg.norm(d)     # arc ends at its point nearest the line, and the
                    ts=np.linspace(0,2*np.pi,721)                   # line turns a little to reach it
                    P=np.c_[cx+a*np.cos(ts),cy+b*np.sin(ts)]; dist=np.abs((P-p)@nl)
                    q=P[np.argmin(dist)]
                    if dist.min()>6.0: continue
                    newline=True
                if np.hypot(*(q-end))>reach: continue
                tq=_par(e[:4],q); sw=A_[5]
                if arc_i==(i-1)%n:                      # arc, jog, line: the arc's end moves to q
                    t0=A_[4]; nsw=(tq-t0)%(2*np.pi) if sw>0 else -((t0-tq)%(2*np.pi))
                    if abs(nsw-sw)>(1.0 if newline else 0.5): continue
                    c[arc_i]=('E',A_[1],q,e,t0,nsw); c[line_i]=('L',q,np.array(L_[2],float))
                    if newline and np.hypot(*(np.array(L_[2],float)-q))<1.5*short: c[arc_i],c[line_i]=A_,L_; continue
                else:                                   # line, jog, arc: the arc's start moves to q
                    t1=A_[4]+sw; nsw=(t1-tq)%(2*np.pi) if sw>0 else -((tq-t1)%(2*np.pi))
                    if abs(nsw-sw)>(1.0 if newline else 0.5): continue
                    c[arc_i]=('E',q,A_[2],e,tq,nsw); c[line_i]=('L',np.array(L_[1],float),q)
                    if newline and np.hypot(*(np.array(L_[1],float)-q))<1.5*short: c[arc_i],c[line_i]=A_,L_; continue
                del c[i]; changed=True; break
            if changed: break
    return c
def condense2(contours,f,gx,ey,X,mitre=None,dot_scale=None,proportional=None,plain=()):
    """narrow by f with no thin core: scale x by f, offset horizontally by gx (vertical strokes back to weight),
    stretch y so that, after the vertical offset by -ey (horizontal strokes to weight), the baseline and x-height are
    back on their lines; dots are kept round"""
    out=[]; dots=[]
    for segs in contours:
        if len(segs)==1 and segs[0][0]=='E':
            cx,cy,a,b=segs[0][3]; dots.append((cx,cy,a,np.sign(segs[0][5]))); continue
        out.append(scale_x([segs],f)[0])
    global MITRE
    k=(X+2*ey)/X
    if DIRN:
        # narrowed and stretched first, then every edge moved along its own normal n by what a stroke of the source
        # weight V needs to become T there: its width is V sqrt(f^2 nx^2 + k^2 ny^2) (normals transform by the inverse
        # transpose), so each side moves (V sqrt(...) - T)/2.  For vertical and horizontal edges that is -gx and ey
        # as before; diagonals no longer lose weight (Series A's z was 0.72 of the stroke)
        V=2*(ey+gx)/(k-f); T=k*V-2*ey
        out=vmap(out,[(-1e9,0.0,1.0,-ey),(0.0,X,k,-ey),(X,1e9,1.0,ey)])
        out=[clean(c) for c in out]
        amounts=_proportional(out,f,k,V,T,ey,X,plain) if (PROPORTIONAL if proportional is None else proportional) else \
                [[(lambda n: -(V*np.sqrt((f*n[0])**2+(k*n[1])**2)-T)/2)]*len(c) for c in out]
        steps=max(1,int(np.ceil(max(abs(gx),abs(ey))/20.0)))
        MITRE=mitre
        try:
            res=[]
            for c,A in zip(out,amounts):
                tg=list(A)
                for _ in range(steps):
                    c,tg=clean(c,tags=tg)
                    c,tg=offset_contour(c,[(lambda n,a=a: a(n)/steps) for a in tg],tags=tg)
                c,tg=clean(c,tags=tg); res.append(c)
            out=[t for t in (_jogs(c) for c in res) if t]
        finally: MITRE=None
    else:
        MITRE=mitre
        try: out=[clean(offset_contour(clean(c),(gx,0.0))) for c in out]
        finally: MITRE=None
        out=vmap(out,[(-1e9,0.0,1.0,-ey),(0.0,X,k,-ey),(X,1e9,1.0,ey)])
        steps=max(1,int(np.ceil(abs(ey)/20.0)))      # in steps of at most 20 units: a feature that a single large offset
        for _ in range(steps):                        # collapses (Series A's a with its stub) is re-joined as it goes
            out=[clean(offset_contour(clean(c),(0.0,-ey/steps))) for c in out]
    for cx,cy,a,sg in dots:
        # dots scale with the stroke (dot_scale = series stroke / source stroke) and keep their top where it was
        r=a*dot_scale if dot_scale else a+(gx-ey)/2
        yc=cy+a-r
        p=np.array([cx*f+r,yc]); out.append([('E',p,p,(cx*f,yc,r,r),0.0,2*np.pi*sg)])
    return out
# ---------------------------------------------------------------------------------------------------- squarer arcs
SQ_BANDS=None   # vmap bands the contours went through: squaring is done before them (set by squarer)
def _bandfn(bands):
    def f(y):
        for i,(y0,y1,k,c) in enumerate(bands):
            if (y<y1 or i==len(bands)-1) and (y>=y0 or i==0): return c+k*y
    return f
def _sq_pt(e,t,n):
    """the point of parameter t on ellipse e = (cx, cy, a, b) moved to the curve |x/a|^2 + |y/b|^n = 1: x kept,
    y = cy + b sign(sin t) |sin t|^(2/n).  With SQ_BANDS the ellipse is a piece of one that went through that vmap:
    it is squared where it was before the map (its band's inverse) and mapped back by the whole map, so the pieces
    of one ellipse split at a band edge stay joined"""
    cx,cy,a,b=e[:4]; s=np.sin(t); S=np.sign(s)*abs(s)**(2.0/n)
    if SQ_BANDS is None or len(e)<5: return np.array([cx+a*np.cos(t),cy+b*S])
    bands,fwd,ev=SQ_BANDS
    ym=cy+b*s                                               # the unsquared point decides the piece's band
    for i,(y0,y1,k,c) in enumerate(bands):
        lo,hi=c+k*y0,c+k*y1
        if (ym<hi or i==len(bands)-1) and (ym>=lo or i==0): break
    dl=e[4]*ev                                              # the vertical erosion's change of the semi-axis (+ for
    cyp=(cy-c)/k; bp=(b-dl)/k                               # arcs travelled anticlockwise, counters), undone
    return np.array([cx+a*np.cos(t),fwd(cyp+bp*S)+dl*S])
def _fit_cubic(P,p0,p3,t0,t3):
    """cubic from p0 to p3 with unit end tangents t0, t3 (handles p0 + al t0, p3 - be t3) fitted to the points P
    (chord-length parameters); returns the four control points and the largest distance from P"""
    d=np.r_[0,np.cumsum(np.linalg.norm(np.diff(P,axis=0),axis=1))]; u=d/d[-1]
    B0=(1-u)**3; B1=3*u*(1-u)**2; B2=3*u*u*(1-u); B3=u**3
    R=P-np.outer(B0+B1,p0)-np.outer(B2+B3,p3)
    A=np.stack([np.outer(B1,t0).ravel(),np.outer(B2,-t3).ravel()],1)
    al,be=np.linalg.lstsq(A,R.ravel(),rcond=None)[0]
    L=np.linalg.norm(p3-p0); al=al if al>1e-3*L else L/3; be=be if be>1e-3*L else L/3
    c=[p0,p0+al*t0,p3-be*t3,p3]
    Q=np.outer(B0,c[0])+np.outer(B1,c[1])+np.outer(B2,c[2])+np.outer(B3,c[3])
    return c,np.linalg.norm(Q-P,axis=1).max()
def _sq_beziers(e,ta,tb,n,pa,pb,tol=0.25,depth=0):
    """cubics for the squared arc from parameter ta to tb, ending exactly at pa and pb"""
    ts=np.linspace(ta,tb,41); P=np.array([_sq_pt(e,t,n) for t in ts]); P[0],P[-1]=pa,pb
    def tan(t,dirn):
        h=1e-5*dirn; q=_sq_pt(e,t+h,n)-_sq_pt(e,t,n); return q/np.linalg.norm(q)
    sg=np.sign(tb-ta)
    c,err=_fit_cubic(P,pa,pb,tan(ta,sg),-tan(tb,-sg))
    if err<=tol or depth>6: return [c]
    tm=(ta+tb)/2; pm=_sq_pt(e,tm,n)
    return _sq_beziers(e,ta,tm,n,pa,pm,tol,depth+1)+_sq_beziers(e,tm,tb,n,pm,pb,tol,depth+1)
def _axis(t): return min(abs(np.sin(t)),1-abs(np.sin(t)))<1e-9
def _sq_line(e,t,n,q,d,win=0.5):
    """the parameter within win radians of the arc's end t where its squared curve crosses the line through q
    along d, the crossing nearest t; t itself if there is none"""
    from scipy.optimize import brentq
    def f(u): w=_sq_pt(e,u,n)-q; return float(d[0]*w[1]-d[1]*w[0])
    us=np.linspace(t-win,t+win,401); vs=[f(u) for u in us]; best=None
    for i in range(400):
        if vs[i]==0: r=us[i]
        elif vs[i]*vs[i+1]<0: r=brentq(f,us[i],us[i+1],xtol=1e-12)
        else: continue
        if best is None or abs(r-t)<abs(best-t): best=r
    return t if best is None else best
def _sq_arcs(e1,t1,e2,t2,n,lo1=None,hi1=None,lo2=None,hi2=None):
    """where the squared curves of two arcs cross, near their ends t1 and t2: their sampled curves (each over its arc
    and 0.6 radians beyond the end) are intersected and the crossing nearest the squared ends kept, then refined"""
    from scipy.optimize import fsolve
    def samp(e,t,lo,hi):
        a_=min(lo,t-0.6) if lo is not None else t-0.6; b_=max(hi,t+0.6) if hi is not None else t+0.6
        u=np.linspace(a_,b_,400); return u,np.array([_sq_pt(e,v,n) for v in u])
    u1,P1=samp(e1,t1,lo1,hi1); u2,P2=samp(e2,t2,lo2,hi2)
    goal=(_sq_pt(e1,t1,n)+_sq_pt(e2,t2,n))/2; best=None
    for i in range(len(P1)-1):
        a0,a1=P1[i],P1[i+1]; da=a1-a0
        d2=P2[1:]-P2[:-1]; w=P2[:-1]-a0
        den=da[0]*d2[:,1]-da[1]*d2[:,0]; ok=np.abs(den)>1e-12
        sa=np.where(ok,(w[:,0]*d2[:,1]-w[:,1]*d2[:,0])/np.where(ok,den,1),-1)
        sb=np.where(ok,(w[:,0]*da[1]-w[:,1]*da[0])/np.where(ok,den,1),-1)
        hit=np.flatnonzero(ok&(sa>=0)&(sa<=1)&(sb>=0)&(sb<=1))
        for j in hit:
            q=a0+sa[j]*da; dist=np.linalg.norm(q-goal)
            if best is None or dist<best[0]: best=(dist,u1[i]+sa[j]*(u1[i+1]-u1[i]),u2[j]+sb[j]*(u2[j+1]-u2[j]))
    if best is None: return (t1,t2)
    r=fsolve(lambda v: _sq_pt(e1,v[0],n)-_sq_pt(e2,v[1],n),[best[1],best[2]],xtol=1e-12)
    if np.linalg.norm(_sq_pt(e1,r[0],n)-_sq_pt(e2,r[1],n))<1e-6: return (r[0],r[1])
    return (best[1],best[2])
def squarer(contours,n,keep_circles=True,bands=None,erosion=0.0):
    """see _squarer.  bands, erosion: the vmap and the vertical erosion (positive) the contours went through"""
    global SQ_BANDS
    SQ_BANDS=(bands,_bandfn(bands),erosion) if bands else None
    try:
        if bands:   # each arc's ellipse carries the sign of its travel (which way the erosion moved it)
            contours=[[(sg[0],sg[1],sg[2],tuple(sg[3][:4])+(float(np.sign(sg[5])),),sg[4],sg[5])
                       if sg[0]=='E' and not (keep_circles and abs(sg[3][2]-sg[3][3])<1e-6) else sg
                       for sg in segs] for segs in contours]
        return _squarer(contours,n,keep_circles)
    finally: SQ_BANDS=None
FILLET=5.0     # squarer: the length of each side of a smoothed join
def _squarer(contours,n,keep_circles=True):
    """every elliptic arc to the curve |x/a|^2 + |y/b|^n = 1 about its own centre (n > 2: straighter sides, blunter
    tops and bottoms, as the condensed capital O; the axis points stay, with their tangents); arcs become cubic
    segments ('B', p0, p3, p1, p2).  Circles (a = b, the dots) are left round.  An arc that ended elsewhere than at
    an axis point is re-ended where its squared curve meets the neighbouring line (whose direction is kept, so cuts
    keep their angles and the e's bar its height) or the neighbouring squared arc."""
    out=[]
    for segs in contours:
        # arcs only a few units long (pieces left between neighbouring arcs of slightly different ellipses) become
        # straight bridges: squaring them about their own centres could fold them back
        segs=[('L',np.array(s[1],float),np.array(s[2],float)) if s[0]=='E' and len(segs)>1 and
              np.linalg.norm(np.array(s[2],float)-np.array(s[1],float))<3.0 else s for s in segs]
        m=len(segs)
        def is_sq(s): return s[0]=='E' and not (keep_circles and abs(s[3][2]-s[3][3])<1e-6)
        T=[[s[4],s[4]+s[5]] if is_sq(s) else None for s in segs]       # arc end parameters, updated below
        V=[None]*m; collapse={}; SNAP={}; done=set()
        # runs of short lines between two arcs (left where arcs of different ellipses met, or vmap split one): the
        # second arc ends at its squared curve's nearest point to the first's squared end, and the run collapses there
        def short(sg): return sg[0]=='L' and np.linalg.norm(np.array(sg[2],float)-np.array(sg[1],float))<12
        from scipy.optimize import minimize_scalar
        for i0 in range(m):
            if not is_sq(segs[i0]) or not short(segs[(i0+1)%m]): continue
            j=(i0+1)%m; run=[]
            while short(segs[j]) and len(run)<m: run.append(j); j=(j+1)%m
            if not is_sq(segs[j]) or j==i0: continue
            pa=_sq_pt(segs[i0][3],T[i0][1],n); tb0=T[j][0]
            r=minimize_scalar(lambda u: np.linalg.norm(_sq_pt(segs[j][3],u,n)-pa),bounds=(tb0-0.3,tb0+0.3),
                              method='bounded',options={'xatol':1e-10})
            T[j][0]=r.x; pm=(pa+_sq_pt(segs[j][3],r.x,n))/2
            for v in [(i0+1)%m]+[(q+1)%m for q in run]: SNAP[(v,0)]=pm; done.add(v)
        # a line tangent to arcs at both ends (the s's spine): it becomes the common tangent of the two squared arcs,
        # so that the joins stay smooth (intersecting the squared curves with the old line left corners, more the
        # squarer the series)
        from scipy.optimize import fsolve
        def sqtan(e,t,h=1e-6): return (_sq_pt(e,t+h,n)-_sq_pt(e,t-h,n))/(2*h)
        def cr(u,v): return u[0]*v[1]-u[1]*v[0]
        for k in range(m):
            i,j=(k-1)%m,(k+1)%m; ln=segs[k]
            if ln[0]!='L' or short(ln) or i==j or not (is_sq(segs[i]) and is_sq(segs[j])): continue
            if k in done or j in done or _axis(T[i][1]) or _axis(T[j][0]): continue
            def tg(s,end):
                if s[0]=='L': v=np.array(s[2],float)-np.array(s[1],float)
                else: v=_etan(s[3][:4],s[4]+(s[5] if end else 0.0))*np.sign(s[5])
                return v/np.linalg.norm(v)
            if tg(segs[i],1)@tg(ln,0)<np.cos(np.radians(2.0)) or tg(ln,1)@tg(segs[j],0)<np.cos(np.radians(2.0)): continue
            ea,eb=segs[i][3],segs[j][3]; u0,v0=T[i][1],T[j][0]
            def F(x):
                d=_sq_pt(eb,x[1],n)-_sq_pt(ea,x[0],n); L=max(np.linalg.norm(d),1e-9)
                return [cr(sqtan(ea,x[0]),d)/np.linalg.norm(sqtan(ea,x[0]))/L,cr(sqtan(eb,x[1]),d)/np.linalg.norm(sqtan(eb,x[1]))/L]
            x=fsolve(F,[u0,v0],xtol=1e-12)
            if max(abs(x[0]-u0),abs(x[1]-v0))>0.4 or max(map(abs,F(x)))>1e-7: continue     # (by the residual:
                                                     # fsolve can report slow progress on a converged solution)
            d=_sq_pt(eb,x[1],n)-_sq_pt(ea,x[0],n)
            if d@(np.array(ln[2],float)-np.array(ln[1],float))<=0: continue
            T[i][1],T[j][0]=x; V[k]=_sq_pt(ea,x[0],n); V[j]=_sq_pt(eb,x[1],n); done.update({k,j})
        for k in range(m):                                   # vertex k: end of segment k-1 = start of segment k
            if k in done: continue
            i,j=(k-1)%m,k; a,b=segs[i],segs[j]
            if is_sq(a) and is_sq(b):
                if not (_axis(T[i][1]) and _axis(T[j][0])): T[i][1],T[j][0]=_sq_arcs(a[3],T[i][1],b[3],T[j][0],n)
                V[k]=_sq_pt(b[3],T[j][0],n)
            elif is_sq(a) or is_sq(b):
                arc,ia,end,ln=(a,i,1,b) if is_sq(a) else (b,j,0,a)
                t=T[ia][end]
                beyond=segs[(j+1)%m] if is_sq(a) else segs[(i-1)%m]              # the segment past the line
                sliver=np.linalg.norm(np.array(ln[2],float)-np.array(ln[1],float))<12
                blen=np.linalg.norm(np.array(beyond[2],float)-np.array(beyond[1],float))
                if sliver and is_sq(beyond):                         # a sliver between two arcs (often one ellipse split
                    bi=(j+1)%m if is_sq(a) else (i-1)%m                # by vmap): both arc ends meet at the mean of their
                    bend=0 if is_sq(a) else 1                          # squared ends, the sliver collapsing
                    # the two ends are re-chosen on the squared curves to meet: the sliver's own arc keeps its end and
                    # the other arc ends at its squared curve's nearest point to it (the squaring magnifies small
                    # differences of parameter near the sides of an ellipse)
                    from scipy.optimize import minimize_scalar
                    pa=_sq_pt(arc[3],T[ia][end],n); tb0=T[bi][bend]
                    r=minimize_scalar(lambda u: np.linalg.norm(_sq_pt(beyond[3],u,n)-pa),bounds=(tb0-0.2,tb0+0.2),
                                      method='bounded',options={'xatol':1e-10})
                    T[bi][bend]=r.x; pm=(pa+_sq_pt(beyond[3],r.x,n))/2
                    SNAP[(ia,end)]=pm; SNAP[(bi,bend)]=pm
                elif sliver and blen<12: pass                          # a chain of slivers: it bridges
                elif not _axis(t):
                    q=np.array(ln[1],float); d=np.array(ln[2],float)-q
                    if np.linalg.norm(d)<12 and beyond[0]=='L':                 # a sliver (a cut's join): use the cut
                        q=np.array(beyond[1],float); d=np.array(beyond[2],float)-q
                        collapse[(k+1)%m if is_sq(a) else (k-1)%m]=k              # and shrink the sliver to nothing
                    T[ia][end]=_sq_line(arc[3],t,n,q,d)
                V[k]=_sq_pt(arc[3],T[ia][end],n)
            else: V[k]=np.array(b[1],float)
        for v,k in collapse.items(): V[v]=V[k]
        for (si,e_),pm in SNAP.items(): V[si if e_==0 else (si+1)%m]=pm
        # joins that were smooth before squaring and are corners after it (two squared arcs, each about its own
        # centre, meet at an angle: where the s's counters end in their small ellipses): each side is backed off by
        # FILLET units and the gap bridged by a cubic along both sides' directions
        FIL={}; END={}
        snapped={si if e_==0 else (si+1)%m for (si,e_) in SNAP}|set(collapse)|done
        def stg(sg,end):
            if sg[0]=='L': v=np.array(sg[2],float)-np.array(sg[1],float)
            else: v=_etan(sg[3][:4],sg[4]+(sg[5] if end else 0.0))*np.sign(sg[5])
            return v/max(np.linalg.norm(v),1e-12)
        for k in range(m):
            i=(k-1)%m; a,b=segs[i],segs[k]
            if k in snapped or m<3 or V[k] is None or not (is_sq(a) or is_sq(b)): continue
            if stg(a,1)@stg(b,0)<np.cos(np.radians(2.0)): continue                       # a corner by design
            if is_sq(a) and is_sq(b) and _axis(T[i][1]) and _axis(T[k][0]): continue      # already smooth
            def back(idx,end):
                """the point FILLET units from vertex k along segment idx, its parameter, and the direction of travel"""
                sg=segs[idx]
                if not is_sq(sg):
                    q0=np.array(sg[1],float) if sg[0]=='L' else None
                    if q0 is None: return None
                    P0,P1=V[idx],V[(idx+1)%m]; L=np.linalg.norm(P1-P0)
                    if L<3*FILLET: return None
                    u=(P1-P0)/L; return (P1-u*FILLET if end else P0+u*FILLET),None,u
                t0,t1=T[idx]; sw=t1-t0
                if abs(sw)<1e-6: return None
                Pv=_sq_pt(sg[3],t1 if end else t0,n)
                lo,hi=0.0,abs(sw)/2
                if np.linalg.norm(_sq_pt(sg[3],(t1-np.sign(sw)*hi) if end else (t0+np.sign(sw)*hi),n)-Pv)<2*FILLET: return None
                for _ in range(60):
                    md=(lo+hi)/2; t=(t1-np.sign(sw)*md) if end else (t0+np.sign(sw)*md)
                    if np.linalg.norm(_sq_pt(sg[3],t,n)-Pv)<FILLET: lo=md
                    else: hi=md
                t=(t1-np.sign(sw)*lo) if end else (t0+np.sign(sw)*lo); h=1e-6*np.sign(sw)
                d=_sq_pt(sg[3],t+h,n)-_sq_pt(sg[3],t-h,n); return _sq_pt(sg[3],t,n),t,d/np.linalg.norm(d)
            ra,rb=back(i,1),back(k,0)
            if ra is None or rb is None: continue
            pa,ta,da=ra; pb,tb,db=rb
            if da@db<np.cos(np.radians(30.0)): continue
            if ta is not None: T[i][1]=ta
            if tb is not None: T[k][0]=tb
            L=np.linalg.norm(pb-pa)/3; FIL[k]=('B',pa,pb,pa+da*L,pb-db*L); END[k]=pa; V[k]=pb
        cs=[]
        for k,s in enumerate(segs):
            if k in FIL: cs.append(FIL[k])
            p0,p1=V[k],END.get((k+1)%m,V[(k+1)%m])
            if not is_sq(s):
                if s[0]=='L' and np.linalg.norm(p1-p0)<1e-9: continue
                cs.append(('L',p0,p1) if s[0]=='L' else (s[0],p0,p1)+tuple(s[3:])); continue
            e=s[3]; t0,t1=T[k]
            cuts=sorted(t for jj in range(-8,9) for t in [jj*np.pi/2] if min(t0,t1)+1e-9<t<max(t0,t1)-1e-9)
            if t1<t0: cuts=cuts[::-1]
            ts=[t0]+cuts+[t1]; pts=[p0]+[_sq_pt(e,t,n) for t in cuts]+[p1]
            for ii in range(len(ts)-1):
                if np.linalg.norm(pts[ii+1]-pts[ii])<1e-6: continue     # nothing left of it
                for c in _sq_beziers(e,ts[ii],ts[ii+1],n,pts[ii],pts[ii+1]):
                    cs.append(('B',c[0],c[3],c[1],c[2]))
        cs=[sg for sg in cs if not (sg[0]=="B" and max(np.linalg.norm(np.array(q)-sg[1]) for q in sg[2:5])<1e-3)]
        out.append(cs)
    return out
