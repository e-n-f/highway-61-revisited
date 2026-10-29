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
    dx,dy=(d,d) if np.isscalar(d) else d
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
        if tags is not None: tags=[t for i,t in enumerate(tags) if i!=bad[0]]
    return (out,tags) if tags is not None else out
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
    n=len(segs); C=[_moved(s,d) for s in segs]
    J=[]; BEV={}
    for i in range(n):
        a,b=segs[i-1],segs[i]; P=a[2]                         # joint between segment i-1 and i
        ua,ub=tangent(a,1),tangent(b,0)
        na=np.array([-ua[1],ua[0]]); nb=np.array([-ub[1],ub[0]])
        smooth=np.degrees(np.arccos(np.clip(ua@ub,-1,1)))<0.5
        if smooth:
            J.append(P+_support((na+nb)/np.linalg.norm(na+nb),d)); continue
        # a corner: intersect the moved carriers near the mitre point
        dd=d if np.isscalar(d) else max(abs(d[0]),abs(d[1]))*(np.sign(d[0]) if d[0]!=0 else np.sign(d[1]))
        m=na+nb; g=P+m/max(np.linalg.norm(m),1e-9)*dd/max(np.cos(np.arccos(np.clip(na@nb,-1,1))/2),0.2)
        q=_intersect(C[i-1],C[i],g)
        q=q if q is not None else g
        if MITRE is not None and np.hypot(*(q-P))>MITRE*abs(dd):
            pa=P+_support(na,d); pb=P+_support(nb,d); BEV[i]=(pa,pb); J.append(None); continue
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

def condense2(contours,f,gx,ey,X,mitre=None,dot_scale=None):
    """narrow by f with no thin core: scale x by f, offset horizontally by gx (vertical strokes back to weight),
    stretch y so that, after the vertical offset by -ey (horizontal strokes to weight), the baseline and x-height are
    back on their lines; dots are kept round"""
    out=[]; dots=[]
    for segs in contours:
        if len(segs)==1 and segs[0][0]=='E':
            cx,cy,a,b=segs[0][3]; dots.append((cx,cy,a,np.sign(segs[0][5]))); continue
        out.append(scale_x([segs],f)[0])
    global MITRE
    MITRE=mitre
    try: out=[clean(offset_contour(clean(c),(gx,0.0))) for c in out]
    finally: MITRE=None
    k=(X+2*ey)/X
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
def _sq_pt(e,t,n):
    """the point of parameter t on ellipse e = (cx, cy, a, b) moved to the curve |x/a|^2 + |y/b|^n = 1: x kept,
    y = cy + b sign(sin t) |sin t|^(2/n)"""
    cx,cy,a,b=e; s=np.sin(t)
    return np.array([cx+a*np.cos(t),cy+b*np.sign(s)*abs(s)**(2.0/n)])
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
def _sq_arcs(e1,t1,e2,t2,n):
    from scipy.optimize import fsolve
    r=fsolve(lambda v: _sq_pt(e1,v[0],n)-_sq_pt(e2,v[1],n),[t1,t2],xtol=1e-12)
    return (r[0],r[1]) if np.linalg.norm(_sq_pt(e1,r[0],n)-_sq_pt(e2,r[1],n))<1e-6 else (t1,t2)
def squarer(contours,n,keep_circles=True):
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
        V=[None]*m; collapse={}
        for k in range(m):                                   # vertex k: end of segment k-1 = start of segment k
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
                if sliver and (is_sq(beyond) or blen<12): pass       # a sliver between arcs or more slivers: it bridges
                elif not _axis(t):
                    q=np.array(ln[1],float); d=np.array(ln[2],float)-q
                    if np.linalg.norm(d)<12 and beyond[0]=='L':                 # a sliver (a cut's join): use the cut
                        q=np.array(beyond[1],float); d=np.array(beyond[2],float)-q
                        collapse[(k+1)%m if is_sq(a) else (k-1)%m]=k              # and shrink the sliver to nothing
                    T[ia][end]=_sq_line(arc[3],t,n,q,d)
                V[k]=_sq_pt(arc[3],T[ia][end],n)
            else: V[k]=np.array(b[1],float)
        for v,k in collapse.items(): V[v]=V[k]
        cs=[]
        for k,s in enumerate(segs):
            p0,p1=V[k],V[(k+1)%m]
            if not is_sq(s):
                if s[0]=='L' and np.linalg.norm(p1-p0)<1e-9: continue
                cs.append(('L',p0,p1) if s[0]=='L' else (s[0],p0,p1)+tuple(s[3:])); continue
            e=s[3]; t0,t1=T[k]
            cuts=sorted(t for jj in range(-8,9) for t in [jj*np.pi/2] if min(t0,t1)+1e-9<t<max(t0,t1)-1e-9)
            if t1<t0: cuts=cuts[::-1]
            ts=[t0]+cuts+[t1]; pts=[p0]+[_sq_pt(e,t,n) for t in cuts]+[p1]
            for ii in range(len(ts)-1):
                for c in _sq_beziers(e,ts[ii],ts[ii+1],n,pts[ii],pts[ii+1]):
                    cs.append(('B',c[0],c[3],c[1],c[2]))
        out.append(cs)
    return out
