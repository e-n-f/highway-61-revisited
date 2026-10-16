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
def _moved(s,d):
    """the moved carrier of a segment: ('L', point, direction) or ('E', ellipse)"""
    if s[0]=='L':
        u=tangent(s,0); n=np.array([-u[1],u[0]])            # left of travel = away from the ink
        return ('L',s[1]+n*d,u)
    cx,cy,a,b=s[3]; k=np.sign(s[5])                         # travelling anticlockwise: the centre is on the left
    return ('E',(cx,cy,a-k*d,b-k*d))
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
def offset_contour(segs,d):
    """offset one contour; segments that the offset makes vanish or run backwards are dropped and their neighbours
    re-joined, until none are left"""
    segs=list(segs)
    for _ in range(len(segs)):
        out=_offset_once(segs,d)
        if len(segs)<=2: return out
        bad=[i for i,(s,o) in enumerate(zip(segs,out)) if _reversed(s,o)]
        if not bad: return out
        segs=[s for i,s in enumerate(segs) if i not in bad[:1]]
    return out
def _reversed(s,o):
    if s[0]=='L':
        return (o[2]-o[1])@(s[2]-s[1])<=0
    if abs(o[5])<1e-9: return True
    if np.sign(o[5])!=np.sign(s[5]): return True
    # an arc whose chord now points against the original chord has turned inside out
    c0=s[2]-s[1]; c1=o[2]-o[1]
    return np.hypot(*c0)>1e-6 and abs(s[5])<np.pi and c0@c1<0
def _offset_once(segs,d):
    n=len(segs); C=[_moved(s,d) for s in segs]
    J=[]
    for i in range(n):
        a,b=segs[i-1],segs[i]; P=a[2]                         # joint between segment i-1 and i
        ua,ub=tangent(a,1),tangent(b,0)
        na=np.array([-ua[1],ua[0]]); nb=np.array([-ub[1],ub[0]])
        smooth=np.degrees(np.arccos(np.clip(ua@ub,-1,1)))<0.5
        if smooth:
            J.append(P+(na+nb)/np.linalg.norm(na+nb)*d); continue
        # a corner: intersect the moved carriers near the mitre point
        m=na+nb; g=P+m/max(np.linalg.norm(m),1e-9)*d/max(np.cos(np.arccos(np.clip(na@nb,-1,1))/2),0.2)
        q=_intersect(C[i-1],C[i],g)
        J.append(q if q is not None else g)
    out=[]
    for i in range(n):
        p0,p1=J[i],J[(i+1)%n]; s=segs[i]
        if s[0]=='L': out.append(('L',p0,p1)); continue
        e=C[i][1]; t0=_par(e,p0); t1=_par(e,p1); k=np.sign(s[5])
        if n==1: out.append(('E',p0,p1,e,t0,s[5])); continue      # a whole ellipse (a dot)
        sw=(t1-t0)%(2*np.pi) if k>0 else -((t0-t1)%(2*np.pi))
        if abs(sw)>abs(s[5])+np.pi: sw=0.0                    # the arc vanished (crossed over itself)
        out.append(('E',p0,p1,e,t0,sw))
    return out
def clean(segs,minlen=0.05):
    """drop segments that collapsed to (nearly) nothing"""
    return [s for s in segs if np.hypot(*(s[2]-s[1]))>minlen or (s[0]=='E' and abs(s[5])>1e-3)]
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
            else:
                for bz in lcoutline.beziers(s): pen.curveTo(T(bz[1]),T(bz[2]),T(bz[3]))
        pen.closePath()
def to_polygons(contours,n=16):
    """flattened rings (for checks and overlays)"""
    rings=[]
    for segs in contours:
        pts=[]
        for s in segs:
            if s[0]=='L': pts.append(s[1])
            else:
                for t in np.linspace(0,1,n,endpoint=False): pts.append(_ept(s[3],s[4]+s[5]*t))
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
