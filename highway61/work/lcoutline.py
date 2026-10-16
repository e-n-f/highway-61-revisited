"""Exact outlines for the lower-case constructions: every curved edge of a construction lies on one of the
quarter-ellipses recorded in lcgeom.REG while it was built, so each polygon edge is classified as a straight edge or as
a piece of a known ellipse, and runs on the same ellipse are drawn as elliptical arcs (cubic Beziers, at most 45 degrees
of the ellipse's parameter per piece, the affine image of the usual circular-arc approximation)."""
import numpy as np
from shapely.geometry.polygon import orient
import lcgeom
TOL=0.03                     # distance (units) of an edge end from its ellipse: sampled vertices are on it, cut
                             # points lie on a chord whose sagitta is far smaller than this
MAXDT=np.radians(1.0)        # an ellipse edge spans at most 1 degree of parameter (sampling is 0.25 degree)
def _prims():
    seen={}; out=[]
    for e in lcgeom.REG:
        k=tuple(np.round(e,6))
        if k not in seen: seen[k]=1; out.append(np.array(e,float))
    return out
def _on(e,p):
    cx,cy,a,b=e
    if a<=1e-6 or b<=1e-6: return None
    u=(p[0]-cx)/a; v=(p[1]-cy)/b; r=np.hypot(u,v)
    g=np.hypot(u/a,v/b)/max(r,1e-12)
    return abs(r-1)/max(g,1e-12), np.arctan2(v,u)
def _edge_prim(P,q,prims,prev):
    order=([prev] if prev is not None else [])+[i for i in range(len(prims)) if i!=prev]
    for i in order:
        a=_on(prims[i],P); b=_on(prims[i],q)
        if a is None or b is None or a[0]>TOL or b[0]>TOL: continue
        dt=(b[1]-a[1]+np.pi)%(2*np.pi)-np.pi
        if 1e-9<abs(dt)<MAXDT: return i
    return None
def _clean(pts):
    out=[p for i,p in enumerate(pts) if np.linalg.norm(p-pts[i-1])>1e-6]
    return out
def segments(ring,prims):
    pts=_clean([np.array(p,float) for p in list(ring)[:-1]]); n=len(pts)
    lab=[]; prev=None
    for i in range(n):
        k=_edge_prim(pts[i],pts[(i+1)%n],prims,prev); lab.append(k); prev=k
    # start at a label change
    s0=next((i for i in range(n) if lab[i]!=lab[i-1] or lab[i] is None),0)
    P=[pts[(s0+i)%n] for i in range(n)]+[pts[s0]]; L=[lab[(s0+i)%n] for i in range(n)]
    segs=[]; i=0
    while i<n:
        if L[i] is None:
            segs.append(('L',P[i],P[i+1])); i+=1; continue
        j=i
        while j+1<n and L[j+1]==L[i]: j+=1
        e=prims[L[i]]; t=[_on(e,P[k])[1] for k in range(i,j+2)]
        sweep=sum((t[k+1]-t[k]+np.pi)%(2*np.pi)-np.pi for k in range(len(t)-1))
        segs.append(('E',P[i],P[j+1],e,t[0],sweep)); i=j+1
    # merge collinear consecutive lines
    out=[]
    for s in segs:
        if out and s[0]=='L' and out[-1][0]=='L':
            a,b,c=out[-1][1],out[-1][2],s[2]; e1,e2=b-a,c-b
            if abs(e1[0]*e2[1]-e1[1]*e2[0])<1e-6*np.linalg.norm(e1)*np.linalg.norm(e2) and e1@e2>0:
                out[-1]=('L',a,c); continue
        out.append(s)
    return out
def beziers(s):
    _,p0,p1,(cx,cy,a,b),t0,sw=s
    m=max(1,int(np.ceil(abs(sw)/np.radians(45)-1e-9))); d=sw/m; k=4/3*np.tan(d/4)
    E=lambda t: np.array([cx+a*np.cos(t),cy+b*np.sin(t)]); D=lambda t: np.array([-a*np.sin(t),b*np.cos(t)])
    out=[]; start=p0
    for i in range(m):
        ta,tb=t0+i*d,t0+(i+1)*d; end=p1 if i==m-1 else E(tb)
        out.append((start,start+k*D(ta),end-k*D(tb),end)); start=end
    return out
def contours(geom,prims):
    geom=orient(geom,sign=-1.0) if geom.geom_type=='Polygon' else type(geom)([orient(g,sign=-1.0) for g in geom.geoms])
    polys=[geom] if geom.geom_type=='Polygon' else list(geom.geoms)
    return [segments(r.coords,prims) for p in polys for r in [p.exterior,*p.interiors]]
def draw(geom,pen,prims=None,dx=0.0,dy=0.0):
    prims=_prims() if prims is None else prims
    T=lambda p:(float(p[0]+dx),float(p[1]+dy))
    for segs in contours(geom,prims):
        pen.moveTo(T(segs[0][1]))
        for k,s in enumerate(segs):
            if s[0]=='L':
                if k<len(segs)-1: pen.lineTo(T(s[2]))
            else:
                for b in beziers(s): pen.curveTo(T(b[1]),T(b[2]),T(b[3]))
        pen.closePath()
def describe(geom,prims=None):
    prims=_prims() if prims is None else prims
    return [(sum(s[0]=='L' for s in c),sum(s[0]=='E' for s in c)) for c in contours(geom,prims)]
