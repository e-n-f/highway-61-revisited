"""Fitting helpers for the lower-case constructions: straight edges and axis-aligned ellipses, fitted robustly to the
outline points of a letter (font units, from lcview.units); results also in squares (SQ units)."""
import numpy as np
from scipy.optimize import least_squares
from lcview import units, SQ
def pts(ch,pred):
    P=np.vstack(units(ch)); return P[pred(P[:,0],P[:,1])]
def ell_dist(P,cx,cy,a,b):
    """approximate signed distance to an axis-aligned ellipse (positive outside)"""
    u=(P[:,0]-cx)/a; v=(P[:,1]-cy)/b; r=np.hypot(u,v)
    g=np.hypot(u/a,v/b)/np.maximum(r,1e-9)   # gradient magnitude of r
    return (r-1)/np.maximum(g,1e-9)
def fit_ellipse(P,init=None,fix=None,trim=3.0):
    """least-squares axis-aligned ellipse; fix: dict of parameters held (cx, cy, a, b); points beyond `trim` units of the
    fit are dropped and the fit repeated (grid-line stubs)"""
    names=['cx','cy','a','b']; fix=fix or {}
    if init is None: init=alg_ellipse(P)
    free=[n for n in names if n not in fix]
    def full(v):
        d=dict(fix); d.update(zip(free,v)); return [d[n] for n in names]
    x0=[init[names.index(n)] for n in free]; Q=P
    for _ in range(4):
        lo=[-np.inf if n in ('cx','cy') else 1.0 for n in free]
        r=least_squares(lambda v: ell_dist(Q,*full(v)),x0,loss='soft_l1',f_scale=1.0,bounds=(lo,np.inf)); x0=r.x
        d=np.abs(ell_dist(P,*full(x0))); Q=P[d<trim]
    p=full(x0); d=ell_dist(Q,*p)
    return dict(zip(names,p)),np.sqrt(np.mean(d**2)),np.abs(d).max(),len(Q)
def alg_ellipse(P):
    """algebraic axis-aligned conic fit A x^2 + C y^2 + D x + E y = 1 (centred coordinates)"""
    m=P.mean(0); X,Y=(P-m).T
    A,C,D,E=np.linalg.lstsq(np.c_[X*X,Y*Y,X,Y],np.ones(len(P)),rcond=None)[0]
    if A<=0 or C<=0: return [m[0],m[1],np.ptp(P[:,0])/2+1,np.ptp(P[:,1])/2+1]
    x0,y0=-D/(2*A),-E/(2*C); k=1+A*x0*x0+C*y0*y0
    return [m[0]+x0,m[1]+y0,np.sqrt(k/A),np.sqrt(k/C)]
def fit_circle(P,init=None,trim=3.0):
    Q=P; x=init if init is not None else alg_ellipse(P)[:3]
    for _ in range(4):
        f=lambda v,Q=Q: np.hypot(Q[:,0]-v[0],Q[:,1]-v[1])-v[2]
        x=least_squares(f,x,loss='soft_l1').x; d=np.abs(np.hypot(P[:,0]-x[0],P[:,1]-x[1])-x[2]); Q=P[d<trim]
    d=np.hypot(Q[:,0]-x[0],Q[:,1]-x[1])-x[2]
    return dict(cx=x[0],cy=x[1],a=x[2],b=x[2]),np.sqrt(np.mean(d**2)),np.abs(d).max(),len(Q)
def fit_line(P,trim=2.0):
    """total-least-squares line: returns point, unit direction, rms"""
    Q=P
    for _ in range(4):
        c=Q.mean(0); _,_,vt=np.linalg.svd(Q-c); u=vt[0]; n=np.array([-u[1],u[0]])
        d=(P-c)@n; Q=P[np.abs(d-np.median(d[np.abs(d)<10]) if (np.abs(d)<10).any() else d)<trim]
    c=Q.mean(0); _,_,vt=np.linalg.svd(Q-c); u=vt[0]; n=np.array([-u[1],u[0]]); d=(Q-c)@n
    return c,u,np.sqrt(np.mean(d**2)),len(Q)
def sq(v): return np.round(np.asarray(v,float)/SQ,3)
def show(name,e,rms,mx,n):
    print(f"{name:10s} cx {e['cx']:7.1f} cy {e['cy']:7.1f} a {e['a']:6.1f} b {e['b']:6.1f}  | sq: cx {e['cx']/SQ:.3f} cy {e['cy']/SQ:.3f} "
          f"a {e['a']/SQ:.3f} b {e['b']/SQ:.3f}  rms {rms:.2f} max {mx:.2f} n {n}")
def fit_quad(P,init,fix=None,trim=3.0):
    """fit a quad oval (lcgeom.Quad) to points; init/fix: dicts over L yl R yr xt T xb B"""
    from lcgeom import Quad
    names=['L','yl','R','yr','xt','T','xb','B']; fix=fix or {}; free=[n for n in names if n not in fix]
    def full(v):
        d=dict(fix); d.update(zip(free,v)); return Quad(*[d[n] for n in names]),d
    x0=[init[n] for n in free]; Q=P
    for _ in range(4):
        r=least_squares(lambda v: full(v)[0].dist(Q),x0,loss='soft_l1',f_scale=1.0); x0=r.x
        dd=np.abs(full(x0)[0].dist(P)); Q=P[dd<trim]
    q,d=full(x0); res=q.dist(Q)
    return d,np.sqrt(np.mean(res**2)),np.abs(res).max(),len(Q),len(P)
def showq(name,d,rms,mx,n,N):
    print(f"{name:12s} "+' '.join(f"{k} {v:6.1f}({v/SQ:.3f})" for k,v in d.items())+f"  rms {rms:.2f} max {mx:.2f} {n}/{N}")
