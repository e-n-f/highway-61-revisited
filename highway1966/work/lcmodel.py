"""Fit a parametric construction (shapely geometry from a parameter vector) to a letter's scan outline.

The loss is two-sided: distances from the scan's outline points to the model's boundary and from the model's boundary
points to the scan's outline, so neither can leave a piece unexplained.  Points within `stub` units of a grid line
crossing are not special-cased; a soft-L1 loss with a 1.5-unit scale keeps the stubs from pulling."""
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, Point, box
from shapely import affinity
from lcview import units, SQ
SPREAD=1.5          # printing gain per side, from the flat baselines and x-height lines of n m r u v w z
def ellipse(cx,cy,a,b,n=720):
    t=np.linspace(0,2*np.pi,n,endpoint=False); return Polygon(np.c_[cx+a*np.cos(t),cy+b*np.sin(t)])
def boundary_pts(geom,step=2.0):
    out=[]
    polys=getattr(geom,'geoms',[geom])
    for p in polys:
        if p.is_empty or p.geom_type!='Polygon': continue
        for r in [p.exterior,*p.interiors]:
            L=r.length; n=max(int(L/step),8); out.append(np.array([r.interpolate(t).coords[0] for t in np.linspace(0,L,n,endpoint=False)]))
    return np.vstack(out) if out else np.zeros((0,2))
def scan_pts(ch,step=2.0):
    cs=units(ch); out=[]
    for c in cs:
        d=np.r_[0,np.cumsum(np.hypot(*np.diff(c,axis=0).T))]; t=np.arange(0,d[-1],step)
        out.append(np.c_[np.interp(t,d,c[:,0]),np.interp(t,d,c[:,1])])
    return np.vstack(out)
def fit(ch,build,p0,names,fixed=None,region=None,verbose=True,ftol=1e-6):
    fixed=fixed or {}; S=scan_pts(ch)
    if region is not None: S=S[region(S[:,0],S[:,1])]
    tS=cKDTree(S); free=[n for n in names if n not in fixed]
    def par(v):
        d=dict(fixed); d.update(zip(free,v)); return d
    def res(v):
        g=build(par(v)).buffer(SPREAD,join_style=2); M=boundary_pts(g)
        if region is not None: M=M[region(M[:,0],M[:,1])]
        if len(M)==0: return np.full(len(S)+NM,1e3)
        M=M[np.linspace(0,len(M)-1,NM).astype(int)]
        dS,_=cKDTree(M).query(S); dM,_=tS.query(M)
        return np.r_[dS,dM*np.sqrt(len(S)/NM)]
    x0=[p0[n] for n in free]; NM=len(S)//2
    r=least_squares(res,x0,loss='soft_l1',f_scale=1.5,diff_step=1e-3,ftol=ftol,xtol=ftol)
    p=par(r.x); g=build(p); M=boundary_pts(g.buffer(SPREAD,join_style=2))
    if region is not None: M=M[region(M[:,0],M[:,1])]
    dS,_=cKDTree(M).query(S); dM,_=tS.query(M)
    st=dict(p95=np.percentile(np.r_[dS,dM],95),med=np.median(np.r_[dS,dM]),max_scan=dS.max(),max_model=dM.max())
    if verbose:
        print(ch,' '.join(f'{n}={p[n]:.1f}({p[n]/SQ:.3f})' for n in names))
        print('   fit: median %.2f p95 %.2f  max scan->model %.1f  model->scan %.1f'%(st['med'],st['p95'],st['max_scan'],st['max_model']))
    return p,g,st
def overlay(ch,g,out,scale=1.2,extra=()):
    from lcview import view
    import cv2
    M=[np.array(r.coords) for p in getattr(g.buffer(SPREAD,join_style=2),'geoms',[g.buffer(SPREAD,join_style=2)]) if p.geom_type=='Polygon' for r in [p.exterior,*p.interiors]]
    view(ch,extra=[(m,(0,0,255)) for m in M]+list(extra),scale=scale,out=out)
