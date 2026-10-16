"""Primitives for the lower-case constructions.

The curved parts are quarter-ellipses with horizontal and vertical axes, joined at their axis points so that neighbours
share a tangent there: a "quad oval" is fixed by its four extreme points, left (L, yl), right (R, yr), top (xt, T) and
bottom (xb, B); its quarter centres are (xt, yl) top-left, (xt, yr) top-right, (xb, yl) bottom-left and (xb, yr)
bottom-right.  An ordinary ellipse is the case xt = xb, yl = yr.  Geometry is built as shapely polygons sampled finely
enough (0.25 degree of parameter) that the outline recovery in lcoutline.py can find the exact arcs again."""
import numpy as np
from shapely.geometry import Polygon, box
N=360          # samples per quarter
REG=[]         # every quarter-ellipse made since the last reset (for lcoutline)
def quarter(cx,cy,a,b,q,n=N):
    """points of one quarter, q in 'TR','TL','BL','BR', counter-clockwise from its first axis point"""
    t0={'TR':0,'TL':90,'BL':180,'BR':270}[q]
    t=np.radians(np.linspace(t0,t0+90,n+1)); return np.c_[cx+a*np.cos(t),cy+b*np.sin(t)]
class Quad:
    def __init__(self,L,yl,R,yr,xt,T,xb,B):
        self.L,self.yl,self.R,self.yr,self.xt,self.T,self.xb,self.B=L,yl,R,yr,xt,T,xb,B
    def quarters(self):
        L,yl,R,yr,xt,T,xb,B=self.L,self.yl,self.R,self.yr,self.xt,self.T,self.xb,self.B
        return {'TR':(xt,yr,R-xt,T-yr),'TL':(xt,yl,xt-L,T-yl),'BL':(xb,yl,xb-L,yl-B),'BR':(xb,yr,R-xb,yr-B)}
    def poly(self):
        Q=self.quarters(); pts=[quarter(*Q[q],q)[:-1] for q in ('TR','TL','BL','BR')]
        REG.extend(Q.values())
        P=Polygon(np.vstack(pts))
        return P if P.is_valid else P.buffer(0)
    def ellipses(self):
        return list(self.quarters().values())
    def dist(self,P):
        """signed distance (approx., positive outside) of points to the oval: each point measured against the quarter
        whose quadrant contains it"""
        Q=self.quarters(); d=np.empty(len(P)); x,y=P[:,0],P[:,1]
        up_l=y>=self.yl; up_r=y>=self.yr; lt=x<=self.xt; lb=x<=self.xb
        for q,m in (('TR',(~lt)&up_r),('TL',lt&up_l),('BL',lb&~up_l),('BR',(~lb)&~up_r)):
            cx,cy,a,b=Q[q]; d[m]=ell_dist(P[m],cx,cy,a,b)
        rest=np.isnan(d)|~(((~lt)&up_r)|(lt&up_l)|(lb&~up_l)|((~lb)&~up_r))
        if rest.any(): d[rest]=np.min([ell_dist(P[rest],*Q[q]) for q in Q],axis=0)
        return d
def ell_dist(P,cx,cy,a,b):
    u=(P[:,0]-cx)/a; v=(P[:,1]-cy)/b; r=np.hypot(u,v)
    g=np.hypot(u/a,v/b)/np.maximum(r,1e-9)
    return (r-1)/np.maximum(g,1e-9)
def ellipse(cx,cy,a,b):
    return Quad(cx-a,cy,cx+a,cy,cx,cy+b,cx,cy-b).poly()
def half(p,u,keep=1):
    """half-plane to the left (keep=1) or right (keep=-1) of the directed line through p along u"""
    u=np.asarray(u,float)/np.hypot(*u); n=np.array([-u[1],u[0]])*keep; p=np.asarray(p,float); L=3000
    return Polygon([p-u*L,p+u*L,p+u*L+n*L,p-u*L+n*L])
def slope_dir(dx,dy):
    v=np.array([dx,dy],float); return v/np.hypot(*v)
