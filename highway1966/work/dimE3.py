"""Series E, page E-3 (M N O P Q R)."""
import numpy as np
from fractions import Fraction as Fr
from dim1961 import Grid
from scipy.optimize import brentq
def inch(s): return float(sum(Fr(p) for p in s.split()))
R8={'A':'1 3/8','B':'7 3/8','C':'4 7/8','D':'6 3/8','E':'2 1/2','F':'1 1/4','G':'6 5/8','H':'3 3/16','I':'2 7/16','J':'1/8',
    'K':'2 1/4','L':'3/16','M':'2 1/8','N':'3 1/8','O':'1 5/8','P':'3/8'}
DIM={k:inch(v)*690/8 for k,v in R8.items()}
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
ring=lambda g,c,ro,ri: g.disk(c,ro)&~g.disk(c,ri)
ang=lambda c,p: np.degrees(np.arctan2(p[1]-c[1],p[0]-c[0]))
def band(g,p,u,w):
    """points between the line through p (direction u) and its parallel w away on the left of u"""
    n=np.array([-u[1],u[0]]); d=(g.X-p[0])*n[0]+(g.Y-p[1])*n[1]
    return (d>=0)&(d<=w) if w>0 else (d<=0)&(d>=w)
def strip(g,p,u,w,toward):
    """band of width w along the line (p,u), on the side containing the point 'toward'"""
    n=np.array([-u[1],u[0]]); sgn=np.sign((toward[0]-p[0])*n[0]+(toward[1]-p[1])*n[1]) or 1
    d=((g.X-p[0])*n[0]+(g.Y-p[1])*n[1])*sgn
    return (d>=0)&(d<=w)
def cross_solve(p,q,w,lo,hi):
    """direction angle of parallel lines through p and q that are w apart"""
    f=lambda t: abs((q-p)[0]*np.sin(t)-(q-p)[1]*np.cos(t))-w
    ts=np.linspace(lo,hi,400); vs=[f(t) for t in ts]
    for i in range(len(ts)-1):
        if vs[i]*vs[i+1]<=0: return brentq(f,ts[i],ts[i+1])
    return None
def letter_M():
    g=Grid(); B,C=DIM['B'],DIM['C']
    stems=((g.X>=0)&(g.X<=A)|(g.X>=B-A)&(g.X<=B))&(g.Y>=0)&(g.Y<=CAP)
    v=np.array([B/2,0.0]); pr=np.array([B-A,C]); u=(pr-v)/np.linalg.norm(pr-v)     # right diagonal's lower edge
    right=strip(g,v,u,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X<=B-A+1)&(g.X>=B/2)
    vl=np.array([B/2,0.0]); ul=np.array([-u[0],u[1]]); left=strip(g,vl,ul,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=B/2)
    return g.bitmap(stems|right|left)
def letter_N():
    g=Grid(); D,E,F=DIM['D'],DIM['E'],DIM['F']
    stems=((g.X>=0)&(g.X<=A)|(g.X>=D-A)&(g.X<=D))&(g.Y>=0)&(g.Y<=CAP)
    p0=np.array([A,CAP-E]); p1=np.array([D-F,0.0]); u=(p1-p0)/np.linalg.norm(p1-p0)  # diagonal's lower edge
    diag=strip(g,p0,u,A,(D,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=D-A+1)
    below=((g.X-p0[0])*u[1]-(g.Y-p0[1])*u[0])>0          # below/left of the diagonal's lower edge
    stems&=~(below&(g.X>D/2))                            # the diagonal runs straight to (D-F, 0), trimming the right stem's foot
    return g.bitmap(stems|diag)
def _O(g):
    W,G=DIM['G'],DIM['G']
    # top/bottom arcs are connecting arcs: the table's H crosses the side arcs (0.44 units over, a 2.3 degree kink);
    # refitted to be tangent with the top kept on the overshoot line (274.3 for 274.9, within the table's 1/16 inch)
    dx=W/2; h=TOP-CAP/2; H=(G*G-dx*dx-h*h)/(2*(G-h))
    cl=np.array([G,CAP/2]); cr=np.array([W-G,CAP/2]); cu=np.array([W/2,TOP-H]); cd=np.array([W/2,BOT+H])
    def area(rs,rt):
        lens=g.disk(cl,rs)&g.disk(cr,rs)
        top=g.sector(cu,ang(cr,cu),ang(cl,cu))&~g.disk(cu,rt)&(g.Y>cu[1])
        bot=g.sector(cd,ang(cl,cd),ang(cr,cd))&~g.disk(cd,rt)&(g.Y<cd[1])
        return lens&~top&~bot
    return area(G,H)&~area(G-A,H-A),W
def letter_O():
    g=Grid(); m,W=_O(g); return g.bitmap(m)
def letter_Q():
    g=Grid(); m,W=_O(g); F,K,L,O,P=DIM['F'],DIM['K'],DIM['L'],DIM['O'],DIM['P']
    TL=np.array([W/2+L,BOT+K]); BL=np.array([W/2+O,BOT-P]); u=(BL-TL)/np.linalg.norm(BL-TL)
    n=np.array([-u[1],u[0]])                      # to the right of the downward axis
    d=(g.X-TL[0])*n[0]+(g.Y-TL[1])*n[1]; t=(g.X-TL[0])*u[0]+(g.Y-TL[1])*u[1]
    tail=(d>=0)&(d<=F)&(t>=0)&(t<=np.linalg.norm(BL-TL))
    return g.bitmap(m|tail)
def _bowl(g,D,I,N,J=0.0):
    c=np.array([D-J-I,(CAP+N)/2])
    outer=((g.X>=0)&(g.X<=c[0])&(g.Y>=N)&(g.Y<=CAP))|(g.disk(c,I)&(g.X>=c[0]))
    inner=((g.X>=A)&(g.X<=c[0])&(g.Y>=N+A)&(g.Y<=CAP-A))|(g.disk(c,I-A)&(g.X>=c[0]))
    return outer&~inner,c
def letter_P():
    g=Grid(); D,I,N=DIM['D'],DIM['I'],DIM['N']
    b,c=_bowl(g,D,I,N); return g.bitmap(b|((g.X>=0)&(g.X<=A)&(g.Y>=0)&(g.Y<=CAP)))
def letter_R():
    g=Grid(); D,I,N,M,J=DIM['D'],DIM['I'],DIM['N'],DIM['M'],DIM['J']
    b,c=_bowl(g,D,I,N,J)
    p=np.array([A+M,N]); q=np.array([D,0.0]); dx,dy=q-p          # leg: upper edge from (A+M, N); lower edge to the corner (D,0); A apart
    th=cross_solve(p,q,A,-np.pi/2+0.01,np.arctan2(dy,dx))
    u=np.array([np.cos(th),np.sin(th)])
    leg=strip(g,p,u,A,q)&(g.Y>=0)&(g.Y<=N+A)
    return g.bitmap(b|leg|((g.X>=0)&(g.X<=A)&(g.Y>=0)&(g.Y<=CAP))),{'leg angle':round(np.degrees(th),1)}
