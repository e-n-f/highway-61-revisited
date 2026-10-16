"""Series E, page E-2 (G H I J K L)."""
import numpy as np
from fractions import Fraction as Fr
from dim1961 import Grid
from scipy.optimize import brentq
def inch(s): return float(sum(Fr(p) for p in s.split()))
R8={'A':'1 3/8','B':'6 3/8','C':'3 1/16','D':'5 7/8','E':'1/8','F':'4 1/8','G':'2','H':'2 3/4','I':'3 1/2','J':'2 1/4','K':'1 1/4',
    'L':'6','M':'2 5/8','N':'1/2','O':'6 1/2'}
R18={'A':'3 3/32','B':'14 11/32','C':'6 57/64','D':'13 7/32','E':'9/32','F':'9 9/32','G':'4 1/2','H':'6 3/16','I':'7 7/8','J':'5 1/16',
     'K':'2 13/16','L':'13 1/2','M':'5 29/32','N':'1 1/8','O':'14 5/8'}
DIM={k:inch(v)*690/18 for k,v in R18.items()}
CHECK={k:round((inch(R18[k])/18-inch(R8[k])/8)*690,2) for k in R18}
A=DIM['A']; CAP=690.0; OV=DIM['E']; TOP=CAP+OV; BOT=-OV
ring=lambda g,c,ro,ri: g.disk(c,ro)&~g.disk(c,ri)
ang=lambda c,p: np.degrees(np.arctan2(p[1]-c[1],p[0]-c[0]))
def letter_G():
    g=Grid(); B,C,D,F,G,H=DIM['B'],DIM['C'],DIM['D'],DIM['F'],DIM['G'],DIM['H']
    cD=np.array([D,CAP/2]); yu=TOP-C; yl=BOT+C
    xc=cD[0]-np.sqrt((D-C)**2-(yu-cD[1])**2); cu=np.array([xc,yu]); cl=np.array([xc,yl])
    yt=TOP-G; xt=cu[0]+np.sqrt(C*C-(yt-cu[1])**2)
    xr=cl[0]+C                                                         # right side: the lower arc's right edge, carried up vertically
    m=ring(g,cD,D,D-A)&g.sector(cD,ang(cD,cu),ang(cD,cl))
    m|=ring(g,cu,C,C-A)&g.sector(cu,ang(cu,(xt,yt)),ang(cD,cu))
    m|=ring(g,cl,C,C-A)&g.sector(cl,ang(cD,cl),0)
    m|=(g.X>=xr-A)&(g.X<=xr)&(g.Y>=cl[1])&(g.Y<=F)                    # jaw: vertical down to the lower arc
    m|=(g.X>=xr-H)&(g.X<=xr)&(g.Y>=F-A)&(g.Y<=F)                      # bar, H long, A thick, top at F
    return g.bitmap(m),{'right side':round(xr,1),'B':round(B,1)}
def letter_H():
    g=Grid(); B,I=DIM['B'],DIM['I']
    m=((g.X>=0)&(g.X<=A)|(g.X>=B-A)&(g.X<=B))&(g.Y>=0)&(g.Y<=CAP); m|=(g.X>=0)&(g.X<=B)&(g.Y>=I)&(g.Y<=I+A)
    return g.bitmap(m)
def letter_I():
    g=Grid(); return g.bitmap((g.X>=0)&(g.X<=A)&(g.Y>=0)&(g.Y<=CAP))
def letter_J():
    g=Grid(); C,J,L=DIM['C'],DIM['J'],DIM['L']
    c=np.array([L-C,BOT+C]); yt=BOT+J; xt=c[0]-np.sqrt(C*C-(yt-c[1])**2)
    m=(g.X>=L-A)&(g.X<=L)&(g.Y>=c[1])&(g.Y<=CAP)
    m|=ring(g,c,C,C-A)&g.sector(c,ang(c,(xt,yt)),0)
    return g.bitmap(m),{'terminal x':round(xt,1)}
def letter_K():
    g=Grid(); O,M,N,J=DIM['O'],DIM['M'],DIM['N'],DIM['J']
    P1=np.array([A,J]); P2=np.array([O-N,CAP])                     # arm's lower edge: stem at height J to the top-right corner
    u=(P2-P1)/np.linalg.norm(P2-P1); nu=np.array([-u[1],u[0]])        # nu points up-left
    yc=CAP-M; crotch=P1+u*((yc-P1[1])/u[1])
    v=crotch-np.array([O,0.0]); v/=np.linalg.norm(v); nv=np.array([-v[1],v[0]])
    if nv[0]>0: nv=-nv
    side=lambda p,d: (g.X-p[0])*d[1]-(g.Y-p[1])*d[0]
    stem=(g.X>=0)&(g.X<=A)&(g.Y>=0)&(g.Y<=CAP)
    arm=(side(P1,u)>=0)&(side(P1+nu*A,u)<=0)&(g.X>=A-1)&(g.Y<=CAP) if side(np.array([0,CAP]),u).mean() else None
    # arm: between its lower edge (through P1,P2) and the parallel edge A further up-left
    s1=side(P1,u); s2=side(P1+nu*A,u)
    arm=(np.sign(s1)!=np.sign(s2))&(g.X>=A-1)&(g.Y<=CAP)
    t1=side(np.array([O,0.0]),v); t2=side(np.array([O,0.0])+nv*A,v)
    leg=(np.sign(t1)!=np.sign(t2))&(g.Y>=0)&(g.Y<=yc)
    return g.bitmap(stem|arm|leg),{'arm angle':round(np.degrees(np.arctan2(u[1],u[0])),2),'crotch':crotch.round(1)}
def letter_L():
    g=Grid(); D,K=DIM['D'],DIM['K']
    return g.bitmap(((g.X>=0)&(g.X<=A)&(g.Y>=0)&(g.Y<=CAP))|((g.X>=0)&(g.X<=D)&(g.Y>=0)&(g.Y<=K)))
