"""Series B, page B-1 (A B C D E F).  Same constructions as Series C (dimC1)."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimC1 import _counter_end
R8={'A':'1','B':'4 1/4','C':'2 1/8','D':'3 3/8','E':'5 3/4','F':'2 3/16','G':'1 1/8','H':'1 1/2','I':'1/2','J':'1/8',
    'K':'16','L':'1 9/16','M':'2','N':'1 3/4','O':'2 15/16','P':'3 7/8','Q':'3','R':'1 1/4'}
R12={'A':'1 1/2','B':'6 3/8','C':'3 3/16','D':'5 1/16','E':'8 5/8','F':'3 9/32','G':'1 11/16','H':'2 1/4','I':'3/4',
     'J':'3/16','K':'24','L':'2 11/32','M':'3','N':'2 5/8','O':'4 13/32','P':'5 13/16','Q':'4 1/2','R':'1 7/8'}
R18={'A':'2 1/4','B':'9 9/16','C':'4 25/32','D':'7 19/32','E':'12 15/16','F':'4 15/16','G':'2 17/32','H':'3 3/8',
     'I':'1 1/8','J':'9/32','K':'36','L':'3 1/2','M':'4 1/2','N':'3 15/16','O':'6 5/8','P':'8 23/32','Q':'6 3/4',
     'R':'2 13/16'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
def letter_A():
    g=Grid(); B,C=DIM['B'],DIM['C']; xl,xr=(B-A)/2,(B+A)/2
    def edge(p0,p1): d=np.array(p1,float)-np.array(p0,float); return np.array(p0,float),d/np.linalg.norm(d)
    pL,dL=edge((0,0),(xl,CAP)); pR,dR=edge((B,0),(xr,CAP))
    side=lambda p,d: (g.X-p[0])*d[1]-(g.Y-p[1])*d[0]
    outer=(side(pL,dL)>=0)&(side(pR,dR)<=0)&(g.Y>=0)&(g.Y<=CAP)
    nL=np.array([dL[1],-dL[0]]); nR=np.array([-dR[1],dR[0]])
    inner=(side(pL+nL*A,dL)>=0)&(side(pR+nR*A,dR)<=0)
    return outer&~(inner&((g.Y<C)|(g.Y>C+A)))
def letter_B():
    """flat-topped bowls: side arcs E centred C below the cap (upper, J in from the right) and F above the baseline
    (lower), corners H; counters are the bowls offset by A.  Middle bar A thick, centred on the crotch height B."""
    g=Grid(); B,C,D,E,F,H,I,J=[DIM[k] for k in 'BCDEFHIJ']
    cu=np.array([D-J-E,CAP-C]); cl=np.array([D-E,F])
    ut,ub=CAP-A,B+A/2; lt,lb=B-A/2,A
    up_cnt,_,_=_counter_end(g,A,ut,ub,cu,E-A,H-A)
    lo_cnt,_,_=_counter_end(g,A,lt,lb,cl,E-A,H-A)
    up_out,_,_=_counter_end(g,0,CAP,ub-A,cu,E,H)
    lo_out,_,_=_counter_end(g,0,lt+A,0,cl,E,H)
    return (up_out|lo_out)&~up_cnt&~lo_cnt,{'counter corner H-A vs I':(round(H-A,1),round(I,1))}
def letter_C():
    """big arc K centred at mid-height (left side at x=0); arcs L tangent to it and to the overshoot lines, ending at
    their rightmost points with horizontal cuts; the lower terminal carried straight up to M above the overshoot line"""
    g=Grid(); D,K,L,M=[DIM[k] for k in 'DKLM']
    cK=np.array([K,CAP/2]); yu=TOP-L; yl=BOT+L
    xc=cK[0]-np.sqrt((K-L)**2-(yu-cK[1])**2); cu=np.array([xc,yu]); cl=np.array([xc,yl])
    m=ring(g,cK,K,K-A)&g.sector(cK,ang(cK,cu),ang(cK,cl))
    m|=ring(g,cu,L,L-A)&g.sector(cu,0,ang(cK,cu))
    m|=ring(g,cl,L,L-A)&g.sector(cl,ang(cK,cl),0)
    m|=rect(g,xc+L-A,xc+L,yl,BOT+M)
    return m,{'tip x':round(xc+L,1),'D':round(D,1),'lower straight':round(BOT+M-yl,1)}
def letter_D():
    """big arc K centred at mid-height (D-K from the left), corner arcs N tangent to it and to the top/bottom lines"""
    g=Grid(); D,K,N=DIM['D'],DIM['K'],DIM['N']
    cU=np.array([D-K,CAP/2]); yu=CAP-N; yl=N
    xT=cU[0]+np.sqrt((K-N)**2-(yu-cU[1])**2); cu=np.array([xT,yu]); cl=np.array([xT,yl])
    a_u,a_l=ang(cU,cu),ang(cU,cl)
    def body(d):
        r=rect(g,d,xT,d,CAP-d)|(g.disk(cU,K-d)&g.sector(cU,a_l,a_u)&(g.X>=xT))
        r|=g.disk(cu,N-d)&(g.X>=xT)&g.sector(cu,a_u,90)
        r|=g.disk(cl,N-d)&(g.X>=xT)&g.sector(cl,-90,a_l)
        return r
    return body(0)&~body(A)
def letter_E():
    g=Grid(); O,P,Q,R=DIM['O'],DIM['P'],DIM['Q'],DIM['R']
    return rect(g,0,A,0,CAP)|rect(g,0,O,CAP-A,CAP)|rect(g,0,A+R,P,P+A)|rect(g,0,Q,0,A)
def letter_F():
    g=Grid(); G,P,Q=DIM['G'],DIM['P'],DIM['Q']
    return rect(g,0,A,0,CAP)|rect(g,0,Q,CAP-A,CAP)|rect(g,0,A+G,P,P+A)
