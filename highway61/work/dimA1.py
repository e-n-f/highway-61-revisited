"""Series A, page A-1 (A B C D E F).  The constructions of Series B (dimB1) with this page's letters."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimC1 import _counter_end
R8={'A':'3/4','B':'3 1/2','C':'2','D':'5','E':'4 1/8','F':'2 1/8','G':'2 1/4','H':'1 1/2','I':'1/16','J':'9/16',
    'K':'1 1/8','L':'3 3/4','M':'3','N':'1 3/8','O':'1 3/4','P':'18','Q':'1/8','R':'2 5/8','S':'1 1/4','T':'2 11/16'}
R12={'A':'1 1/8','B':'5 1/4','C':'3','D':'7 1/2','E':'6 3/16','F':'3 3/16','G':'3 3/8','H':'2 1/4','I':'3/32',
     'J':'27/32','K':'1 11/16','L':'5 5/8','M':'4 1/2','N':'2 1/16','O':'2 5/8','P':'27','Q':'3/16','R':'3 15/16',
     'S':'1 7/8','T':'4 1/32'}
R18={'A':'1 11/16','B':'7 7/8','C':'4 1/2','D':'11 1/4','E':'9 9/32','F':'4 25/32','G':'5 1/16','H':'3 3/8',
     'I':'1/8','J':'1 1/4','K':'2 17/32','L':'8 7/16','M':'6 3/4','N':'3 3/32','O':'3 15/16','P':'40 1/2','Q':'9/32',
     'R':'5 29/32','S':'2 13/16','T':'6 1/16'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['Q']; TOP=CAP+OV; BOT=-OV
def letter_A():
    """apex flat A wide, feet B apart; bar A thick with its lower edge C above the baseline"""
    g=Grid(); B,C=DIM['B'],DIM['C']; xl,xr=(B-A)/2,(B+A)/2
    def edge(p0,p1): d=np.array(p1,float)-np.array(p0,float); return np.array(p0,float),d/np.linalg.norm(d)
    pL,dL=edge((0,0),(xl,CAP)); pR,dR=edge((B,0),(xr,CAP))
    side=lambda p,d: (g.X-p[0])*d[1]-(g.Y-p[1])*d[0]
    outer=(side(pL,dL)>=0)&(side(pR,dR)<=0)&(g.Y>=0)&(g.Y<=CAP)
    nL=np.array([dL[1],-dL[0]]); nR=np.array([-dR[1],dR[0]])
    inner=(side(pL+nL*A,dL)>=0)&(side(pR+nR*A,dR)<=0)
    return outer&~(inner&((g.Y<C)|(g.Y>C+A)))
def letter_B():
    """flat-topped bowls: side arcs D centred F below the cap (upper, I in from the right) and E centred G above the
    baseline (lower), corners H; counters are the bowls offset by A (corners H-A, which equals A).  The middle bar is
    A thick, centred on the crotch height E."""
    g=Grid(); D,E,F,G,H,I,J,M=[DIM[k] for k in 'DEFGHIJM']
    cu=np.array([M-I-D,CAP-F]); cl=np.array([M-E,G])
    ut,ub=CAP-A,E+A/2; lt,lb=E-A/2,A
    up_cnt,_,_=_counter_end(g,A,ut,ub,cu,D-A,H-A)
    lo_cnt,_,_=_counter_end(g,A,lt,lb,cl,E-A,H-A)
    up_out,_,_=_counter_end(g,0,CAP,ub-A,cu,D,H)
    lo_out,_,_=_counter_end(g,0,lt+A,0,cl,E,H)
    m=(up_out|lo_out)&~up_cnt&~lo_cnt
    return m,{'counter corner H-A vs A':(round(H-A,1),round(A,1)),'J':round(J,1)}
def letter_C():
    """big arc P centred at mid-height (left side at x=0); arcs N tangent to it and to the overshoot lines, ending at
    their rightmost points; the terminals run straight on to O below the upper overshoot line and C above the lower
    one, cut horizontally"""
    g=Grid(); C,M,N,O,P=[DIM[k] for k in 'CMNOP']
    cP=np.array([P,CAP/2]); yu=TOP-N; yl=BOT+N
    xc=cP[0]-np.sqrt((P-N)**2-(yu-cP[1])**2); cu=np.array([xc,yu]); cl=np.array([xc,yl])
    m=ring(g,cP,P,P-A)&g.sector(cP,ang(cP,cu),ang(cP,cl))
    m|=ring(g,cu,N,N-A)&g.sector(cu,0,ang(cP,cu))
    m|=ring(g,cl,N,N-A)&g.sector(cl,ang(cP,cl),0)
    m|=rect(g,xc+N-A,xc+N,TOP-O,yu)|rect(g,xc+N-A,xc+N,yl,BOT+C)
    return m,{'tip x':round(xc+N,1),'M':round(M,1)}
def letter_D():
    """big arc P centred at mid-height (M-P from the left), corner arcs H tangent to it and to the top/bottom lines"""
    g=Grid(); M,P,H=DIM['M'],DIM['P'],DIM['H']
    cU=np.array([M-P,CAP/2]); yu=CAP-H; yl=H
    xT=cU[0]+np.sqrt((P-H)**2-(yu-cU[1])**2); cu=np.array([xT,yu]); cl=np.array([xT,yl])
    a_u,a_l=ang(cU,cu),ang(cU,cl)
    def body(d):
        r=rect(g,d,xT,d,CAP-d)|(g.disk(cU,P-d)&g.sector(cU,a_l,a_u)&(g.X>=xT))
        r|=g.disk(cu,H-d)&(g.X>=xT)&g.sector(cu,a_u-0.5,90)       # sectors overlap the big arc's slightly: no seam
        r|=g.disk(cl,H-d)&(g.X>=xT)&g.sector(cl,-90,a_l+0.5)
        return r
    return body(0)&~body(A)
def letter_E():
    """top bar R, bottom bar T; middle bar S beyond the stem, its lower edge L above the baseline"""
    g=Grid(); L,R,S,T=DIM['L'],DIM['R'],DIM['S'],DIM['T']
    return rect(g,0,A,0,CAP)|rect(g,0,R,CAP-A,CAP)|rect(g,0,A+S,L,L+A)|rect(g,0,T,0,A)
def letter_F():
    """top bar T; middle bar K beyond the stem, its lower edge L above the baseline"""
    g=Grid(); K,L,T=DIM['K'],DIM['L'],DIM['T']
    return rect(g,0,A,0,CAP)|rect(g,0,T,CAP-A,CAP)|rect(g,0,A+K,L,L+A)
