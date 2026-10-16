"""Series C, page C-1 (A B C D E F)."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
R8={'A':'1 1/8','B':'5','C':'2 3/16','D':'4 3/8','E':'2 1/32','F':'1 3/4','G':'2 3/4','H':'2 7/8','I':'1/8','J':'4 3/16',
    'K':'5/8','L':'1','M':'10','N':'2 1/4','O':'3 7/8','P':'4','Q':'3 3/4','R':'1 11/16','S':'1 5/8'}
R12={'A':'1 11/16','B':'7 1/2','C':'3 9/32','D':'6 9/16','E':'3 3/64','F':'2 5/8','G':'4 1/8','H':'4 5/16','I':'3/16',
     'J':'6 9/32','K':'15/16','L':'1 1/2','M':'15','N':'3 3/8','O':'5 13/16','P':'6','Q':'5 5/8','R':'2 17/32','S':'2 7/16'}
R18={'A':'2 17/32','B':'11 1/4','C':'4 15/16','D':'9 27/32','E':'4 9/16','F':'3 15/16','G':'6 3/16','H':'6 15/32',
     'I':'9/32','J':'9 7/16','K':'1 13/32','L':'2 1/4','M':'22 1/2','N':'5 1/16','O':'8 23/32','P':'9','Q':'8 7/16',
     'R':'3 13/16','S':'3 21/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['I']; TOP=CAP+OV; BOT=-OV
def letter_A():
    import dimD1
    g=Grid(); B,C=DIM['B'],DIM['C']; xl,xr=(B-A)/2,(B+A)/2
    def edge(p0,p1): d=np.array(p1,float)-np.array(p0,float); return np.array(p0,float),d/np.linalg.norm(d)
    pL,dL=edge((0,0),(xl,CAP)); pR,dR=edge((B,0),(xr,CAP))
    side=lambda p,d: (g.X-p[0])*d[1]-(g.Y-p[1])*d[0]
    outer=(side(pL,dL)>=0)&(side(pR,dR)<=0)&(g.Y>=0)&(g.Y<=CAP)
    nL=np.array([dL[1],-dL[0]]); nR=np.array([-dR[1],dR[0]])
    inner=(side(pL+nL*A,dL)>=0)&(side(pR+nR*A,dR)<=0)
    return outer&~(inner&((g.Y<C)|(g.Y>C+A)))
def _counter_end(g,x0,top,bot,c,R,k):
    """flat-topped counter from x0 to the right: side arc R centred at c, corners k tangent to the flats and to R"""
    def corner(y,sgn):                      # corner centre: k from the flat, R-k from c
        yc=y-sgn*k; return np.array([c[0]+np.sqrt((R-k)**2-(yc-c[1])**2),yc])
    ct,cb=corner(top,1),corner(bot,-1)
    m=rect(g,x0,2000,bot,top)&(g.disk(c,R)|(g.X<=c[0]))
    # round the corners: beyond each corner arc, in the wedge between the flat and the tangent point with R
    m&=~(g.sector(ct,ang(c,ct),90)&~g.disk(ct,k))
    m&=~(g.sector(cb,-90,ang(c,cb))&~g.disk(cb,k))
    return m,ct,cb
def letter_B():
    """counters: flat top and bottom, side arcs G-A (upper, centred E below the cap) and H-A (lower, centred C above the
    baseline) with corners K; the outer bowls are the counters offset by A (side arcs G, H; corners F = K+A).  The
    middle bar is A thick, centred on the crotch height J.  Upper bowl I in from the right."""
    g=Grid(); C,D,E,F,G,H,I,J,K,L=[DIM[k] for k in 'CDEFGHIJKL']
    cg=np.array([D-I-G,CAP-E]); ch=np.array([D-H,C])
    ut,ub=CAP-A,J+A/2; lt,lb=J-A/2,A
    up_cnt,_,_=_counter_end(g,A,ut,ub,cg,G-A,K)
    lo_cnt,_,_=_counter_end(g,A,lt,lb,ch,H-A,K)
    up_out,_,_=_counter_end(g,0,CAP,ub-A,cg,G,F)
    lo_out,_,_=_counter_end(g,0,lt+A,0,ch,H,F)
    m=(up_out|lo_out)&~up_cnt&~lo_cnt
    return m,{'F vs K+A':(round(F,1),round(K+A,1))}
def letter_C():
    """big arc M centred at mid-height (left side at x=0); arcs E tangent to it and to the overshoot lines, ending at
    their rightmost points with horizontal cuts; the lower terminal is carried straight up to N above the overshoot line"""
    g=Grid(); D,E,M,N=[DIM[k] for k in 'DEMN']
    cM=np.array([M,CAP/2]); yu=TOP-E; yl=BOT+E
    xc=cM[0]-np.sqrt((M-E)**2-(yu-cM[1])**2); cu=np.array([xc,yu]); cl=np.array([xc,yl])
    m=ring(g,cM,M,M-A)&g.sector(cM,ang(cM,cu),ang(cM,cl))
    m|=ring(g,cu,E,E-A)&g.sector(cu,0,ang(cM,cu))
    m|=ring(g,cl,E,E-A)&g.sector(cl,ang(cM,cl),0)
    m|=rect(g,xc+E-A,xc+E,yl,BOT+N)
    return m,{'tip x':round(xc+E,1),'D':round(D,1),'lower straight':round(BOT+N-yl,1)}
def letter_D():
    """big arc M centred at mid-height (D-M from the left), corner arcs N tangent to it and to the top/bottom lines"""
    g=Grid(); D,M,N=DIM['D'],DIM['M'],DIM['N']
    cU=np.array([D-M,CAP/2]); yu=CAP-N; yl=N
    xT=cU[0]+np.sqrt((M-N)**2-(yu-cU[1])**2); cu=np.array([xT,yu]); cl=np.array([xT,yl])
    a_u,a_l=ang(cU,cu),ang(cU,cl)
    def body(d):
        r=rect(g,d,xT,d,CAP-d)|(g.disk(cU,M-d)&g.sector(cU,a_l,a_u)&(g.X>=xT))
        r|=g.disk(cu,N-d)&(g.X>=xT)&g.sector(cu,a_u,90)
        r|=g.disk(cl,N-d)&(g.X>=xT)&g.sector(cl,-90,a_l)
        return r
    return body(0)&~body(A)
def letter_E():
    g=Grid(); O,P,Q,R=DIM['O'],DIM['P'],DIM['Q'],DIM['R']
    return rect(g,0,A,0,CAP)|rect(g,0,O,CAP-A,CAP)|rect(g,0,A+R,Q,Q+A)|rect(g,0,P,0,A)
def letter_F():
    g=Grid(); P,Q,S=DIM['P'],DIM['Q'],DIM['S']
    return rect(g,0,A,0,CAP)|rect(g,0,P,CAP-A,CAP)|rect(g,0,A+S,Q,Q+A)
