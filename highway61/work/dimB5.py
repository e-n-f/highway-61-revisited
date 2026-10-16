"""Series B, page B-5 (Y Z 1 2 3 4)."""
import numpy as np
from scipy.optimize import brentq, fsolve
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimD5 import line_inter, poly
from dimE4 import arc
R8={'A':'1','B':'4 1/4','C':'3 5/8','D':'3 3/8','E':'3/16','F':'1 1/8','G':'3/4','H':'5/8','I':'1/2','J':'2 7/8',
    'K':'1 11/16','L':'3 1/2','M':'7 1/4','N':'1 19/32','O':'1/8','P':'4 1/2','Q':'4 3/8','R':'1 5/8','S':'1 13/16',
    'T':'2 3/8','U':'3/8','V':'15/16','W':'3 7/8','X':'2'}
R12={'A':'1 1/2','B':'6 3/8','C':'5 7/16','D':'5 1/16','E':'9/32','F':'1 11/16','G':'1 1/8','H':'15/16','I':'3/4',
     'J':'4 5/16','K':'2 17/32','L':'5 1/4','M':'10 7/8','N':'2 25/64','O':'3/16','P':'6 3/4','Q':'6 9/16','R':'2 7/16',
     'S':'2 23/32','T':'3 9/16','U':'9/16','V':'1 13/32','W':'5 13/16','X':'3'}
R18={'A':'2 1/4','B':'9 9/16','C':'8 5/32','D':'7 19/32','E':'7/16','F':'2 17/32','G':'1 11/16','H':'1 13/32',
     'I':'1 1/8','J':'6 15/32','K':'3 13/16','L':'7 7/8','M':'16 5/16','N':'3 37/64','O':'9/32','P':'10 1/8',
     'Q':'9 27/32','R':'3 21/32','S':'4 3/32','T':'5 11/32','U':'27/32','V':'2 1/8','W':'8 23/32','X':'4 1/2'}
# M: the 8-inch entry is cut off in the scan; 7 1/4 agrees with the 12- and 18-inch rows
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['O']; TOP=CAP+OV; BOT=-OV
def letter_Y():
    B,C=DIM['B'],DIM['C']; c=B/2; sl,sr=c-A/2,c+A/2
    p0=np.array([0,CAP]); p1=np.array([sl,C]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*A; xt=line_inter(qi,d,(0,CAP),(1,0))[0]; v=line_inter(qi,d,(c,0),(0,1))
    return poly([(0,CAP),(xt,CAP),tuple(v),(B-xt,CAP),(B,CAP),(sr,C),(sr,0),(sl,0),(sl,C)])
def letter_Z():
    """diagonal's edges through (D, cap-A) and (0, F), A apart; top bar starts E in from the left"""
    D,E,F=DIM['D'],DIM['E'],DIM['F']; G=A
    P1=np.array([D,CAP-G]); P2=np.array([0,F]); dx,dy=P1-P2
    th=brentq(lambda t: dx*np.sin(t)-dy*np.cos(t)-A, np.arctan2(dy,dx), np.pi/2); u=np.array([np.cos(th),np.sin(th)])
    r_bot=line_inter(P1,u,(0,A),(1,0)); l_top=line_inter(P2,u,(0,CAP-A),(1,0))
    return poly([(E,CAP),(D,CAP),(D,CAP-G),tuple(r_bot),(D,A),(D,0),(0,0),(0,F),tuple(l_top),(E,CAP-A)])
def fig1():
    """flag: flat top G long, concave arc J down to the flag's left end I below the cap; flag H deep and H wide"""
    G,H,I,J=DIM['G'],DIM['H'],DIM['I'],DIM['J']
    sl=H; sr=H+A; a=np.array([0,CAP-I]); b=np.array([sr-G,CAP])
    m=(a+b)/2; ch=b-a; L=np.linalg.norm(ch); h=np.sqrt(J*J-(L/2)**2); perp=np.array([-ch[1],ch[0]])/L
    c=m+perp*h if (m+perp*h)[1]>m[1] else m-perp*h
    a0=np.arctan2(a[1]-c[1],a[0]-c[0]); a1=np.arctan2(b[1]-c[1],b[0]-c[0])
    pts=[tuple(c+J*np.array([np.cos(t),np.sin(t)])) for t in np.linspace(a1,a0,200)]
    return poly([(0,CAP-I-H),(sl,CAP-I-H),(sl,0),(sr,0),(sr,CAP)]+pts)
SPINE='lifted'
def fig2():
    """upper bowl: top arc K on the overshoot line and right arc L on its centre line, touching at the right side;
    terminal: K's left end, cut horizontally.  Spine: arc M at x=M (vertical through the bottom-left corner), its centre
    raised just enough (4.6) for L to be externally tangent to it; this fits the drawing better than a spine centred on
    the baseline (bowl 4.6 low) or on the O line (bowl 6 high)."""
    g=Grid(); D,K,L,M,O_=[DIM[k] for k in 'DKLMO']
    if SPINE=='baseline':  cM=np.array([M,0.0]); y=np.sqrt((L+M-A)**2-(M-D+L)**2)
    elif SPINE=='O':       cM=np.array([M,O_]); y=O_+np.sqrt((L+M-A)**2-(M-D+L)**2)
    elif SPINE=='lifted':  # bowl top on the overshoot line; spine centre at x=M, raised just enough for tangency
        y=TOP-K; cM=np.array([M,y-np.sqrt((L+M-A)**2-(M-D+L)**2)])
    else:                                                     # bowl top on the overshoot line, spine centre x by tangency
        y=TOP-K; cM=np.array([D-L+np.sqrt((L+M-A)**2-y*y),0.0])
    cL=np.array([D-L,y]); cK=np.array([D-K,y])
    m=ring(g,cM,M,M-A)&g.sector(cM,ang(cM,cL),180)&(g.Y>=0)
    m|=ring(g,cL,L,L-A)&g.sector(cL,ang(cL,cM),0)
    m|=ring(g,cK,K,K-A)&g.sector(cK,0,180)
    m|=rect(g,0,D,0,A)&g.disk(cM,M)
    return m,{'bowl top vs TOP':(round(y+K,1),round(TOP,1)),'spine centre':cM.round(1),'left edge at y=0':round(cM[0]-np.sqrt(M*M-cM[1]**2),1)}
def fig3():
    """upper bowl: top arc N on the overshoot line, right arc Q centred S below it (E in from the right), N internally
    tangent to Q; lower bowl: bottom arc R on the overshoot line, right arc Q centred T above it (touching the right),
    R internally tangent to Q.  K joins the lower Q to the notch on the upper Q at height P; bar A thick from x=F,
    centred on the notch height; fillets H at the bar.  Terminals: N's and R's left ends, cut horizontally."""
    D,E,F,H,K,N,P,Q,R,S,T=[DIM[k] for k in 'DEFHKNPQRST']
    cQu=np.array([D-E-Q,TOP-S]); yN=TOP-N; cN=np.array([cQu[0]+np.sqrt((Q-N)**2-(yN-cQu[1])**2),yN])
    cQl=np.array([D-Q,BOT+T]); yR=BOT+R; cR=np.array([cQl[0]+np.sqrt((Q-R)**2-(yR-cQl[1])**2),yR])
    yn=BOT+P; bar_bot,bar_top=yn-A/2,yn+A/2
    npt=np.array([cQu[0]+np.sqrt(Q*Q-(yn-cQu[1])**2),yn])
    cK=np.array(fsolve(lambda v:[np.hypot(*(v-cQl))-(Q-K),np.hypot(*(v-npt))-K],[cQl[0]+230,250.0]))
    yG=bar_top+H; cG=np.array([cQu[0]+np.sqrt((Q-A-H)**2-(yG-cQu[1])**2),yG])
    yF=bar_bot-H; cF=np.array([cQl[0]+np.sqrt((Q-A-H)**2-(yF-cQl[1])**2),yF])
    a_QN,a_QR,a_QK=ang(cQu,cN),ang(cQl,cR),ang(cQl,cK)
    pts=[arc(cN,N,180,a_QN,False),arc(cQu,Q,a_QN,ang(cQu,npt),False),arc(cK,K,ang(cK,npt),a_QK,False),
         arc(cQl,Q,a_QK,a_QR,False),arc(cR,R,a_QR,-180,False),arc(cR,R-A,-180,a_QR,True),
         arc(cQl,Q-A,a_QR,ang(cQl,cF),True),arc(cF,H,ang(cQl,cF),90,True),[(F,bar_bot),(F,bar_top)],
         arc(cG,H,-90,ang(cQu,cG),True),arc(cQu,Q-A,ang(cQu,cG),a_QN,True),arc(cN,N-A,a_QN,180,True)]
    P_=np.vstack([np.atleast_2d(p) for p in pts])
    return Region(Polygon(P_).buffer(0)),{'notch':npt.round(1),'D-H':round(D-H,1),'terminals x':(round(cN[0]-N,1),round(cR[0]-R,1))}
def fig4():
    """stem A at the right (right edge W-U); top flat F wide; diagonal V wide from the top flat's left corner to the left
    edge A above the crossbar's bottom; crossbar V thick, bottom at X, running U past the stem"""
    F,U,V,W,X=DIM['F'],DIM['U'],DIM['V'],DIM['W'],DIM['X']
    sr=W-U; sl=sr-A; yb=X+V; tl=sr-F
    p0=np.array([0,X+A]); p1=np.array([tl,CAP]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*V; ib=line_inter(qi,d,(0,yb),(1,0)); ia=line_inter(qi,d,(sl,0),(0,1))
    outer=[(0,X),(sl,X),(sl,0),(sr,0),(sr,X),(W,X),(W,yb),(sr,yb),(sr,CAP),(tl,CAP),(0,X+A)]
    return poly(outer,[[tuple(ib),(sl,yb),tuple(ia)]])
