"""Series C, page C-5 (Y Z 1 2 3 4)."""
import numpy as np
from scipy.optimize import brentq
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
from dimD5 import line_inter, poly
R8={'A':'1 1/8','B':'5','C':'3 5/8','D':'4 3/8','E':'5/16','F':'1 1/4','G':'3/4','H':'1/2','I':'7/8','J':'1/8','K':'3/16',
    'L':'3','M':'6','N':'1 15/16','O':'1/4','P':'1 13/16','Q':'6 5/8','R':'2 3/16','S':'1 3/4','T':'1 7/8','U':'15/16',
    'V':'13/16','W':'1 1/16'}
R12={'A':'1 11/16','B':'7 1/2','C':'5 7/16','D':'6 9/16','E':'15/32','F':'1 7/8','G':'1 1/8','H':'3/4','I':'1 5/16',
     'J':'3/16','K':'9/32','L':'4 1/2','M':'9','N':'2 15/16','O':'3/8','P':'2 23/32','Q':'9 15/16','R':'3 9/32','S':'2 5/8',
     'T':'2 13/16','U':'1 13/32','V':'1 7/32','W':'1 19/32'}
R18={'A':'2 17/32','B':'11 1/4','C':'8 5/32','D':'9 27/32','E':'11/16','F':'2 13/16','G':'1 11/16','H':'1 1/8','I':'1 31/32',
     'J':'9/32','K':'7/16','L':'6 3/4','M':'13 1/2','N':'4 7/16','O':'9/16','P':'4 1/16','Q':'14 29/32','R':'4 15/16',
     'S':'3 15/16','T':'4 7/32','U':'2 1/8','V':'1 13/16','W':'2 3/8'}
# the 18-inch row prints D as 9 7/8; the 8- and 12-inch rows (4 3/8, 6 9/16) and the C-1..C-4 widths give 9 27/32
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
def letter_Y():
    B,C=DIM['B'],DIM['C']; c=B/2; sl,sr=c-A/2,c+A/2
    p0=np.array([0,CAP]); p1=np.array([sl,C]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*A; xt=line_inter(qi,d,(0,CAP),(1,0))[0]; v=line_inter(qi,d,(c,0),(0,1))
    return poly([(0,CAP),(xt,CAP),tuple(v),(B-xt,CAP),(B,CAP),(sr,C),(sr,0),(sl,0),(sl,C)])
def letter_Z():
    """diagonal's edges through (D, cap-G) and (0, F), A apart; top bar starts E in from the left"""
    D,E,F,G=DIM['D'],DIM['E'],DIM['F'],DIM['G']
    P1=np.array([D,CAP-G]); P2=np.array([0,F]); dx,dy=P1-P2
    th=brentq(lambda t: dx*np.sin(t)-dy*np.cos(t)-A, np.arctan2(dy,dx), np.pi/2); u=np.array([np.cos(th),np.sin(th)])
    r_bot=line_inter(P1,u,(0,A),(1,0)); l_top=line_inter(P2,u,(0,CAP-A),(1,0))
    return poly([(E,CAP),(D,CAP),(D,CAP-G),tuple(r_bot),(D,A),(D,0),(0,0),(0,F),tuple(l_top),(E,CAP-A)])
def fig1():
    """flag: flat top G long, concave arc F down to the flag's left end H below the cap; flag I deep, H wide"""
    G,F,H,I=DIM['G'],DIM['F'],DIM['H'],DIM['I']
    sl=H; sr=H+A; a=np.array([0,CAP-H]); b=np.array([sr-G,CAP])
    m=(a+b)/2; ch=b-a; L=np.linalg.norm(ch); h=np.sqrt(F*F-(L/2)**2); perp=np.array([-ch[1],ch[0]])/L
    c=m+perp*h if (m+perp*h)[1]>m[1] else m-perp*h
    a0=np.arctan2(a[1]-c[1],a[0]-c[0]); a1=np.arctan2(b[1]-c[1],b[0]-c[0])
    arc=[tuple(c+F*np.array([np.cos(t),np.sin(t)])) for t in np.linspace(a1,a0,200)]
    return poly([(0,CAP-H-I),(sl,CAP-H-I),(sl,0),(sr,0),(sr,CAP)]+arc)
def fig2():
    """spine: arc M centred on the baseline at x=M; bowl R touching the right side, externally tangent to the spine
    (which sets its height); terminal: the bowl's left end, cut horizontally at its centre height"""
    g=Grid(); D,M,R=DIM['D'],DIM['M'],DIM['R']
    cM=np.array([M,0.0]); cR=np.array([D-R,np.sqrt((R+M-A)**2-(M-D+R)**2)])
    m=ring(g,cM,M,M-A)&g.sector(cM,ang(cM,cR),180)&(g.Y>=0)
    m|=ring(g,cR,R,R-A)&g.sector(cR,ang(cR,cM),180)
    m|=rect(g,0,D,0,A)&g.disk(cM,M)
    return m,{'bowl top vs TOP':(round(cR[1]+R,1),round(TOP,1)),'terminal x':round(cR[0]-R,1)}
NOTCH='on L'
def fig3():
    """upper bowl: top arc N (top on the overshoot line, K in from the right) and side arc L on the same centre line,
    touching at the right; lower bowl: bottom arc R and side arc Q likewise, touching the right side.  T joins Q to the
    notch on L at the notch height D (internally tangent to Q).  Counters: N-A, L-A and Q-A, R-A, each with a fillet G
    at the bar.  Bar A thick from x=S, centred on the notch height D above the
    bottom overshoot line.  Upper terminal at (O, top-P), cut along the radius; lower one the bowl's left end, cut
    horizontally.  Traced as one contour."""
    from dimE4 import arc
    D,G,K,L,N,O,P,Q,R,S,T,U=[DIM[k] for k in 'DGKLNOPQRSTU']
    cN=np.array([D-K-N,TOP-N]); cL=np.array([D-K-L,cN[1]])
    cR=np.array([D-R,BOT+R]); cQ=np.array([D-Q,cR[1]])
    yn=BOT+D; bar_bot,bar_top=yn-A/2,yn+A/2
    if NOTCH=='on L':                                             # T through L at the notch height D
        npt=np.array([cL[0]+np.sqrt(L*L-(yn-cL[1])**2),yn])
    else:                                                         # T through the table's notch (D-U, D)
        npt=np.array([D-U,yn])
    def f(v): return [np.hypot(*(v-cQ))-(Q-T),np.hypot(*(v-npt))-T]
    from scipy.optimize import fsolve
    cT=np.array(fsolve(f,[cQ[0]+400,250.0]))
    yG=bar_top+G; cG=np.array([cL[0]+np.sqrt((L-A-G)**2-(yG-cL[1])**2),yG])
    d=cL-cT; dist=np.linalg.norm(d); a=(T*T-L*L+dist*dist)/(2*dist); h=np.sqrt(T*T-a*a)
    base=cT+d/dist*a; nv=np.array([d[1],-d[0]])/dist; notch=base+nv*h if (base+nv*h)[0]>base[0] else base-nv*h
    yF=bar_bot-G; cF=np.array([cQ[0]+np.sqrt((Q-A-G)**2-(yF-cQ[1])**2),yF]); a_QF=ang(cQ,cF)  # lower fillet: Q-A to the bar
    tip=np.array([O,TOP-P]); a_tip=ang(cN,tip); a_TQ=ang(cQ,cT); a_LG=ang(cL,cG)
    pts=[arc(cN,N,a_tip,0,False),arc(cL,L,0,ang(cL,notch),False),arc(cT,T,ang(cT,notch),a_TQ,False),
         arc(cQ,Q,a_TQ,0,False),arc(cR,R,0,-180,False),arc(cR,R-A,-180,0,True),arc(cQ,Q-A,0,a_QF,True),
         arc(cF,G,a_QF,90,True),[(S,bar_bot),(S,bar_top)],arc(cG,G,-90,a_LG,True),arc(cL,L-A,a_LG,0,True),
         arc(cN,N-A,0,a_tip,True)]
    P_=np.vstack([np.atleast_2d(p) for p in pts])
    return Region(Polygon(P_).buffer(0)),{'notch':notch.round(1),'table notch':(round(D-U,1),round(yn,1)),'cT':cT.round(1),
              'tip on N (dist-N)':round(np.linalg.norm(tip-cN)-N,1),'lower terminal x':round(cR[0]-R,1)}
def fig4():
    """as Series D: stem A at the right (right edge D); diagonal W wide from the stem's top-left corner to the left edge
    V above the crossbar; crossbar W thick, bottom at N, running H past the stem"""
    D,H,N,V,W=DIM['D'],DIM['H'],DIM['N'],DIM['V'],DIM['W']
    sr=D; sl=D-A; yb=N+W
    p0=np.array([0,N+V]); p1=np.array([sl,CAP]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*W; ib=line_inter(qi,d,(0,yb),(1,0)); ia=line_inter(qi,d,(sl,0),(0,1))
    outer=[(0,N),(sl,N),(sl,0),(sr,0),(sr,N),(D+H,N),(D+H,yb),(sr,yb),(sr,CAP),(sl,CAP),(0,N+V)]
    return poly(outer,[[tuple(ib),(sl,yb),tuple(ia)]])
