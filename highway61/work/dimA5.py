"""Series A, page A-5 (Y Z 1 2 3 4).  Y Z 1 3 4 follow Series B (dimB5) with this page's letters; the 2 has a
two-arc spine of its own."""
import numpy as np
from scipy.optimize import fsolve
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimD5 import line_inter, poly
from dimE4 import arc
R8={'A':'3/4','B':'3 1/2','C':'3 1/8','D':'3','E':'3/16','F':'5/8','G':'9/16','H':'1/2','I':'1 3/4','J':'1/8',
    'K':'11/16','L':'1 1/2','M':'5 3/4','N':'8 3/4','O':'1/4','P':'1 3/8','Q':'1 7/8','R':'1 7/16','S':'1 1/4',
    'T':'5 1/8','U':'7 1/4','V':'2 1/2','W':'4 3/8'}
R12={'A':'1 1/8','B':'5 1/4','C':'4 11/16','D':'4 1/2','E':'9/32','F':'15/16','G':'27/32','H':'3/4','I':'2 5/8',
     'J':'3/16','K':'1 1/32','L':'2 1/4','M':'8 5/8','N':'13 1/8','O':'3/8','P':'2 1/16','Q':'2 13/16','R':'2 5/32',
     'S':'1 7/8','T':'7 11/16','U':'10 7/8','V':'3 3/4','W':'6 9/16'}
R18={'A':'1 11/16','B':'7 7/8','C':'7 1/32','D':'6 3/4','E':'7/16','F':'1 13/32','G':'1 1/4','H':'1 1/8',
     'I':'3 15/16','J':'9/32','K':'1 9/16','L':'3 3/8','M':'12 15/16','N':'19 11/16','O':'9/16','P':'3 3/32',
     'Q':'4 7/32','R':'3 1/4','S':'2 13/16','T':'11 17/32','U':'16 5/16','V':'5 5/8','W':'9 27/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
def letter_Y():
    """width B; stem below C"""
    B,C=DIM['B'],DIM['C']; c=B/2; sl,sr=c-A/2,c+A/2
    p0=np.array([0,CAP]); p1=np.array([sl,C]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*A; xt=line_inter(qi,d,(0,CAP),(1,0))[0]; v=line_inter(qi,d,(c,0),(0,1))
    return poly([(0,CAP),(xt,CAP),tuple(v),(B-xt,CAP),(B,CAP),(sr,C),(sr,0),(sl,0),(sl,C)])
def letter_Z():
    """width D; diagonal's edges through (D, cap-F) and (0, F), A apart; top bar starts E in from the left"""
    D,E,F=DIM['D'],DIM['E'],DIM['F']
    from scipy.optimize import brentq
    P1=np.array([D,CAP-F]); P2=np.array([0,F]); dx,dy=P1-P2
    th=brentq(lambda t: dx*np.sin(t)-dy*np.cos(t)-A, np.arctan2(dy,dx), np.pi/2); u=np.array([np.cos(th),np.sin(th)])
    r_bot=line_inter(P1,u,(0,A),(1,0)); l_top=line_inter(P2,u,(0,CAP-A),(1,0))
    return poly([(E,CAP),(D,CAP),(D,CAP-F),tuple(r_bot),(D,A),(D,0),(0,0),(0,F),tuple(l_top),(E,CAP-A)])
def fig1():
    """flag: flat top G long, concave arc I down to the flag's left end H below the cap; flag F deep and F wide"""
    G,H,I,F=DIM['G'],DIM['H'],DIM['I'],DIM['F']
    sl=F; sr=F+A; a=np.array([0,CAP-H]); b=np.array([sr-G,CAP])
    m=(a+b)/2; ch=b-a; L=np.linalg.norm(ch); h=np.sqrt(I*I-(L/2)**2); perp=np.array([-ch[1],ch[0]])/L
    c=m+perp*h if (m+perp*h)[1]>m[1] else m-perp*h
    a0=np.arctan2(a[1]-c[1],a[0]-c[0]); a1=np.arctan2(b[1]-c[1],b[0]-c[0])
    pts=[tuple(c+I*np.array([np.cos(t),np.sin(t)])) for t in np.linspace(a1,a0,200)]
    return poly([(0,CAP-H-F),(sl,CAP-H-F),(sl,0),(sr,0),(sr,CAP)]+pts)
def fig2():
    """width D.  Bowl: arc L tangent to the top overshoot line and to the right side; terminal at its left end (the
    left side), carried straight down to I below the overshoot line.  Spine: arc M centred A below the top overshoot
    line, internally tangent to L (its right edge), running into arc N-A (N centred O below the baseline, on the
    right; its left edge N reaches the baseline at the left side).  M and N-A are tangent: the table leaves them 1
    unit apart, so N's centre is moved along its centre line to make them touch.  Base bar A thick."""
    g=Grid(); D,I,L,M,N,O=[DIM[k] for k in 'DILMNO']
    cL=np.array([D-L,TOP-L])
    yM=TOP-A; cM=np.array([cL[0]-np.sqrt((M-L)**2-(yM-cL[1])**2),yM])
    yN=-O; dy=cM[1]-yN; cN=np.array([cM[0]+np.sqrt((M+N-A)**2-dy*dy),yN])
    x0=cN[0]-np.sqrt(N*N-yN*yN)                                   # N's left edge at the baseline
    m=ring(g,cM,M,M-A)&g.sector(cM,ang(cM,cN),ang(cM,cL))
    m|=ring(g,cN,N,N-A)&g.sector(cN,ang(cN,cM),200)&(g.Y>=0)
    m|=ring(g,cL,L,L-A)&g.sector(cL,ang(cM,cL),180)
    m|=rect(g,cL[0]-L,cL[0]-L+A,TOP-I,cL[1])
    m|=rect(g,0,D,0,A)&g.disk(cN,N)
    return m,{'N centre x (table: from tangency)':round(cN[0],1),'left edge at y=0':round(x0,1),
              'N centre for left edge at 0':round(np.sqrt(N*N-yN*yN),1)}
NOTCH_FROM='BOT'
P_REFIT=True
def fig3():
    """width D.  Upper bowl: top arc P on the overshoot line, internally tangent to the right arc T, which is centred Q
    below the overshoot line and J in from the right; lower bowl: bottom arc R on the overshoot line, internally
    tangent to the right arc U, centred V above the overshoot line and touching the right side.  The notch is at
    (D-F, W above the bottom overshoot line); arcs S join it to T and to U (internally tangent to them).  Inside, the counters end in arcs S-A (= H)
    concentric with those, so the stroke stays A round the notch; the bar runs between them from x=S.  Terminals: P's
    left end cut at its centre height (P refitted 118.6 -> 122.8 so that it also reaches the left side, as drawn); R's left end carried straight up to I above the overshoot line."""
    D,F,H,I,P,Q,R,S,T,U,V,W=[DIM[k] for k in 'DFHIPQRSTUVW']
    J=OV
    cT=np.array([D-J-T,TOP-Q])
    if P_REFIT:                       # P tangent to the overshoot line, the left side and T: |(P, TOP-P) - cT| = T-P
        from scipy.optimize import brentq
        P=brentq(lambda p: np.hypot(p-cT[0],TOP-p-cT[1])-(T-p),20,T-1)
    yP=TOP-P; cP=np.array([cT[0]+np.sqrt((T-P)**2-(yP-cT[1])**2),yP])
    cU=np.array([D-U,BOT+V]); yR=BOT+R; cR=np.array([cU[0]+np.sqrt((U-R)**2-(yR-cU[1])**2),yR])
    npt=np.array([D-F,W if NOTCH_FROM=='baseline' else BOT+W])
    cSu=np.array(fsolve(lambda v:[np.hypot(*(v-cT))-(T-S),np.hypot(*(v-npt))-S],[npt[0]-S*0.7,npt[1]+S*0.7]))
    cSl=np.array(fsolve(lambda v:[np.hypot(*(v-cU))-(U-S),np.hypot(*(v-npt))-S],[npt[0]-S*0.7,npt[1]-S*0.7]))
    r=S-A; bar_top=cSu[1]-r; bar_bot=cSl[1]+r
    a_TP,a_UR=ang(cT,cP),ang(cU,cR); a_TS,a_US=ang(cT,cSu),ang(cU,cSl)
    xl=cR[0]-R; yI=BOT+I
    pts=[arc(cP,P,180,a_TP,False),arc(cT,T,a_TP,a_TS,False),arc(cSu,S,a_TS,ang(cSu,npt),False),
         arc(cSl,S,ang(cSl,npt),a_US,False),arc(cU,U,a_US,a_UR,False),arc(cR,R,a_UR,-180,False),
         [(xl,yI),(xl+A,yI)],arc(cR,R-A,180,a_UR,True),arc(cU,U-A,a_UR,a_US,True),arc(cSl,r,a_US,90,True),
         [(S,bar_bot),(S,bar_top)],arc(cSu,r,-90,a_TS,True),arc(cT,T-A,a_TS,a_TP,True),arc(cP,P-A,a_TP,180,True)]
    P_=np.vstack([np.atleast_2d(p) for p in pts])
    return Region(Polygon(P_).buffer(0)),{'bar thickness':round(bar_top-bar_bot,1),'S-A vs H':(round(r,1),round(H,1)),
        'P':(round(P,1),round(DIM['P'],1)),'terminals x':(round(cP[0]-P,1),round(xl,1))}
def fig4():
    """D to the stem's right edge; stem A; top flat A wide; diagonal K wide from the top flat's left corner
    to the left edge K above the crossbar's bottom; crossbar K thick, bottom at Q, running H past the stem"""
    H,K,Q,D=DIM['H'],DIM['K'],DIM['Q'],DIM['D']
    sr=D; sl=sr-A; yb=Q+K; tl=sl
    p0=np.array([0,Q+K]); p1=np.array([tl,CAP]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*K; ib=line_inter(qi,d,(0,yb),(1,0)); ia=line_inter(qi,d,(sl,0),(0,1))
    outer=[(0,Q),(sl,Q),(sl,0),(sr,0),(sr,Q),(D+H,Q),(D+H,yb),(sr,yb),(sr,CAP),(tl,CAP),(0,Q+K)]
    return poly(outer,[[tuple(ib),(sl,yb),tuple(ia)]])
