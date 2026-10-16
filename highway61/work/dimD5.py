"""Series D, page D-5 (Y Z 1 2 3 4)."""
import numpy as np
from scipy.optimize import brentq, fsolve
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
R8={'A':'1 1/4','B':'6 3/4','C':'3 1/2','D':'5 3/8','E':'5/16','F':'1 1/16','G':'3/4','H':'1','I':'1 1/8','J':'11/16',
    'K':'1/2','L':'2 5/8','M':'3 1/8','N':'4 7/8','O':'2 1/8','P':'2 5/16','Q':'1/8','R':'1 5/8','S':'2 7/16',
    'T':'2 11/16','U':'2 3/4','V':'1 3/16','W':'2','X':'4 1/2','Y':'5 1/8'}
R12={'A':'1 7/8','B':'10 1/8','C':'5 1/4','D':'8 1/16','E':'15/32','F':'1 19/32','G':'1 1/8','H':'1 1/2','I':'1 11/16',
     'J':'1 1/32','K':'3/4','L':'3 15/16','M':'4 11/16','N':'7 5/16','O':'3 3/16','P':'3 15/32','Q':'3/16','R':'2 7/16',
     'S':'3 21/32','T':'4 1/32','U':'4 1/8','V':'1 25/32','W':'3','X':'6 3/4','Y':'7 11/16'}
R18={'A':'2 13/16','B':'15 3/16','C':'7 7/8','D':'12 3/32','E':'11/16','F':'2 3/8','G':'1 11/16','H':'2 1/4','I':'2 17/32',
     'J':'1 9/16','K':'1 1/8','L':'5 29/32','M':'7 1/32','N':'10 31/32','O':'4 25/32','P':'5 3/16','Q':'9/32','R':'3 21/32',
     'S':'5 1/2','T':'6 1/16','U':'6 3/16','V':'2 11/16','W':'4 1/2','X':'10 1/8','Y':'11 17/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['Q']; TOP=CAP+OV; BOT=-OV
def line_inter(p,d,q,e):
    d,e=np.asarray(d,float),np.asarray(e,float); M=np.array([d,-e]).T; s,_=np.linalg.solve(M,np.asarray(q,float)-np.asarray(p,float)); return np.asarray(p,float)+np.asarray(d,float)*s
def poly(pts,holes=()):
    g=Polygon([tuple(map(float,p)) for p in pts])
    for h in holes: g=g.difference(Polygon([tuple(map(float,p)) for p in h]))
    return Region(g)
def letter_Y():
    B,C=DIM['B'],DIM['C']; c=B/2; sl,sr=c-A/2,c+A/2
    p0=np.array([0,CAP]); p1=np.array([sl,C]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*A
    xt=line_inter(qi,d,(0,CAP),(1,0))[0]; v=line_inter(qi,d,(c,0),(0,1))
    return poly([(0,CAP),(xt,CAP),tuple(v),(B-xt,CAP),(B,CAP),(sr,C),(sr,0),(sl,0),(sl,C)])
def letter_Z():
    """diagonal's edges through (D, cap-G) and (0, F), A apart; top bar starts E in from the left"""
    D,E,F,G=DIM['D'],DIM['E'],DIM['F'],DIM['G']
    P1=np.array([D,CAP-G]); P2=np.array([0,F]); dx,dy=P1-P2
    th=brentq(lambda t: dx*np.sin(t)-dy*np.cos(t)-A, np.arctan2(dy,dx), np.pi/2)
    u=np.array([np.cos(th),np.sin(th)])
    r_bot=line_inter(P1,u,(0,A),(1,0)); l_top=line_inter(P2,u,(0,CAP-A),(1,0))
    return poly([(E,CAP),(D,CAP),(D,CAP-G),tuple(r_bot),(D,A),(D,0),(0,0),(0,F),tuple(l_top),(E,CAP-A)]),{'angle':round(np.degrees(th),2)}
def fig1():
    """flag: flat top H long, then a concave arc I down to the flag's left end K below the cap; flag H deep, J wide"""
    H,I,J,K=DIM['H'],DIM['I'],DIM['J'],DIM['K']
    sl=J; sr=J+A; a=np.array([0,CAP-K]); b=np.array([sr-H,CAP])
    m=(a+b)/2; ch=b-a; L=np.linalg.norm(ch); h=np.sqrt(I*I-(L/2)**2); perp=np.array([-ch[1],ch[0]])/L
    c=m+perp*h if (m+perp*h)[1]>m[1] else m-perp*h              # centre above the chord (concave flag)
    a0=np.arctan2(a[1]-c[1],a[0]-c[0]); a1=np.arctan2(b[1]-c[1],b[0]-c[0])
    ts=np.linspace(a1,a0 if a0<a1 else a0-2*np.pi,200) if False else np.linspace(a1,a0,200)
    arc=[tuple(c+I*np.array([np.cos(t),np.sin(t)])) for t in ts]
    return poly([(0,CAP-K-H),(sl,CAP-K-H),(sl,0),(sr,0),(sr,CAP)]+arc),{'centre':c.round(1)}
def fig2():
    """spine: arc N centred on the baseline at x=N (vertical at the bottom-left corner); bowl P touches the right side and
    is externally tangent to the spine, which fixes its height; top arc M touches the overshoot line, internally tangent
    to P; arc L internally tangent to M and ending at the terminal (left edge, O below the overshoot line)."""
    g=Grid(); D,L,M,N,O,P=[DIM[k] for k in 'DLMNOP']
    cN=np.array([N,0.0])
    cP=np.array([D-P,np.sqrt((P+N-A)**2-(N-D+P)**2)])
    yM=TOP-M; cM=np.array([cP[0]-np.sqrt((M-P)**2-(cP[1]-yM)**2),yM])
    yt=TOP-O; tip=np.array([0.0,yt])
    def f(v):                                                   # L's centre: L from the tip, M-L from M's centre
        return [np.hypot(*(v-tip))-L, np.hypot(*(v-cM))-(M-L)]
    cL=np.array(fsolve(f,[cM[0],cM[1]+40]))
    m=ring(g,cN,N,N-A)&g.sector(cN,ang(cN,cP),180)&(g.Y>=0)
    m|=ring(g,cP,P,P-A)&g.sector(cP,ang(cP,cN),ang(cM,cP))
    m|=ring(g,cM,M,M-A)&g.sector(cM,ang(cM,cP),ang(cM,cL))
    m|=ring(g,cL,L,L-A)&g.sector(cL,ang(cM,cL),ang(cL,tip))
    m|=rect(g,0,D,0,A)&g.disk(cN,N)
    return m,{'cP':cP.round(1),'cM':cM.round(1),'cL':cL.round(1)}
def fig3():
    """upper bowl S touching the top overshoot line E in from the right; lower bowl T touching the bottom overshoot line
    and the right side.  Inner: S-A and T-A, with V (level with T's centre, tangent to T-A) under the bar and a fillet K
    above it.  Bar A thick from x=U, its bottom on V; notch where the outer arcs cross (reported against X).
    Terminals R below / O above the overshoot lines, cut along the radius."""
    g=Grid(); D,E,K,O,R,S,T,U,V,X=[DIM[k] for k in 'DEKORSTUVX']
    cS=np.array([D-E-S,TOP-S]); cT=np.array([D-T,BOT+T])
    cV=np.array([cT[0]+(T-A-V),cT[1]])                         # tangent to T-A on the right
    bar_bot=cV[1]+V; bar_top=bar_bot+A
    yk=bar_top+K; cK=np.array([cS[0]+np.sqrt((S-A-K)**2-(yk-cS[1])**2),yk])
    # notch: crossing of the outer circles on the right
    d=cS-cT; dist=np.linalg.norm(d); a=(T*T-S*S+dist*dist)/(2*dist); h=np.sqrt(T*T-a*a)
    base=cT+d/dist*a; nv=np.array([d[1],-d[0]])/dist; notch=base+nv*h if (base+nv*h)[0]>base[0] else base-nv*h
    yt=TOP-R; xt=cS[0]-np.sqrt(S*S-(yt-cS[1])**2)
    yb=BOT+O; xb=cT[0]-np.sqrt(T*T-(yb-cT[1])**2)
    a_t=ang(cS,(xt,yt)); a_b=ang(cT,(xb,yb))
    up=g.disk(cS,S)&(g.Y>=notch[1]); lo=g.disk(cT,T)&(g.Y<=notch[1])
    ink=(up&g.sector(cS,-90,a_t))|(lo&g.sector(cT,a_b,90))|(rect(g,U,D,bar_bot,bar_top)&(up|lo))
    up_cnt=(g.disk(cS,S-A)&(g.Y>=bar_top)&~((g.X>=cK[0])&(g.Y<=cK[1])))|g.disk(cK,K)
    lo_cnt=(g.disk(cT,T-A)&(g.Y<=bar_bot)&~((g.X>=cV[0])&(g.Y>=cV[1])))|g.disk(cV,V)
    m=ink&~up_cnt&~lo_cnt
    return m,{'notch':notch.round(1),'BOT+X':round(BOT+X,1),'bar':(round(bar_bot,1),round(bar_top,1)),'tips x':(round(xt,1),round(xb,1)),'cK':cK.round(1)}
def fig4():
    """stem A at the right (right edge Y); diagonal I wide from the stem's top-left corner to the left edge G above the
    crossbar; crossbar I thick, bottom at W, running G past the stem"""
    G,I,W,Y=DIM['G'],DIM['I'],DIM['W'],DIM['Y']
    sr=Y; sl=Y-A; yb=W+I
    p0=np.array([0,W+G]); p1=np.array([sl,CAP]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*I
    ib=line_inter(qi,d,(0,yb),(1,0)); ia=line_inter(qi,d,(sl,0),(0,1))
    outer=[(0,W),(sl,W),(sl,0),(sr,0),(sr,W),(Y+G,W),(Y+G,yb),(sr,yb),(sr,CAP),(sl,CAP),(0,W+G)]
    return poly(outer,[[tuple(ib),(sl,yb),tuple(ia)]])
