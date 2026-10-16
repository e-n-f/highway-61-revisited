"""Series F, page F-5 (Y Z 1 2 3 4).  Dimensions from the 8-inch row (agrees with the 12-inch row)."""
import numpy as np
from scipy.optimize import brentq, fsolve
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimD5 import line_inter, poly
from dimE4 import arc, internal_tangent
R8={'A':'1 1/2','B':'9','C':'3 1/4','D':'7 3/8','E':'1/4','F':'1 1/8','G':'1','H':'1 3/4','I':'1/2','J':'1/8',
    'K':'2 3/16','L':'1 5/16','M':'8','N':'4 3/4','O':'5 1/4','P':'3 7/8','Q':'3 11/16','R':'4 1/2','S':'3/4',
    'T':'2','U':'1 7/16'}
R12={'A':'2 1/4','B':'13 1/2','C':'4 7/8','D':'11 1/16','E':'3/8','F':'1 11/16','G':'1 1/2','H':'2 5/8','I':'3/4',
     'J':'3/16','K':'3 9/32','L':'1 31/32','M':'12','N':'7 1/8','O':'7 7/8','P':'5 13/16','Q':'5 17/32','R':'6 3/4',
     'S':'1 1/8','T':'3','U':'2 5/32'}
DIM=units(R8,8); CHECK=check({8:R8,12:R12})
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
deg=lambda v: np.degrees(np.arctan2(v[1],v[0]))
def letter_Y():
    B,C=DIM['B'],DIM['C']; c=B/2; sl,sr=c-A/2,c+A/2
    p0=np.array([0,CAP]); p1=np.array([sl,C]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*A; xt=line_inter(qi,d,(0,CAP),(1,0))[0]; v=line_inter(qi,d,(c,0),(0,1))
    return poly([(0,CAP),(xt,CAP),tuple(v),(B-xt,CAP),(B,CAP),(sr,C),(sr,0),(sl,0),(sl,C)])
def letter_Z():
    """diagonal's edges through (D, cap-A) and (0, A), A apart; top bar starts E in from the left"""
    D,E=DIM['D'],DIM['E']; F=G=A
    P1=np.array([D,CAP-G]); P2=np.array([0,F]); dx,dy=P1-P2
    th=brentq(lambda t: dx*np.sin(t)-dy*np.cos(t)-A, np.arctan2(dy,dx), np.pi/2); u=np.array([np.cos(th),np.sin(th)])
    r_bot=line_inter(P1,u,(0,A),(1,0)); l_top=line_inter(P2,u,(0,CAP-A),(1,0))
    return poly([(E,CAP),(D,CAP),(D,CAP-G),tuple(r_bot),(D,A),(D,0),(0,0),(0,F),tuple(l_top),(E,CAP-A)])
def fig1():
    """flag: flat top F long, concave arc H down to the flag's left end G below the cap; flag F deep and F wide"""
    F,G,H=DIM['F'],DIM['G'],DIM['H']
    sl=F; sr=F+A; a=np.array([0,CAP-G]); b=np.array([sr-F,CAP])
    m=(a+b)/2; ch=b-a; L=np.linalg.norm(ch); h=np.sqrt(H*H-(L/2)**2); perp=np.array([-ch[1],ch[0]])/L
    c=m+perp*h if (m+perp*h)[1]>m[1] else m-perp*h
    a0=np.arctan2(a[1]-c[1],a[0]-c[0]); a1=np.arctan2(b[1]-c[1],b[0]-c[0])
    pts=[tuple(c+H*np.array([np.cos(t),np.sin(t)])) for t in np.linspace(a1,a0,200)]
    return poly([(0,CAP-G-F),(sl,CAP-G-F),(sl,0),(sr,0),(sr,CAP)]+pts)
def fig2():
    """top arc O centred on the line x=Q, its top on the overshoot line; upper terminal L below it, cut along the
    radius; bowl K touching the right side, internally tangent to O; spine arc N centred on the baseline at x=N
    (vertical at the bottom-left corner); the spine's edges are the internal tangents between the bowl and the spine
    (they stand 10.7 apart)"""
    D,K,L,N,O,Q=[DIM[k] for k in 'DKLNOQ']
    cO=np.array([Q,TOP-O]); cK=np.array([D-K,cO[1]+np.sqrt((O-K)**2-(D-K-Q)**2)]); cN=np.array([N,0.0])
    yt=TOP-L; xt=cO[0]-np.sqrt(O*O-(yt-cO[1])**2); a_t=ang(cO,(xt,yt)); a_OK=ang(cO,cK)
    T1,T2,n1=internal_tangent(cK,K,cN,N-A,lambda s: max(s,key=lambda v:v[0]))
    T3,T4,n2=internal_tangent(cN,N,cK,K-A,lambda s: min(s,key=lambda v:v[0]))
    xf=cN[0]-np.sqrt((N-A)**2-A*A)
    P=[arc(cO,O,a_t,a_OK,False),arc(cK,K,a_OK,deg(n1),False),[T2],arc(cN,N-A,deg(-n1),ang(cN,(xf,A)),True),
       [(D,A),(D,0),(0,0)],arc(cN,N,180,deg(n2),False),[T4],arc(cK,K-A,deg(-n2),a_OK,True),arc(cO,O-A,a_OK,a_t,True)]
    poly_=np.vstack([np.atleast_2d(p) for p in P])
    return Region(Polygon(poly_).buffer(0)),{'cK':cK.round(1),'tip':(round(xt,1),round(yt,1)),
        'bowl/spine gap':round(np.linalg.norm(cK-cN)-(K+N-A),1)}
def fig3():
    """outer: arcs O on the axis x=D-Q, top and bottom on the overshoot lines; upper bowl T (I in from the right) and
    lower bowl K (touching the right side), each internally tangent to its O; they cross at the notch.  Bar A thick from
    x=P, centred on the notch height R.  Counters: O-A, then connecting arcs (table S, F; refitted) tangent to O-A and the
    bar, A in from the bowls' right extremes.  Terminals:
    upper F below the top overshoot line, lower A above the bottom one, cut along the radius."""
    D,F,I,K,O,P,Q,R,S,T=[DIM[k] for k in 'DFIKOPQRST']
    xa=D-Q; cOu=np.array([xa,TOP-O]); cOl=np.array([xa,BOT+O])
    cT=np.array([D-I-T,cOu[1]+np.sqrt((O-T)**2-(D-I-T-xa)**2)])
    cK=np.array([D-K,cOl[1]-np.sqrt((O-K)**2-(D-K-xa)**2)])
    yn=BOT+R; bb,bt=yn-A/2,yn+A/2
    # counter ends: connecting arcs tangent to the bar, internally tangent to O-A, and A in from the bowl's right extreme.
    # The table's S and F cannot be tangent to both (they only fit a straight-sided slot); the drawing shows the curve,
    # so the radii are refitted (as connecting arcs): solve (xr-r-xa)^2 + (yb +/- r - yO)^2 = (O-A-r)^2
    def fillet(xr,ybar,cO,up):
        s_=1 if up else -1
        f=lambda r: (xr-r-xa)**2+(ybar+s_*r-cO[1])**2-(O-A-r)**2
        from scipy.optimize import brentq as bq
        r=bq(f,1.0,O-A-1); return np.array([xr-r,ybar+s_*r]),r
    cS,S=fillet(D-I-A,bt,cOu,True); cF,F_=fillet(D-A,bb,cOl,False)
    d=cT-cK; dist=np.linalg.norm(d); a=(K*K-T*T+dist*dist)/(2*dist); h=np.sqrt(K*K-a*a)
    base=cK+d/dist*a; nv=np.array([d[1],-d[0]])/dist; notch=base+nv*h if (base+nv*h)[0]>base[0] else base-nv*h
    yu=TOP-F; tu=np.array([xa-np.sqrt(O*O-(yu-cOu[1])**2),yu]); au=ang(cOu,tu)
    yl=BOT+A; tl=np.array([xa-np.sqrt(O*O-(yl-cOl[1])**2),yl]); al=ang(cOl,tl)
    P_=[arc(cOu,O,au,ang(cOu,cT),False),arc(cT,T,ang(cOu,cT),ang(cT,notch),False),arc(cK,K,ang(cK,notch),ang(cOl,cK),False),
        arc(cOl,O,ang(cOl,cK),al,False),arc(cOl,O-A,al,ang(cOl,cF),True),arc(cF,F_,ang(cOl,cF),90,True),[(P,bb),(P,bt)],
        arc(cS,S,-90,ang(cOu,cS),True),arc(cOu,O-A,ang(cOu,cS),au,True)]
    pts=np.vstack([np.atleast_2d(p) for p in P_])
    return Region(Polygon(pts).buffer(0)),{'notch':notch.round(1),'table (D-F, BOT+R)':(round(D-F,1),round(yn,1)),
        'tips':(tu.round(1),tl.round(1)),'counter ends S, F refit':(round(S,1),round(F_,1)),'table S, F':(round(DIM['S'],1),round(F,1))}
def fig4():
    """stem A at the right (right edge M-G), its top the flat; diagonal U wide from the stem's top-left corner to the
    left edge at the crossbar's top; crossbar U thick, bottom A above the baseline, running G past the stem"""
    G,M,U=DIM['G'],DIM['M'],DIM['U']
    sr=M-G; sl=sr-A; yb=A+U
    p0=np.array([0,A+U]); p1=np.array([sl,CAP]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n
    qi=p0+n*U; ib=line_inter(qi,d,(0,yb),(1,0)); ia=line_inter(qi,d,(sl,0),(0,1))
    outer=[(0,A),(sl,A),(sl,0),(sr,0),(sr,A),(M,A),(M,yb),(sr,yb),(sr,CAP),(sl,CAP),(0,A+U)]
    return poly(outer,[[tuple(ib),(sl,yb),tuple(ia)]])
