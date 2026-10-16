"""Series D, page D-4 (S T U V W X)."""
import numpy as np
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip, cross_solve
from dimE4 import arc, internal_tangent
R8={'A':'1 1/4','B':'5 3/8','C':'2 7/16','D':'2','E':'3 1/8','F':'2 1/2','G':'5/16','H':'1 5/8','I':'1/8','J':'4 7/8',
    'K':'1 1/8','L':'2 11/16','M':'6','N':'1','O':'7','P':'5 1/8'}
R12={'A':'1 7/8','B':'8 1/16','C':'3 21/32','D':'3','E':'4 11/16','F':'3 3/4','G':'15/32','H':'2 7/16','I':'3/16',
     'J':'7 5/16','K':'1 11/16','L':'4 1/32','M':'9','N':'1 1/2','O':'10 1/2','P':'7 11/16'}
R18={'A':'2 13/16','B':'12 3/32','C':'5 1/2','D':'4 1/2','E':'7 1/32','F':'5 5/8','G':'11/16','H':'3 21/32','I':'9/32',
     'J':'10 31/32','K':'2 17/32','L':'6 3/64','M':'13 1/2','N':'2 1/4','O':'15 3/4','P':'11 17/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['I']; TOP=CAP+OV; BOT=-OV
deg=lambda v: np.degrees(np.arctan2(v[1],v[0]))
def letter_S():
    """top arc C and upper-left bowl D centred on one horizontal (they touch on the left, G in from the left edge);
    lower-right bowl C and bottom arc E centred on one vertical (they touch at the bottom), the bowl touching the right
    edge.  Terminals: upper H below the top overshoot line, lower D above the bottom one, cut along the radius.
    The spine is the pair of internal tangents between the two bowls."""
    B,C,D,E,G,H=[DIM[k] for k in 'BCDEGH']
    cT=np.array([G+C,TOP-C]); cD=np.array([G+D,TOP-C])
    cL=np.array([B-C,BOT+C]); cB=np.array([B-C,BOT+E])
    yt=TOP-H; xt=cT[0]+np.sqrt(C*C-(yt-cT[1])**2); a_t=deg(np.array([xt,yt])-cT)
    yb=BOT+D; xb=cB[0]-np.sqrt(E*E-(yb-cB[1])**2); a_b=deg(np.array([xb,yb])-cB)
    T1,T2,n1=internal_tangent(cD,D,cL,C-A,lambda s: min(s,key=lambda v:v[1]))
    T3,T4,n2=internal_tangent(cL,C,cD,D-A,lambda s: max(s,key=lambda v:v[1]))
    P=[arc(cT,C,a_t,180,True), arc(cD,D,180,deg(n1),True), [T2], arc(cL,C-A,deg(-n1),-90,False),
       arc(cB,E-A,-90,a_b,False), arc(cB,E,a_b,-90,True), arc(cL,C,-90,deg(n2),True), [T4],
       arc(cD,D-A,deg(-n2),180,False), arc(cT,C-A,180,a_t,False)]
    poly=np.vstack([np.atleast_2d(p) for p in P])
    return Region(Polygon(poly).buffer(0)),{'upper tip x':round(xt,1),'B-G':round(B-G,1),'lower tip x':round(xb,1)}
def letter_T():
    g=Grid(); J,K=DIM['J'],DIM['K']
    return rect(g,0,J,CAP-K,CAP)|rect(g,J/2-A/2,J/2+A/2,0,CAP)
def letter_U():
    g=Grid(); B,L=DIM['B'],DIM['L']; c=np.array([B/2,BOT+L])
    m=rect(g,0,A,c[1],CAP)|rect(g,B-A,B,c[1],CAP)|(ring(g,c,L,L-A)&(g.Y<=c[1]))
    return m,{'L vs B/2':(round(L,1),round(B/2,1))}
def letter_V():
    g=Grid(); M,N=DIM['M'],DIM['N']
    l0,l1=np.array([0,CAP]),np.array([(M-N)/2,0.0]); r0,r1=np.array([M,CAP]),np.array([(M+N)/2,0.0])
    ul=(l1-l0)/np.linalg.norm(l1-l0); ur=(r1-r0)/np.linalg.norm(r1-r0)
    m=(strip(g,l0,ul,A,(M/2,CAP))&(g.X<=M/2))|(strip(g,r0,ur,A,(M/2,CAP))&(g.X>=M/2))
    return m&(g.Y>=0)&(g.Y<=CAP)
def letter_W():
    """feet K wide, F apart in the middle; P: height of the middle counter's apex"""
    g=Grid(); O,K,F,P=DIM['O'],DIM['K'],DIM['F'],DIM['P']
    x1=(O-2*K-F)/2
    s1=strip(g,np.array([0,CAP]),(np.array([x1,0])-np.array([0,CAP]))/np.hypot(x1,CAP),A,(O/2,CAP))
    a2=np.array([x1+K,0.0]); ap=np.array([O/2,P])
    u2=(ap-a2)/np.linalg.norm(ap-a2); a3=np.array([O-x1-K,0.0]); u3=np.array([-u2[0],u2[1]])
    s2=strip(g,a2,u2,A,(0,CAP))&(g.X<=O/2)
    s3=strip(g,a3,u3,A,(O,CAP))&(g.X>=O/2)
    s4=strip(g,np.array([O,CAP]),(np.array([O-x1,0])-np.array([O,CAP]))/np.hypot(x1,CAP),A,(O/2,CAP))
    m=(s1&(g.X<=x1+K))|s2|s3|(s4&(g.X>=O-x1-K))
    # the steep inner strokes would otherwise poke past the feet's outer corners: keep within the outer strokes' edges
    inside=strip(g,np.array([0,CAP]),np.array([x1,-CAP])/np.hypot(x1,CAP),1e4,(O/2,CAP))
    inside&=strip(g,np.array([O,CAP]),np.array([-x1,-CAP])/np.hypot(x1,CAP),1e4,(O/2,CAP))
    return m&inside&(g.Y>=0)&(g.Y<=CAP)
def letter_X():
    """each stroke from a top corner I in from the side to the opposite bottom corner, A wide"""
    g=Grid(); B,I=DIM['B'],DIM['I']
    p=np.array([I,CAP]); q=np.array([B,0.0])
    th=cross_solve(p,q,A,-np.pi/2+0.01,-0.01); u=np.array([np.cos(th),np.sin(th)])
    s1=strip(g,p,u,A,q)
    p2=np.array([B-I,CAP]); q2=np.array([0.0,0]); u2=np.array([-u[0],u[1]])
    s2=strip(g,p2,u2,A,q2)
    return (s1|s2)&(g.Y>=0)&(g.Y<=CAP),{'angle':round(np.degrees(th),1)}
