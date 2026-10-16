"""Series F, page F-4 (S T U V W X)."""
import numpy as np
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip, cross_solve
from dimE4 import arc, internal_tangent
R8={'A':'1 1/2','B':'7 3/8','C':'6 1/4','D':'5 3/4','E':'15/16','F':'1 3/16','G':'1/2','H':'2 1/4','I':'2 1/8',
    'J':'1/8','K':'6 3/4','L':'1 3/8','M':'3 11/16','N':'8 1/4','O':'1 3/4','P':'9 1/4','Q':'5','R':'2 3/4',
    'S':'1 5/8','T':'8','U':'1/4'}
R12={'A':'2 1/4','B':'11 1/16','C':'9 3/8','D':'8 5/8','E':'1 13/32','F':'1 25/32','G':'3/4','H':'3 3/8','I':'3 3/16',
     'J':'3/16','K':'10 1/8','L':'2 1/16','M':'5 17/32','N':'12 3/8','O':'2 5/8','P':'13 7/8','Q':'7 1/2','R':'4 1/8',
     'S':'2 7/16','T':'12','U':'3/8'}
R18={'A':'3 3/8','B':'16 19/32','C':'14 1/16','D':'12 15/16','E':'2 1/8','F':'2 11/16','G':'1 1/8','H':'5 1/16',
     'I':'4 25/32','J':'9/32','K':'15 3/16','L':'3 3/32','M':'8 19/64','N':'18 9/16','O':'3 15/16','P':'20 13/16',
     'Q':'11 1/4','R':'6 3/16','S':'3 21/32','T':'18','U':'9/16'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
deg=lambda v: np.degrees(np.arctan2(v[1],v[0]))
def letter_S():
    """as Series E: top arc D centred on the axis with its top on the overshoot line, bottom arc C likewise; upper-left
    bowl I, G in from the left, and lower-right bowl H, touching the right side, internally tangent to them; the spine
    is the pair of internal tangents between the bowls.  Terminals E below / F above the overshoot lines, cut along
    the radius"""
    B,C,D,E,F,G,H,I=[DIM[k] for k in 'BCDEFGHI']
    cD=np.array([B/2,TOP-D]); cC=np.array([B/2,BOT+C])
    cI=np.array([G+I,cD[1]+np.sqrt((D-I)**2-(B/2-G-I)**2)])
    cH=np.array([B-H,cC[1]-np.sqrt((C-H)**2-(B-H-B/2)**2)])
    yt=TOP-E; xt=cD[0]+np.sqrt(D*D-(yt-cD[1])**2); a_t=deg(np.array([xt,yt])-cD)
    yb=BOT+F; xb=cC[0]-np.sqrt(C*C-(yb-cC[1])**2); a_b=deg(np.array([xb,yb])-cC)
    a_DI=deg(cI-cD); a_CH=deg(cH-cC)
    T1,T2,n1=internal_tangent(cI,I,cH,H-A,lambda s: min(s,key=lambda v:v[1]))
    T3,T4,n2=internal_tangent(cH,H,cI,I-A,lambda s: max(s,key=lambda v:v[1]))
    P=[arc(cD,D,a_t,a_DI,True), arc(cI,I,a_DI,deg(n1),True), [T2], arc(cH,H-A,deg(-n1),a_CH,False),
       arc(cC,C-A,a_CH,a_b,False), arc(cC,C,a_b,a_CH,True), arc(cH,H,a_CH,deg(n2),True), [T4],
       arc(cI,I-A,deg(-n2),a_DI,False), arc(cD,D-A,a_DI,a_t,False)]
    poly=np.vstack([np.atleast_2d(p) for p in P])
    return Region(Polygon(poly).buffer(0)),{'upper tip x':round(xt,1),'B-G':round(B-G,1),'lower tip x':round(xb,1)}
def letter_T():
    g=Grid(); K,L=DIM['K'],DIM['L']
    return rect(g,0,K,CAP-L,CAP)|rect(g,K/2-A/2,K/2+A/2,0,CAP)
def letter_U():
    g=Grid(); B,M=DIM['B'],DIM['M']; c=np.array([B/2,BOT+M])
    return rect(g,0,A,c[1],CAP)|rect(g,B-A,B,c[1],CAP)|(ring(g,c,M,M-A)&(g.Y<=c[1])),{'M vs B/2':(round(M,1),round(B/2,1))}
def letter_V():
    g=Grid(); M,N=DIM['N'],DIM['O']
    l0,l1=np.array([0,CAP]),np.array([(M-N)/2,0.0]); r0,r1=np.array([M,CAP]),np.array([(M+N)/2,0.0])
    ul=(l1-l0)/np.linalg.norm(l1-l0); ur=(r1-r0)/np.linalg.norm(r1-r0)
    m=(strip(g,l0,ul,A,(M/2,CAP))&(g.X<=M/2))|(strip(g,r0,ur,A,(M/2,CAP))&(g.X>=M/2))
    return m&(g.Y>=0)&(g.Y<=CAP)
def letter_W():
    """feet S wide, R apart in the middle; Q: height of the middle counter's apex"""
    g=Grid(); O,K,F,P=DIM['P'],DIM['S'],DIM['R'],DIM['Q']
    x1=(O-2*K-F)/2
    s1=strip(g,np.array([0,CAP]),(np.array([x1,0])-np.array([0,CAP]))/np.hypot(x1,CAP),A,(O/2,CAP))
    a2=np.array([x1+K,0.0]); ap=np.array([O/2,P])
    u2=(ap-a2)/np.linalg.norm(ap-a2); a3=np.array([O-x1-K,0.0]); u3=np.array([-u2[0],u2[1]])
    s2=strip(g,a2,u2,A,(0,CAP))&(g.X<=O/2)
    s3=strip(g,a3,u3,A,(O,CAP))&(g.X>=O/2)
    s4=strip(g,np.array([O,CAP]),(np.array([O-x1,0])-np.array([O,CAP]))/np.hypot(x1,CAP),A,(O/2,CAP))
    m=(s1&(g.X<=x1+K))|s2|s3|(s4&(g.X>=O-x1-K))
    inside=strip(g,np.array([0,CAP]),np.array([x1,-CAP])/np.hypot(x1,CAP),1e4,(O/2,CAP))
    inside&=strip(g,np.array([O,CAP]),np.array([-x1,-CAP])/np.hypot(x1,CAP),1e4,(O/2,CAP))
    return m&inside&(g.Y>=0)&(g.Y<=CAP)
def letter_X():
    g=Grid(); B,I=DIM['T'],DIM['U']
    p=np.array([I,CAP]); q=np.array([B,0.0])
    th=cross_solve(p,q,A,-np.pi/2+0.01,-0.01); u=np.array([np.cos(th),np.sin(th)])
    s1=strip(g,p,u,A,q)
    p2=np.array([B-I,CAP]); q2=np.array([0.0,0]); u2=np.array([-u[0],u[1]])
    return (s1|strip(g,p2,u2,A,q2))&(g.Y>=0)&(g.Y<=CAP),{'angle':round(np.degrees(th),1)}
