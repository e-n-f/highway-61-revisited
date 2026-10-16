"""Series C, page C-4 (S T U V W X)."""
import numpy as np
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip, cross_solve
from dimE4 import arc, internal_tangent
R8={'A':'1 1/8','B':'4 3/8','C':'2','D':'2 3/16','E':'3 1/4','F':'2 1/2','G':'4 7/8','H':'3/8','I':'4','J':'1/8',
    'K':'1 1/16','L':'7/8','M':'6','N':'4 3/16','O':'1','P':'1 1/2','Q':'4 11/16'}
R12={'A':'1 11/16','B':'6 9/16','C':'3','D':'3 9/32','E':'4 7/8','F':'3 3/4','G':'7 5/16','H':'9/16','I':'6','J':'3/16',
     'K':'1 19/32','L':'1 5/16','M':'9','N':'6 9/32','O':'1 1/2','P':'2 1/4','Q':'7 1/32'}
R18={'A':'2 17/32','B':'9 27/32','C':'4 1/2','D':'4 15/16','E':'7 5/16','F':'5 5/8','G':'10 31/32','H':'27/32','I':'9',
     'J':'9/32','K':'2 3/8','L':'1 31/32','M':'13 1/2','N':'9 7/16','O':'2 1/4','P':'3 3/8','Q':'10 9/16'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
deg=lambda v: np.degrees(np.arctan2(v[1],v[0]))
def letter_S():
    """upper bowl C (top on the overshoot line, H in from the left); lower bowl D (bottom on the overshoot line, touching
    the right side) and arc E on the same centre line, touching it at the right; the spine is the pair of internal
    tangents C / E-A and E / C-A.  Terminals cut horizontally: the upper at its centre height on the right, the lower
    carried straight up on the left to F above the overshoot line."""
    B,C,D,E,F,H=[DIM[k] for k in 'BCDEFH']
    cC=np.array([H+C,TOP-C]); cD=np.array([B-D,BOT+D]); cE=np.array([B-E,cD[1]])
    T1,T2,n1=internal_tangent(cC,C,cE,E-A,lambda s: min(s,key=lambda v:v[0]+v[1]))
    T3,T4,n2=internal_tangent(cE,E,cC,C-A,lambda s: max(s,key=lambda v:v[0]+v[1]))
    xl=cD[0]-D; yF=BOT+F
    P=[arc(cC,C,0,deg(n1),True),[T2],arc(cE,E-A,deg(-n1),0,False),arc(cD,D-A,0,-180,False),
       [(cD[0]-D+A,yF),(xl,yF)],arc(cD,D,180,360,True),arc(cE,E,0,deg(n2),True),[T4],
       arc(cC,C-A,deg(-n2),0,False)]
    poly=np.vstack([np.atleast_2d(p) for p in P])
    return Region(Polygon(poly).buffer(0)),{'right extremes C/D/E':(round(cC[0]+C,1),round(cD[0]+D,1),round(cE[0]+E,1)),
                                           'B':round(B,1),'lower terminal left x':round(xl,1)}
def letter_T():
    g=Grid(); I,K=DIM['I'],DIM['K']
    return rect(g,0,I,CAP-K,CAP)|rect(g,I/2-A/2,I/2+A/2,0,CAP)
def letter_U():
    g=Grid(); B,D=DIM['B'],DIM['D']; c=np.array([B/2,BOT+D])
    return rect(g,0,A,c[1],CAP)|rect(g,B-A,B,c[1],CAP)|(ring(g,c,D,D-A)&(g.Y<=c[1])),{'D vs B/2':(round(D,1),round(B/2,1))}
def letter_V():
    g=Grid(); M,N=DIM['G'],DIM['L']
    l0,l1=np.array([0,CAP]),np.array([(M-N)/2,0.0]); r0,r1=np.array([M,CAP]),np.array([(M+N)/2,0.0])
    ul=(l1-l0)/np.linalg.norm(l1-l0); ur=(r1-r0)/np.linalg.norm(r1-r0)
    m=(strip(g,l0,ul,A,(M/2,CAP))&(g.X<=M/2))|(strip(g,r0,ur,A,(M/2,CAP))&(g.X>=M/2))
    return m&(g.Y>=0)&(g.Y<=CAP)
def letter_W():
    """feet O wide, P apart in the middle; N: height of the middle counter's apex"""
    g=Grid(); O,K,F,P=DIM['M'],DIM['O'],DIM['P'],DIM['N']
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
    g=Grid(); B,I=DIM['Q'],DIM['J']
    p=np.array([I,CAP]); q=np.array([B,0.0])
    th=cross_solve(p,q,A,-np.pi/2+0.01,-0.01); u=np.array([np.cos(th),np.sin(th)])
    s1=strip(g,p,u,A,q)
    p2=np.array([B-I,CAP]); q2=np.array([0.0,0]); u2=np.array([-u[0],u[1]])
    return (s1|strip(g,p2,u2,A,q2))&(g.Y>=0)&(g.Y<=CAP),{'angle':round(np.degrees(th),1)}
