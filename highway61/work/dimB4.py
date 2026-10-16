"""Series B, page B-4 (S T U V W X)."""
import numpy as np
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip, cross_solve
from dimE4 import arc
R8={'A':'1','B':'3 3/8','C':'3 3/16','D':'4 5/8','E':'1 19/32','F':'1 11/16','G':'2 1/8','H':'1 7/8','I':'1/8',
    'J':'3/16','K':'3','L':'3 11/16','M':'7/8','N':'5','O':'4 1/2','P':'1 1/2','Q':'1/16'}
R12={'A':'1 1/2','B':'5 1/16','C':'4 25/32','D':'6 15/16','E':'2 25/64','F':'2 17/32','G':'3 3/16','H':'2 13/16',
     'I':'3/16','J':'9/32','K':'4 1/2','L':'5 17/32','M':'1 5/16','N':'7 1/2','O':'6 3/4','P':'2 1/4','Q':'3/32'}
R18={'A':'2 1/4','B':'7 19/32','C':'7 5/32','D':'10 13/32','E':'3 37/64','F':'3 51/64','G':'4 25/32','H':'4 7/32',
     'I':'9/32','J':'7/16','K':'6 3/4','L':'8 5/16','M':'1 31/32','N':'11 1/4','O':'10 1/8','P':'3 3/8','Q':'1/8'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['I']; TOP=CAP+OV; BOT=-OV
def letter_S():
    """upper: top arc E (on the overshoot line) and left arc C on one centre line, touching at the left extreme J in
    from the left; lower: bottom arc F (on the overshoot line) and right arc D on one centre line, touching the right
    side.  The spine is C running into D-A (and C-A into D): the table's C leaves them overlapping by 1.6 units, an
    8.6 degree kink, so C is refitted to make them tangent.  Terminals cut horizontally, carried straight to H below the
    top overshoot line (upper, on the right) and G above the bottom one (lower, on the left)."""
    B,D,E,F,G,H,J=[DIM[k] for k in 'BDEFGHJ']
    yu=TOP-E; yl=BOT+F
    cE=np.array([J+E,yu]); cF=np.array([B-F,yl]); cD=np.array([B-D,yl])
    dy=yu-yl; k=J+D-B                                          # C's centre x = J+C; distance to cD = C+D-A
    C=((D-A)**2-k*k-dy*dy)/(2*(k-(D-A)))
    cC=np.array([J+C,yu])
    w=cD-cC; Tp=cC+C*w/np.linalg.norm(w); Tq=cC+(C-A)*w/np.linalg.norm(w)   # touching points of the spine edges
    a_CD=ang(cC,cD); a_DC=ang(cD,cC)
    yt=TOP-H; yb=BOT+G
    P=[[(cE[0]+E,yt)],arc(cE,E,0,180,True),arc(cC,C,180,a_CD,True),arc(cD,D-A,a_DC,0,False),arc(cF,F-A,0,-180,False),
       [(cF[0]-F+A,yb),(cF[0]-F,yb)],arc(cF,F,180,360,True),arc(cD,D,0,a_DC,True),arc(cC,C-A,a_CD,180,False),
       arc(cE,E-A,180,0,False),[(cE[0]+E-A,yt)]]
    poly=np.vstack([np.atleast_2d(p) for p in P])
    return Region(Polygon(poly).buffer(0)),{'C refit':round(C,1),'C table':round(DIM['C'],1),
        'upper terminal right x':round(cE[0]+E,1),'lower left x':round(cF[0]-F,1),'B':round(B,1)}
def letter_T():
    g=Grid(); K=DIM['K']
    return rect(g,0,K,CAP-A,CAP)|rect(g,K/2-A/2,K/2+A/2,0,CAP)
def letter_U():
    g=Grid(); B,F=DIM['B'],DIM['F']; c=np.array([B/2,BOT+F])
    return rect(g,0,A,c[1],CAP)|rect(g,B-A,B,c[1],CAP)|(ring(g,c,F,F-A)&(g.Y<=c[1])),{'F vs B/2':(round(F,1),round(B/2,1))}
def letter_V():
    g=Grid(); M,N=DIM['L'],DIM['M']
    l0,l1=np.array([0,CAP]),np.array([(M-N)/2,0.0]); r0,r1=np.array([M,CAP]),np.array([(M+N)/2,0.0])
    ul=(l1-l0)/np.linalg.norm(l1-l0); ur=(r1-r0)/np.linalg.norm(r1-r0)
    m=(strip(g,l0,ul,A,(M/2,CAP))&(g.X<=M/2))|(strip(g,r0,ur,A,(M/2,CAP))&(g.X>=M/2))
    return m&(g.Y>=0)&(g.Y<=CAP)
def letter_W():
    """feet M wide, P apart in the middle; O: height of the middle counter's apex"""
    g=Grid(); O,K,F,P=DIM['N'],DIM['M'],DIM['P'],DIM['O']
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
    g=Grid(); B,I=DIM['L'],DIM['Q']
    p=np.array([I,CAP]); q=np.array([B,0.0])
    th=cross_solve(p,q,A,-np.pi/2+0.01,-0.01); u=np.array([np.cos(th),np.sin(th)])
    s1=strip(g,p,u,A,q)
    p2=np.array([B-I,CAP]); q2=np.array([0.0,0]); u2=np.array([-u[0],u[1]])
    return (s1|strip(g,p2,u2,A,q2))&(g.Y>=0)&(g.Y<=CAP),{'angle':round(np.degrees(th),1)}
