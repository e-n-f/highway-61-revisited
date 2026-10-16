"""Series A, page A-4 (S T U V W X).  The constructions of Series B (dimB4) with this page's letters."""
import numpy as np
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip, cross_solve
from dimE4 import arc
R8={'A':'3/4','B':'3','C':'4','D':'4 3/4','E':'2','F':'1 3/4','G':'1 1/2','H':'1 13/32','I':'1/8','J':'3/16',
    'K':'2 11/16','L':'11/16','M':'3 1/2','N':'4 3/8','O':'1','P':'4 1/4','Q':'1/16'}
R12={'A':'1 1/8','B':'4 1/2','C':'6','D':'7 1/8','E':'3','F':'2 5/8','G':'2 1/4','H':'2 7/64','I':'3/16','J':'9/32',
     'K':'4 1/32','L':'1 1/32','M':'5 1/4','N':'6 9/16','O':'1 1/2','P':'6 3/8','Q':'3/32'}
R18={'A':'1 11/16','B':'6 3/4','C':'9','D':'10 11/16','E':'4 1/2','F':'3 15/16','G':'3 3/8','H':'3 5/32','I':'9/32',
     'J':'7/16','K':'6 1/16','L':'1 9/16','M':'7 7/8','N':'9 27/32','O':'2 1/4','P':'9 9/16','Q':'1/8'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['I']; TOP=CAP+OV; BOT=-OV
def letter_S():
    """upper: top arc H (on the overshoot line) and the big arc C on one centre line, touching at the left extreme J in
    from the left; lower: bottom arc G (on the overshoot line) and the big arc D on one centre line, touching the right
    side (B).  The spine is C running into D-A (and C-A into D), C refitted to make them tangent if the table does
    not.  Terminals cut horizontally, carried straight to F below the top overshoot line (upper, on the right) and E
    above the bottom one (lower, on the left)."""
    B,D,E,F,G,H,J=[DIM[k] for k in 'BDEFGHJ']
    yu=TOP-H; yl=BOT+G
    cH=np.array([J+H,yu]); cG=np.array([B-G,yl]); cD=np.array([B-D,yl])
    dy=yu-yl; k=J+D-B                                          # C's centre x = J+C; distance to cD = C+D-A
    C=((D-A)**2-k*k-dy*dy)/(2*(k-(D-A)))
    cC=np.array([J+C,yu])
    a_CD=ang(cC,cD); a_DC=ang(cD,cC)
    yt=TOP-F; yb=BOT+E
    P=[[(cH[0]+H,yt)],arc(cH,H,0,180,True),arc(cC,C,180,a_CD,True),arc(cD,D-A,a_DC,0,False),arc(cG,G-A,0,-180,False),
       [(cG[0]-G+A,yb),(cG[0]-G,yb)],arc(cG,G,180,360,True),arc(cD,D,0,a_DC,True),arc(cC,C-A,a_CD,180,False),
       arc(cH,H-A,180,0,False),[(cH[0]+H-A,yt)]]
    poly=np.vstack([np.atleast_2d(p) for p in P])
    return Region(Polygon(poly).buffer(0)),{'C refit':round(C,1),'C table':round(DIM['C'],1),
        'upper terminal right x':round(cH[0]+H,1),'lower left x':round(cG[0]-G,1),'B':round(B,1)}
def letter_T():
    """bar K wide and L thick"""
    g=Grid(); K,L=DIM['K'],DIM['L']
    return rect(g,0,K,CAP-L,CAP)|rect(g,K/2-A/2,K/2+A/2,0,CAP)
def letter_U():
    """width B; bottom arc G centred on the axis, tangent to the bottom overshoot line"""
    g=Grid(); B,G=DIM['B'],DIM['G']; c=np.array([B/2,BOT+G])
    return rect(g,0,A,c[1],CAP)|rect(g,B-A,B,c[1],CAP)|(ring(g,c,G,G-A)&(g.Y<=c[1])),{'G vs B/2':(round(G,1),round(B/2,1))}
def letter_V():
    """width M, foot L"""
    g=Grid(); M,N=DIM['M'],DIM['L']
    l0,l1=np.array([0,CAP]),np.array([(M-N)/2,0.0]); r0,r1=np.array([M,CAP]),np.array([(M+N)/2,0.0])
    ul=(l1-l0)/np.linalg.norm(l1-l0); ur=(r1-r0)/np.linalg.norm(r1-r0)
    m=(strip(g,l0,ul,A,(M/2,CAP))&(g.X<=M/2))|(strip(g,r0,ur,A,(M/2,CAP))&(g.X>=M/2))
    return m&(g.Y>=0)&(g.Y<=CAP)
def letter_W():
    """width N; feet L wide, O apart in the middle; P: height of the middle counter's apex"""
    g=Grid(); O,K,F,P=DIM['N'],DIM['L'],DIM['O'],DIM['P']
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
    """width B; the strokes' outer edges start Q in from the top corners"""
    g=Grid(); B,I=DIM['B'],DIM['Q']
    p=np.array([I,CAP]); q=np.array([B,0.0])
    th=cross_solve(p,q,A,-np.pi/2+0.01,-0.01); u=np.array([np.cos(th),np.sin(th)])
    s1=strip(g,p,u,A,q)
    p2=np.array([B-I,CAP]); q2=np.array([0.0,0]); u2=np.array([-u[0],u[1]])
    return (s1|strip(g,p2,u2,A,q2))&(g.Y>=0)&(g.Y<=CAP),{'angle':round(np.degrees(th),1)}
