"""Series E, page E-4 (S T U V W X) - S deferred."""
import numpy as np
from fractions import Fraction as Fr
from dim1961 import Grid
from dimE3 import strip, cross_solve
def inch(s): return float(sum(Fr(p) for p in s.split()))
R8={'A':'1 3/8','B':'6 3/8','C':'4','D':'2 1/8','E':'3 7/8','F':'1 9/16','G':'1 1/4','H':'3/8','I':'2 5/16','J':'1/8',
    'K':'5 7/8','L':'3 3/16','M':'7 1/4','N':'1 1/2','O':'8 3/8','P':'2 1/2','Q':'5','R':'6 7/8','S':'1/4'}
DIM={k:inch(v)*690/8 for k,v in R8.items()}
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
ring=lambda g,c,ro,ri: g.disk(c,ro)&~g.disk(c,ri)
def letter_T():
    g=Grid(); K,G=DIM['K'],DIM['G']
    return g.bitmap(((g.X>=0)&(g.X<=K)&(g.Y>=CAP-G)&(g.Y<=CAP))|((g.X>=K/2-A/2)&(g.X<=K/2+A/2)&(g.Y>=0)&(g.Y<=CAP)))
def letter_U():
    g=Grid(); B,L=DIM['B'],DIM['L']; c=np.array([B/2,BOT+L])
    m=((g.X>=0)&(g.X<=A)|(g.X>=B-A)&(g.X<=B))&(g.Y>=c[1])&(g.Y<=CAP)
    m|=ring(g,c,L,L-A)&(g.Y<=c[1])
    return g.bitmap(m),{'L vs B/2':(round(L,1),round(B/2,1))}
def letter_V():
    g=Grid(); M,N=DIM['M'],DIM['N']
    l0,l1=np.array([0,CAP]),np.array([(M-N)/2,0.0]); r0,r1=np.array([M,CAP]),np.array([(M+N)/2,0.0])
    ul=(l1-l0)/np.linalg.norm(l1-l0); ur=(r1-r0)/np.linalg.norm(r1-r0)
    m=(strip(g,l0,ul,A,(M/2,CAP))&(g.X<=M/2))|(strip(g,r0,ur,A,(M/2,CAP))&(g.X>=M/2))
    return g.bitmap(m&(g.Y>=0)&(g.Y<=CAP))
def letter_W():
    g=Grid(); O,N,P,Q=DIM['O'],DIM['N'],DIM['P'],DIM['Q']
    x1=(O-2*N-P)/2
    s1=strip(g,np.array([0,CAP]),(np.array([x1,0])-np.array([0,CAP]))/np.hypot(x1,CAP),A,(O/2,CAP))
    a2=np.array([x1+N,0.0]); ap=np.array([O/2,Q])          # Q: top of the middle counter (the inner apex)
    u2=(ap-a2)/np.linalg.norm(ap-a2); a3=np.array([O-x1-N,0.0]); u3=np.array([-u2[0],u2[1]])
    s2=strip(g,a2,u2,A,(0,CAP))&(g.X<=O/2)
    s3=strip(g,a3,u3,A,(O,CAP))&(g.X>=O/2)
    s4=strip(g,np.array([O,CAP]),(np.array([O-x1,0])-np.array([O,CAP]))/np.hypot(x1,CAP),A,(O/2,CAP))
    m=(s1&(g.X<=x1+N))|s2|s3|(s4&(g.X>=O-x1-N))
    return g.bitmap(m&(g.Y>=0)&(g.Y<=CAP))
def letter_X():
    g=Grid(); R,S=DIM['R'],DIM['S']
    p=np.array([S,CAP]); q=np.array([R,0.0])                       # one stroke: top-left corner (S,cap) to bottom-right corner (R,0), A wide
    th=cross_solve(p,q,A,-np.pi/2+0.01,-0.01); u=np.array([np.cos(th),np.sin(th)])
    s1=strip(g,p,u,A,q)
    p2=np.array([R-S,CAP]); q2=np.array([0.0,0]); u2=np.array([-u[0],u[1]])
    s2=strip(g,p2,u2,A,q2)
    return g.bitmap((s1|s2)&(g.Y>=0)&(g.Y<=CAP)),{'angle':round(np.degrees(th),1)}

def arc(c,r,a0,a1,ccw,step=1.0):
    a0=np.radians(a0); a1=np.radians(a1)
    if ccw:
        while a1<a0: a1+=2*np.pi
    else:
        while a1>a0: a1-=2*np.pi
    n=max(2,int(abs(a1-a0)*r/step)); t=np.linspace(a0,a1,n)
    return np.column_stack([c[0]+r*np.cos(t),c[1]+r*np.sin(t)])
def internal_tangent(c1,r1,c2,r2,pick):
    w=c2-c1; d=np.linalg.norm(w); al=np.arctan2(w[1],w[0]); be=np.arccos((r1+r2)/d)
    sols=[np.array([np.cos(al+s*be),np.sin(al+s*be)]) for s in (1,-1)]
    n=pick(sols); return c1+r1*n, c2-r2*n, n
def letter_S():
    B,C,D,E,F,G,H,I=[DIM[k] for k in 'BCDEFGHI']
    cE=np.array([B/2,TOP-E]); cC=np.array([B/2,BOT+C])
    cD=np.array([H+D, cE[1]+np.sqrt((E-D)**2-(B/2-H-D)**2)])          # upper-left bowl: H in from the left, tangent to the top arc
    cI=np.array([B-I, cC[1]-np.sqrt((C-I)**2-(B-I-B/2)**2)])          # lower-right bowl: touches the right edge, tangent to the bottom arc
    deg=lambda v: np.degrees(np.arctan2(v[1],v[0]))
    yt=TOP-G; xt=cE[0]+np.sqrt(E*E-(yt-cE[1])**2); a_t=deg(np.array([xt,yt])-cE)
    yb=BOT+F; xb=cC[0]-np.sqrt(C*C-(yb-cC[1])**2); a_b=deg(np.array([xb,yb])-cC)
    a_ED=deg(cD-cE); a_CI=deg(cI-cC)
    T1,T2,n1=internal_tangent(cD,D,cI,I-A,lambda s: min(s,key=lambda v:v[1]))     # spine, lower-left edge
    T3,T4,n2=internal_tangent(cI,I,cD,D-A,lambda s: max(s,key=lambda v:v[1]))     # spine, upper-right edge
    P=[arc(cE,E,a_t,a_ED,True), arc(cD,D,a_ED,deg(n1),True), [T2], arc(cI,I-A,deg(-n1),a_CI,False),
       arc(cC,C-A,a_CI,a_b,False), arc(cC,C,a_b,a_CI,True), arc(cI,I,a_CI,deg(n2),True), [T4],
       arc(cD,D-A,deg(-n2),a_ED,False), arc(cE,E-A,a_ED,a_t,False)]
    poly=np.vstack([np.atleast_2d(p) for p in P])
    from dim1961 import raster
    return raster([[tuple(p) for p in poly]]),{'cD':cD.round(1),'cI':cI.round(1),'top tip x':round(xt,1),'bottom tip x':round(xb,1)}
