"""Series E, page E-1 (A B C D E F) from the 1961 dimensioned drawings."""
import numpy as np
from fractions import Fraction as Fr
from dim1961 import Grid
def inch(s): return float(sum(Fr(p) for p in s.split()))
R8={'A':'1 3/8','B':'8','C':'1 7/8','D':'6 3/8','E':'2 3/8','F':'2 1/8','G':'1','H':'1 1/4','I':'1/4','J':'4 1/4','K':'5 7/8',
    'L':'3 5/16','M':'1 3/4','N':'1/8','O':'3 3/4','P':'5 3/4','Q':'3 5/8','R':'3 1/2','S':'2 1/2'}
R18={'A':'3 3/32','B':'18','C':'4 7/32','D':'14 11/32','E':'5 11/32','F':'4 25/32','G':'2 1/4','H':'2 13/16','I':'9/16','J':'9 9/16',
     'K':'13 7/32','L':'7 29/64','M':'3 15/16','N':'9/32','O':'8 7/16','P':'12 15/16','Q':'8 5/32','R':'7 7/8','S':'5 5/8'}
DIM={k:inch(v)*690/18 for k,v in R18.items()}
CHECK={k:round((inch(R18[k])/18-inch(R8[k])/8)*690,2) for k in R18}
A=DIM['A']; CAP=690.0; N=DIM['N']; TOP=CAP+N; BOT=-N
ring=lambda g,c,ro,ri: g.disk(c,ro)&~g.disk(c,ri)
def letter_A():
    g=Grid(); B,C=DIM['B'],DIM['C']; xl,xr=(B-A)/2,(B+A)/2
    def edge(p0,p1): d=np.array(p1)-np.array(p0); return np.array(p0,float),d/np.linalg.norm(d)
    def side(p,d,X,Y): return (X-p[0])*d[1]-(Y-p[1])*d[0]          # >0 right of the direction
    pL,dL=edge((0,0),(xl,CAP)); pR,dR=edge((B,0),(xr,CAP))
    outer=(side(pL,dL,g.X,g.Y)>=0)&(side(pR,dR,g.X,g.Y)<=0)&(g.Y>=0)&(g.Y<=CAP)
    # inner edges: A in from the outer edges (perpendicular)
    nL=np.array([dL[1],-dL[0]]); nR=np.array([-dR[1],dR[0]])
    inner=(side(pL+nL*A,dL,g.X,g.Y)>=0)&(side(pR+nR*A,dR,g.X,g.Y)<=0)
    m=outer&~(inner&((g.Y<C)|(g.Y>C+A)))
    return g.bitmap(m)
KU=None; KL=None   # corner radii of the upper and lower counters' ends (None: fitted, BK below); the table gives none
def letter_B():
    """counters: bars A thick, the middle bar centred on J; each counter's end a side arc, G (upper) or H (lower),
    centred at the counter's mid-height, joined to its flat top and bottom by corner arcs tangent to both (radius BK:
    the table has no corner dimension, so it is fitted to the drawing, which shows the ends running smoothly into the
    flats).  (A counter is shorter than its side arc's diameter, so the arc alone met the flats at about 30 degrees.)"""
    g=Grid(); D,E,F,G,H,I,J=DIM['D'],DIM['E'],DIM['F'],DIM['G'],DIM['H'],DIM['I'],DIM['J']
    cu=np.array([D-I-F,CAP-F]); cl=np.array([D-E,E])                     # outer bowls touch the top / baseline
    outer=((g.X>=0)&(g.X<=cu[0])&(g.Y>=J)&(g.Y<=CAP))|((g.X>=0)&(g.X<=cl[0])&(g.Y>=0)&(g.Y<=J))
    outer|=(g.disk(cu,F)&(g.X>=cu[0])&(g.Y>=J-200))|(g.disk(cl,E)&(g.X>=cl[0])&(g.Y<=J+200))
    outer&=(g.Y>=0)&(g.Y<=CAP)
    ut,ub=CAP-A,J+A/2; lt,lb=J-A/2,A
    ciu=np.array([D-I-A-G,(ut+ub)/2]); cil=np.array([D-A-H,(lt+lb)/2])
    def end(top,bot,c,R,k):
        """counter from x = A rightwards: side arc R centred at c, corners k tangent to the flats and to R"""
        def corner(y,sgn):
            yc=y-sgn*k; return np.array([c[0]+np.sqrt((R-k)**2-(yc-c[1])**2),yc])
        ct,cb=corner(top,1),corner(bot,-1)
        ang=lambda p,q: np.degrees(np.arctan2(q[1]-p[1],q[0]-p[0]))
        m=(g.X>=A)&(g.Y>bot)&(g.Y<top)&(g.disk(c,R)|(g.X<=c[0]))
        m&=~(g.sector(ct,ang(c,ct),90)&~g.disk(ct,k))
        m&=~(g.sector(cb,-90,ang(c,cb))&~g.disk(cb,k))
        return m
    cnt=end(ut,ub,ciu,G,BK[0] if KU is None else KU)|end(lt,lb,cil,H,BK[1] if KL is None else KL)
    return g.bitmap(outer&~cnt)
BK=(70.0,85.0)   # fitted to the drawing (page 180): upper counter half-height 72.8, lower 94.3
def letter_C():
    g=Grid(); D,K,L,M=DIM['D'],DIM['K'],DIM['L'],DIM['M']
    cK=np.array([K,CAP/2]); yu=TOP-L; yl=BOT+L
    xL=cK[0]-np.sqrt((K-L)**2-(yu-cK[1])**2); cu=np.array([xL,yu]); cl=np.array([xL,yl])
    a_ku=np.degrees(np.arctan2(cu[1]-cK[1],cu[0]-cK[0])); a_kl=np.degrees(np.arctan2(cl[1]-cK[1],cl[0]-cK[0]))
    yt=TOP-M; xt=cu[0]+np.sqrt(L*L-(yt-cu[1])**2); a_t=np.degrees(np.arctan2(yt-cu[1],xt-cu[0]))
    yb=BOT+M; xb=cl[0]+np.sqrt(L*L-(yb-cl[1])**2); a_b=np.degrees(np.arctan2(yb-cl[1],xb-cl[0]))
    m=ring(g,cK,K,K-A)&g.sector(cK,a_ku,a_kl)
    m|=ring(g,cu,L,L-A)&g.sector(cu,a_t,a_ku)
    m|=ring(g,cl,L,L-A)&g.sector(cl,a_kl,a_b)
    return g.bitmap(m),{'terminal tip x':round(xt,1),'D':round(D,1)}
def letter_D():
    g=Grid(); D,K,O=DIM['D'],DIM['K'],DIM['O']
    cK=np.array([D-K,CAP/2]); yu=CAP-O; yl=O
    xO=cK[0]+np.sqrt((K-O)**2-(yu-cK[1])**2); cu=np.array([xO,yu]); cl=np.array([xO,yl])
    a_u=np.degrees(np.arctan2(cu[1]-cK[1],cu[0]-cK[0])); a_l=np.degrees(np.arctan2(cl[1]-cK[1],cl[0]-cK[0]))
    outer=((g.X>=0)&(g.X<=xO)&(g.Y>=0)&(g.Y<=CAP))|(g.disk(cK,K)&g.sector(cK,a_l,a_u))|(g.disk(cu,O)&(g.X>=xO)&g.sector(cu,0,90))|(g.disk(cl,O)&(g.X>=xO)&g.sector(cl,-90,0))
    outer|=(g.disk(cK,K)&(g.X>=xO))&(g.Y<=cu[1])&(g.Y>=cl[1])
    inner=((g.X>=A)&(g.X<=xO)&(g.Y>=A)&(g.Y<=CAP-A))|(g.disk(cK,K-A)&g.sector(cK,a_l,a_u)&(g.X>=A))|(g.disk(cu,O-A)&(g.X>=xO)&g.sector(cu,0,90))|(g.disk(cl,O-A)&(g.X>=xO)&g.sector(cl,-90,0))
    inner|=(g.disk(cK,K-A)&(g.X>=xO))&(g.Y<=cu[1])&(g.Y>=cl[1])
    return g.bitmap(outer&~inner)
def letter_E():
    g=Grid(); P,E,K,Q=DIM['P'],DIM['E'],DIM['K'],DIM['Q']
    m=((g.X>=0)&(g.X<=A)&(g.Y>=0)&(g.Y<=CAP))|((g.X>=0)&(g.X<=P)&(g.Y>=CAP-A)&(g.Y<=CAP))
    m|=(g.X>=0)&(g.X<=A+E)&(g.Y>=Q)&(g.Y<=Q+A); m|=(g.X>=0)&(g.X<=K)&(g.Y>=0)&(g.Y<=A)
    return g.bitmap(m)
def letter_F():
    g=Grid(); K,S,R=DIM['K'],DIM['S'],DIM['R']
    m=((g.X>=0)&(g.X<=A)&(g.Y>=0)&(g.Y<=CAP))|((g.X>=0)&(g.X<=K)&(g.Y>=CAP-A)&(g.Y<=CAP))
    m|=(g.X>=0)&(g.X<=A+S)&(g.Y>=R)&(g.Y<=R+A)
    return g.bitmap(m)
