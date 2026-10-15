"""Series D, page D-1 (A B C D E F) from the 1961 dimensioned drawings.  Written directly against the vector Grid."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
R8={'A':'1 1/4','B':'6 11/16','C':'2 1/4','D':'5 3/8','E':'3 1/2','F':'3/16','G':'4 3/16','H':'2','I':'2 3/16','J':'1 5/16',
    'K':'1 1/4','L':'1/4','M':'3 1/2','N':'2 5/8','O':'2 9/16','P':'1 5/8','Q':'1 3/4','R':'6 5/16','S':'1/8','T':'2 7/8',
    'U':'5 5/8','V':'2 5/16','W':'4 7/8'}
R12={'A':'1 7/8','B':'10 1/32','C':'3 3/8','D':'8 1/16','E':'5 1/4','F':'9/32','G':'6 9/32','H':'3','I':'3 9/32','J':'1 31/32',
     'K':'1 7/8','L':'3/8','M':'5 1/4','N':'3 15/16','O':'3 27/32','P':'2 7/16','Q':'2 5/8','R':'9 15/32','S':'3/16','T':'4 5/16',
     'U':'8 7/16','V':'3 15/32','W':'7 5/16'}
R18={'A':'2 13/16','B':'15 1/16','C':'5 1/16','D':'12 3/32','E':'7 7/8','F':'7/16','G':'9 7/16','H':'4 1/2','I':'4 15/16',
     'J':'2 15/16','K':'2 13/16','L':'9/16','M':'7 7/8','N':'5 29/32','O':'5 3/4','P':'3 21/32','Q':'3 15/16','R':'14 3/16',
     'S':'9/32','T':'6 15/32','U':'12 21/32','V':'5 3/16','W':'10 31/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['S']; TOP=CAP+OV; BOT=-OV
ring=lambda g,c,ro,ri: g.disk(c,ro)&~g.disk(c,ri)
ang=lambda c,p: np.degrees(np.arctan2(p[1]-c[1],p[0]-c[0]))
def rect(g,x0,x1,y0,y1): return (g.X>=x0)&(g.X<=x1)&(g.Y>=y0)&(g.Y<=y1)
def letter_A():
    g=Grid(); B,C=DIM['B'],DIM['C']; xl,xr=(B-A)/2,(B+A)/2
    def edge(p0,p1): d=np.array(p1,float)-np.array(p0,float); return np.array(p0,float),d/np.linalg.norm(d)
    def side(p,d): return (g.X-p[0])*d[1]-(g.Y-p[1])*d[0]
    pL,dL=edge((0,0),(xl,CAP)); pR,dR=edge((B,0),(xr,CAP))
    outer=(side(pL,dL)>=0)&(side(pR,dR)<=0)&(g.Y>=0)&(g.Y<=CAP)
    nL=np.array([dL[1],-dL[0]]); nR=np.array([-dR[1],dR[0]])
    pi,pj=pL+nL*A,pR+nR*A
    inner=((g.X-pi[0])*dL[1]-(g.Y-pi[1])*dL[0]>=0)&((g.X-pj[0])*dR[1]-(g.Y-pj[1])*dR[0]<=0)
    return outer&~(inner&((g.Y<C)|(g.Y>C+A)))
def letter_B():
    """outer bowls H (upper) and I (lower) centred on one vertical, tangent to the cap line and the baseline; they cross at
    the crotch (height G).  Counters end in arcs K (upper) and J (lower), A in from the outer extremes (upper one F in)."""
    g=Grid(); D,E,F,G,H,I,J,K=[DIM[k] for k in 'DEFGHIJK']
    cu=np.array([D-F-H,CAP-H]); cl=np.array([D-I,I])
    # crotch: where the two outer circles cross
    dy=cu[1]-cl[1]; yi=cl[1]+(dy*dy+I*I-H*H)/(2*dy)
    outer=rect(g,0,cu[0],0,CAP)|(g.disk(cu,H)&(g.X>=cu[0])&(g.Y>=yi))|(g.disk(cl,I)&(g.X>=cl[0])&(g.Y<=yi))
    ut,ub=CAP-A,E+A; lt,lb=E,A
    # counter arcs are centred level with their outer arcs' centres, which makes them tangent to the middle bar
    # (CAP-H-K = E+A, I+J = E); they meet the top and bottom bars at a slight angle, as drawn
    ku=np.array([D-F-A-K,cu[1]]); kl=np.array([D-A-J,cl[1]])
    cnt=(rect(g,A,ku[0],ub,ut)|(g.disk(ku,K)&(g.X>=ku[0])))&(g.Y>=ub)&(g.Y<=ut)
    cnt|=(rect(g,A,kl[0],lb,lt)|(g.disk(kl,J)&(g.X>=kl[0])))&(g.Y>=lb)&(g.Y<=lt)
    return outer&~cnt,{'crotch y':round(yi,1),'G':round(G,1),'centres x':(round(cu[0],1),round(cl[0],1))}
def letter_C():
    """big arc R centred at mid-height (left side touches x=0); upper arc N and lower arc O touch the overshoot lines and
    are internally tangent to R; terminals at P below / Q above the overshoot lines, cut along the radius"""
    g=Grid(); D,N,O,P,Q,R=[DIM[k] for k in 'DNOPQR']
    cR=np.array([R,CAP/2])
    yn=TOP-N; cN=np.array([R-np.sqrt((R-N)**2-(yn-cR[1])**2),yn])
    yo=BOT+O; cO=np.array([R-np.sqrt((R-O)**2-(yo-cR[1])**2),yo])
    yt=TOP-P; xt=cN[0]+np.sqrt(N*N-(yt-cN[1])**2)
    yb=BOT+Q; xb=cO[0]+np.sqrt(O*O-(yb-cO[1])**2)
    m=ring(g,cR,R,R-A)&g.sector(cR,ang(cR,cN),ang(cR,cO))
    m|=ring(g,cN,N,N-A)&g.sector(cN,ang(cN,(xt,yt)),ang(cR,cN))
    m|=ring(g,cO,O,O-A)&g.sector(cO,ang(cR,cO),ang(cO,(xb,yb)))
    return m,{'upper tip x':round(xt,1),'lower tip x':round(xb,1),'D':round(D,1)}
def letter_D():
    """big arc U centred at mid-height (D-U from the left), corner arcs T tangent to it and to the top/bottom lines"""
    g=Grid(); D,T,U=DIM['D'],DIM['T'],DIM['U']
    cU=np.array([D-U,CAP/2]); yu=CAP-T; yl=T
    xT=cU[0]+np.sqrt((U-T)**2-(yu-cU[1])**2); cu=np.array([xT,yu]); cl=np.array([xT,yl])
    a_u,a_l=ang(cU,cu),ang(cU,cl)
    def body(d):
        r=rect(g,d,xT,d,CAP-d)|(g.disk(cU,U-d)&g.sector(cU,a_l,a_u)&(g.X>=xT))
        r|=g.disk(cu,T-d)&(g.X>=xT)&g.sector(cu,a_u,90)
        r|=g.disk(cl,T-d)&(g.X>=xT)&g.sector(cl,-90,a_l)
        return r
    return body(0)&~body(A)
def letter_E():
    g=Grid(); L,M,V,W=DIM['L'],DIM['M'],DIM['V'],DIM['W']
    return rect(g,0,A,0,CAP)|rect(g,0,W-L,CAP-A,CAP)|rect(g,0,A+V,M,M+A)|rect(g,0,W,0,A)
def letter_F():
    g=Grid(); M,V,W=DIM['M'],DIM['V'],DIM['W']
    return rect(g,0,A,0,CAP)|rect(g,0,W,CAP-A,CAP)|rect(g,0,A+V,M,M+A)
