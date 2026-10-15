"""Series E capitals/figures rebuilt from the 1961 dimensioned drawings (page E-5)."""
import numpy as np,cv2,pickle
from fractions import Fraction as Fr
H_IN=18.0; U=690.0/H_IN                    # font units per inch at 18-inch letter height
# 18-inch row of the table, read from the page (inches)
T18={'A':'3 3/32','B':'18','C':'7 5/16','D':'14 11/32','E':'9/16','F':'3 21/32','G':'4 1/2','H':'2 17/32','I':'1 31/32',
     'J':'2 1/4','K':'10 11/16','L':'10 31/32','M':'5 1/16','N':'2 13/16','O':'7 5/32','P':'8 7/16','Q':'9','R':'10 1/8',
     'S':'27/32','T':'9/32','U':'16 19/32','V':'2 15/16'}
# 8-inch row, for a consistency check
T8={'A':'1 3/8','B':'8','C':'3 1/4','D':'6 3/8','E':'1/4','F':'1 5/8','G':'2','H':'1 1/8','I':'7/8','J':'1','K':'4 3/4','L':'4 7/8',
    'M':'2 1/4','N':'1 1/4','O':'3 3/16','P':'3 3/4','Q':'4','R':'4 1/2','S':'3/8','T':'1/8','U':'7 3/8','V':'1 5/16'}
def inch(s): return float(sum(Fr(p) for p in s.split()))
DIM={k:inch(v)*U for k,v in T18.items()}
def check():
    out=[]
    for k in T18:
        a=inch(T18[k])/18; b=inch(T8[k])/8; out.append((k,round(a*690,1),round(b*690,1),round((a-b)*690,1)))
    return out
A=DIM['A']; CAP=690.0
def line_inter(p,d,q,e):
    M=np.array([d,-e]).T; s,_=np.linalg.solve(M,q-p); return p+d*s
def Y():
    B,C=DIM['B'],DIM['C']; c=B/2; sl,sr=c-A/2,c+A/2
    p0=np.array([0,CAP]); p1=np.array([sl,C]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])   # n points right/up into the arm
    if n[0]<0: n=-n
    qi=p0+n*A                                                   # inner edge of left arm
    xt=line_inter(qi,d,np.array([0,CAP]),np.array([1,0]))[0]
    v=line_inter(qi,d,np.array([c,0]),np.array([0,1]))
    return [[(0,CAP),(xt,CAP),tuple(v),(B-xt,CAP),(B,CAP),(sr,C),(sr,0),(sl,0),(sl,C)]]
def Z():
    D,E,F=DIM['D'],DIM['E'],DIM['F']
    P1=np.array([D,CAP-F]); P2=np.array([0,F]); dx,dy=P1-P2
    # angle of the diagonal so that the two parallel edges through P1 and P2 are A apart (perpendicular)
    from scipy.optimize import brentq
    th=brentq(lambda t: dx*np.sin(t)-dy*np.cos(t)-A, np.arctan2(dy,dx), np.pi/2)
    u=np.array([np.cos(th),np.sin(th)])
    r_bot=line_inter(P1,u,np.array([0,A]),np.array([1,0]))     # right edge meets bottom bar top
    l_top=line_inter(P2,u,np.array([0,CAP-A]),np.array([1,0])) # left edge meets top bar underside
    return [[(E,CAP),(D,CAP),(D,CAP-F),tuple(r_bot),(D,A),(D,0),(0,0),(0,F),tuple(l_top),(E,CAP-A)]], np.degrees(th)
def one():
    G,H,I,J=DIM['G'],DIM['H'],DIM['I'],DIM['J']
    sl=J; sr=J+A; a=np.array([0,CAP-I]); b=np.array([sr-H,CAP])
    m=(a+b)/2; ch=b-a; L=np.linalg.norm(ch); h=np.sqrt(G*G-(L/2)**2); perp=np.array([-ch[1],ch[0]])/L
    cA=m+perp*h; cB=m-perp*h; c=cA if (cA[0]<m[0] and cA[1]>m[1]) else cB            # centre up and to the left (concave flag)
    a0=np.arctan2(a[1]-c[1],a[0]-c[0]); a1=np.arctan2(b[1]-c[1],b[0]-c[0])
    if a1<a0: a1+=2*np.pi
    if a1-a0>np.pi: a0,a1=a1,a0+2*np.pi
    arc=[tuple(c+G*np.array([np.cos(t),np.sin(t)])) for t in np.linspace(a0,a1,40)]
    # arc runs from a to b along the short way; order it from a to b
    if np.hypot(*(np.array(arc[0])-a))>1: arc=arc[::-1]
    return [[(0,CAP-I-H),(sl,CAP-I-H),(sl,0),(sr,0),(sr,CAP)]+arc[::-1][:-1]+[(0,CAP-I)]], c
def four():
    Uw,F,V,H=DIM['U'],DIM['F'],DIM['V'],DIM['H']
    sr=Uw-H; sl=sr-A; yb=F+V
    p0=np.array([0,yb]); p1=np.array([sl,CAP]); d=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([d[1],-d[0]])
    if n[0]<0: n=-n                                             # into the counter (down-right)
    qi=p0+n*A
    ib=line_inter(qi,d,np.array([0,yb]),np.array([1,0])); ia=line_inter(qi,d,np.array([sl,0]),np.array([0,1]))
    outer=[(0,F),(sl,F),(sl,0),(sr,0),(sr,F),(Uw,F),(Uw,yb),(sr,yb),(sr,CAP),(sl,CAP),(0,yb)]
    hole=[tuple(ib),(sl,yb),tuple(ia)]
    return [outer,hole]
def raster(polys,scale=1.0,pad=40):
    allp=np.vstack([np.array(p) for p in polys]); x0=allp[:,0].min(); y1=allp[:,1].max()
    W=int(allp[:,0].max()-x0+2*pad); Hh=int(y1-min(0,allp[:,1].min())+2*pad)
    img=np.zeros((Hh,W),np.uint8)
    pts=[np.round(np.column_stack([np.array(p)[:,0]-x0+pad, pad+y1-np.array(p)[:,1]])*16).astype(np.int32) for p in polys]
    cv2.fillPoly(img,pts[:1],1,shift=4)
    for h in pts[1:]: cv2.fillPoly(img,[h],0,shift=4)
    return img, pad, pad+y1                                     # bitmap, x of letter's left edge, row of baseline

# ---------------- curved figures: unions of circular rings (compass construction) ----------------
TT=DIM['T']; TOP=CAP+TT; BOT=-TT
class Grid:
    def __init__(s,x0=-60,x1=780,y0=-80,y1=760):
        s.x0,s.y1=x0,y1; xs=np.arange(x0,x1)+0.5; ys=y1-np.arange(0,y1-y0)-0.5
        s.X,s.Y=np.meshgrid(xs,ys)
    def disk(s,c,r): return (s.X-c[0])**2+(s.Y-c[1])**2<=r*r
    def ang(s,c): return np.degrees(np.arctan2(s.Y-c[1],s.X-c[0]))
    def sector(s,c,a0,a1):
        a=s.ang(c); a0%=360; a1%=360; a=a%360
        return (a>=a0)&(a<=a1) if a0<=a1 else (a>=a0)|(a<=a1)
    def bitmap(s,m): return m.astype(np.uint8), -s.x0, s.y1
def two():
    g=Grid(); D,K,M,L,N=DIM['D'],DIM['K'],DIM['M'],DIM['L'],DIM['N']
    xv=D-DIM['O']                                             # the vertical construction line
    cK=np.array([xv,TOP-K])                                   # top arc: top of the curve at cap height + overshoot T
    cM=np.array([D-M, cK[1]+np.sqrt((K-M)**2-(D-M-xv)**2)])  # bowl: touches the right side, internally tangent to K
    cL=np.array([L, cM[1]-np.sqrt((M+L-A)**2-(L-cM[0])**2)])  # spine: starts at (0,0), externally tangent to the bowl
    yt=TOP-N; xt=cK[0]-np.sqrt(K*K-(yt-cK[1])**2)             # terminal: N below the overshoot line, cut along the radius
    a_t=np.degrees(np.arctan2(yt-cK[1],xt-cK[0])); a_km=np.degrees(np.arctan2(cM[1]-cK[1],cM[0]-cK[0]))
    a_ml=np.degrees(np.arctan2(cL[1]-cM[1],cL[0]-cM[0])); a_lm=np.degrees(np.arctan2(cM[1]-cL[1],cM[0]-cL[0]))
    ring=lambda c,ro,ri: g.disk(c,ro)&~g.disk(c,ri)
    m=ring(cK,K,K-A)&g.sector(cK,a_km,a_t)
    m|=ring(cM,M,M-A)&g.sector(cM,a_ml,a_km)
    m|=ring(cL,L,L-A)&g.sector(cL,a_lm,180)
    m|=(g.Y>=0)&(g.Y<=A)&(g.X<=D)&g.disk(cL,L)
    info={'cK':cK.round(1),'cM':cM.round(1),'cL':cL.round(1),'terminal':(round(xt,1),round(yt,1))}
    return g.bitmap(m),info
def three():
    g=Grid(); D,M,N,Q,G,I,S,F=DIM['D'],DIM['M'],DIM['N'],DIM['Q'],DIM['G'],DIM['I'],DIM['S'],DIM['F']
    xv=D-DIM['O']
    xp=DIM['P']                                   # the tongue's left end: an implied square cap at x = P (the I and N centres)
    yt=TOP-N; RT=((S-xv)**2+(TOP-yt)**2)/(2*(TOP-yt)); cT=np.array([xv,TOP-RT])    # top arc through (S, top-N), top on the overshoot line
    cG=np.array([D-S-G, cT[1]+np.sqrt((RT-G)**2-(D-S-G-xv)**2)])                   # upper bowl outer: right side S in from D, tangent to the top arc
    cI=np.array([D-S-A-I, cT[1]+np.sqrt((RT-A-I)**2-(D-S-A-I-xv)**2)])              # upper bowl inner (radius I), tangent to the top inner arc
    cQ=np.array([xv,BOT+Q])                                                         # lower bowl bottom on the overshoot line
    cM=np.array([D-M, cQ[1]-np.sqrt((Q-M)**2-(D-M-xv)**2)])                         # lower bowl outer right (radius M), tangent to Q
    cN=np.array([D-A-N, cQ[1]-np.sqrt((Q-A-N)**2-(D-A-N-xv)**2)])                   # lower bowl inner right (radius N), tangent to Q-A
    ang=lambda c,p: np.degrees(np.arctan2(p[1]-c[1],p[0]-c[0]))
    a_t=ang(cT,(S,yt)); a_tg=ang(cT,cG); a_ti=ang(cT,cI)
    yF=BOT+F; xF=cQ[0]-np.sqrt(Q*Q-(yF-cQ[1])**2); a_lt=ang(cQ,(xF,yF)); a_qm=ang(cQ,cM); a_qn=ang(cQ,cN)   # F measured from the overshoot line, like N
    tongue_top=cI[1]-I; tongue_bot=cN[1]+N
    outer=(g.disk(cT,RT)&g.sector(cT,a_tg,a_t))|(g.disk(cG,G)&g.sector(cG,-150,ang(cG,cG+(cG-cT))))
    outer|=(g.disk(cQ,Q)&g.sector(cQ,a_lt,a_qm))|(g.disk(cM,M)&((g.X>=xv)|(g.Y<=cM[1])))
    outer|=(g.X>=xp)&(g.X<=cI[0]+I)&(g.Y>=tongue_bot)&(g.Y<=tongue_top)
    # upper counter: the top inner arc, then the inner bowl arc I; the tongue's left side is a plain vertical at the construction line
    counter=(g.disk(cT,RT-A)&g.sector(cT,a_ti,a_t)&((g.X<xp)|(g.Y>cI[1])))|g.disk(cI,I)
    counter|=(g.disk(cQ,Q-A)&g.sector(cQ,a_lt,a_qn))|g.disk(cN,N)
    cap=(g.X>=xv)&(g.X<=xp)&(g.Y>=tongue_bot)&(g.Y<=tongue_top)   # rectangular extension of the bar to the construction line
    m=(outer&~counter)|cap
    m&=~((g.X<xv)&(g.Y>tongue_bot)&(g.Y<cG[1]))          # the middle opening: nothing left of the bar's end
    m&=~((g.X<xp)&~cap&(g.Y>tongue_bot)&(g.Y<cG[1]))
    m&=~((g.X<xp)&(g.X>=xv-40)&(g.Y>cN[1])&(g.Y<=tongue_bot)&~g.disk(cQ,Q-A)==False) if False else m
    m&=~((g.X<xp)&(g.X>=xv-40)&(g.Y>cN[1])&(g.Y<=tongue_bot+1))   # the lower counter reaches the cap cleanly
    info={'top_radius':round(RT,1),'cT':cT.round(1),'cG':cG.round(1),'cI':cI.round(1),'cQ':cQ.round(1),'cM':cM.round(1),'cN':cN.round(1),
          'tongue':(round(tongue_bot,1),round(tongue_top,1)),'lower terminal x':round(xF,1)}
    return g.bitmap(m),info
