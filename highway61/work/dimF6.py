"""Series F, page F-6 (figures 5 6 7 8 9 0)."""
import numpy as np
from shapely import affinity
from shapely.geometry import Polygon
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
from dimE4 import arc
from dimE6 import oval
from dimD6 import pm
R8={'A':'1 1/2','B':'7 3/8','C':'8 1/4','D':'5 3/4','E':'6','F':'2 5/8','G':'3/4','H':'1/4','I':'1/8','J':'4 3/4',
    'K':'5','L':'2 1/4','M':'1 3/8','N':'4','O':'3','P':'4 1/2','Q':'3/8','R':'2 3/8','S':'5 5/8','T':'3 3/4',
    'U':'9 1/8','V':'2 7/8','W':'2','X':'1','Y':'3 11/16','Z':'7 5/8'}
R12={'A':'2 1/4','B':'11 1/16','C':'12 3/8','D':'8 5/8','E':'9','F':'3 15/16','G':'1 1/8','H':'3/8','I':'3/16',
     'J':'7 1/8','K':'7 1/2','L':'3 3/8','M':'2 1/16','N':'6','O':'4 1/2','P':'6 3/4','Q':'9/16','R':'3 9/16',
     'S':'8 7/16','T':'5 5/8','U':'13 11/16','V':'4 5/16','W':'3','X':'1 1/2','Y':'5 17/32','Z':'11 7/16'}
R18={'A':'3 3/8','B':'16 19/32','C':'18 9/16','D':'12 15/16','E':'13 1/2','F':'5 29/32','G':'1 11/16','H':'9/16',
     'I':'9/32','J':'10 11/16','K':'11 1/4','L':'5 1/16','M':'3 3/32','N':'9','O':'6 3/4','P':'10 1/8','Q':'27/32',
     'R':'5 11/32','S':'12 21/32','T':'8 7/16','U':'20 17/32','V':'6 15/32','W':'4 1/2','X':'2 1/4','Y':'8 5/16',
     'Z':'17 5/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['I']; TOP=CAP+OV; BOT=-OV
def fig5():
    """bowl: top arc E centred on the axis with its top D above the bottom overshoot line; right arc F touching the right
    side, internally tangent to E; bottom arc C on the overshoot line, internally tangent to F.  Stem's left edge from
    (H, cap-N) to (2H, cap), A wide, its foot running to the notch on E-A, H lower; top bar to G short of B; lower
    terminal tip on the left edge, cut along the radius"""
    g=Grid(); B,C,D,E,F,G,H,N=[DIM[k] for k in 'BCDEFGHN']
    cE=np.array([B/2,BOT+D-E])
    cF=np.array([B-F,cE[1]+np.sqrt((E-F)**2-(B-F-B/2)**2)])
    yC=BOT+C; cC=np.array([cF[0]-np.sqrt((C-F)**2-(yC-cF[1])**2),yC])
    yn=CAP-N-H; xn=cE[0]-np.sqrt((E-A)**2-(yn-cE[1])**2)
    p0=np.array([H,CAP-N]); p1=np.array([2*H,CAP]); u=(p1-p0)/np.linalg.norm(p1-p0)
    foot=lambda: ((g.Y-p0[1])*(xn-p0[0])-(g.X-p0[0])*(yn-p0[1]))>=0
    xt=0.0; yt=cC[1]-np.sqrt(C*C-(xt-cC[0])**2)
    a_EF,a_FC=ang(cE,cF),ang(cF,cC)
    m=ring(g,cE,E,E-A)&g.sector(cE,a_EF,180)&(foot()|(g.X>=xn))&(g.Y>=cE[1])
    m|=ring(g,cF,F,F-A)&g.sector(cF,ang(cC,cF),a_EF)
    m|=ring(g,cC,C,C-A)&g.sector(cC,ang(cC,(xt,yt)),ang(cC,cF))
    m|=strip(g,p0,u,A,(B,CAP))&(g.Y<=CAP)&foot()
    m&=~(strip(g,p0,u,1e4,(-1e3,CAP))&~strip(g,p0,u,0.0,(-1e3,CAP))&(g.Y>=yn))
    m|=rect(g,2*H,B-G,CAP-A,CAP)
    return m,{'cF':cF.round(1),'cC':cC.round(1),'tip y':round(yt,1),'notch':(round(xn,1),round(yn,1))}
TAIL='J'
K_FROM_BASELINE=True
def bowl6():
    """four-centre oval B wide, K tall from the bottom overshoot line: top/bottom arcs D, side arcs refitted to be
    tangent to them while touching x=0 and B on the centre line (table L)"""
    B,D,K=DIM['B'],DIM['D'],DIM['K']
    Kh=K-BOT if K_FROM_BASELINE else K                  # bowl height: from the bottom overshoot line up to K
    yc=BOT+Kh/2; h=D-Kh/2; K=Kh
    # side centre (L, yc), top centre (B/2, BOT+K-D): (B/2-L)^2 + h^2 = (D-L)^2
    L=(D*D-(B/2)**2-h*h)/(2*D-B)
    ct=np.array([B/2,BOT+K-D]); cb=np.array([B/2,BOT+D]); cl=np.array([L,yc]); cr=np.array([B-L,yc])
    return pm(oval(ct,D,cl,L,cb,D,cr,L)),pm(oval(ct,D-A,cl,L-A,cb,D-A,cr,L-A)),cl,L
def fig6():
    """tail: an arc (radius J) internally tangent to the bowl's left side arc and reaching the cap line at B-M; inner
    edge concentric; tip cut along the radius"""
    from scipy.optimize import fsolve
    g=Grid(); B,M=DIM['B'],DIM['M']; r=DIM[TAIL]
    out,inn,cl,L=bowl6()
    tip=np.array([B-M,CAP])
    c=np.array(fsolve(lambda v:[np.hypot(*(v-cl))-(r-L),np.hypot(*(v-tip))-r],[cl[0]+r-L,cl[1]+50]))
    tail=ring(g,c,r,r-A)&g.sector(c,ang(c,tip),ang(c,cl))        # from the tip round to the tangent point
    return (out|tail)&~inn,{'side L refit':round(L,1),'L table':round(DIM['L'],1),'tail centre':c.round(1)}
def fig9():
    m,info=fig6()
    return Region(affinity.rotate(m.geom,180,origin=(DIM['B']/2,CAP/2))),info
def fig7():
    """stem: ring A wide whose right edge (radius C) runs from (B-O, 0) to the top bar's corner (B, cap-A)"""
    g=Grid(); B,C,O=DIM['B'],DIM['C'],DIM['O']
    p=np.array([B-O,0.0]); q=np.array([B,CAP-A]); mid=(p+q)/2; ch=q-p; Lc=np.linalg.norm(ch)
    perp=np.array([ch[1],-ch[0]])/Lc; c=mid+perp*np.sqrt(C*C-(Lc/2)**2)
    if c[0]<mid[0]: c=mid-perp*np.sqrt(C*C-(Lc/2)**2)
    return (ring(g,c,C+A,C)&(g.Y>=0)&(g.Y<=CAP-A)&(g.X<=B))|rect(g,0,B,CAP-A,CAP),{'centre':c.round(1)}
def _oval_tangent_y(ca,ra,xs,rs,above):
    """height of a side arc's centre (at x=xs) so that the side arc rs is internally tangent to the end arc (ca, ra)"""
    dy=np.sqrt((ra-rs)**2-(xs-ca[0])**2); return ca[1]+dy if above else ca[1]-dy
def fig8():
    """two four-centre ovals.  Upper: top arc S on the overshoot line and side arcs W touching Q in from the sides
    (centred where tangent to S); the counter's bottom is E (inner), tangent to the counter sides W-A, and the bowl's
    bottom its offset E+A.  Lower: bottom arc U on the overshoot line and side arcs R touching the sides (tangent to
    U); the counter's top is T (inner), tangent to R-A, the bowl's top T+A.  Counters concentric (A).  The notch is
    where the bowls cross."""
    B,E,P,Q,R,S,T,U,W,X=[DIM[k] for k in 'BEPQRSTUWX']
    xc=B/2
    cS=np.array([xc,TOP-S]); yw=_oval_tangent_y(cS,S,Q+W,W,True)
    cE=np.array([xc,yw+np.sqrt((E-W+A)**2-(Q+W-xc)**2)])
    cU=np.array([xc,BOT+U]); yr=_oval_tangent_y(cU,U,R,R,False)
    cT=np.array([xc,yr-np.sqrt((T-R+A)**2-(R-xc)**2)])
    wl,wr=np.array([Q+W,yw]),np.array([B-Q-W,yw]); rl,rr=np.array([R,yr]),np.array([B-R,yr])
    up_o=pm(oval(cS,S,wl,W,cE,E+A,wr,W)); up_i=pm(oval(cS,S-A,wl,W-A,cE,E,wr,W-A))
    lo_o=pm(oval(cT,T+A,rl,R,cU,U,rr,R)); lo_i=pm(oval(cT,T,rl,R-A,cU,U-A,rr,R-A))
    m=(up_o|lo_o)&~up_i&~lo_i
    bd=up_o.geom.exterior.intersection(lo_o.geom.exterior)
    pts=sorted([(p.x,p.y) for p in getattr(bd,'geoms',[bd])],key=lambda p:-p[0])
    return m,{'notch':tuple(round(v,1) for v in pts[0]) if pts else None,'table (B-X, BOT+P)':(round(B-X,1),round(BOT+P,1)),
              'cE':cE.round(1),'cT':cT.round(1)}
def fig0():
    """four-centre oval Z wide: side arcs J centred on the far side at mid-height (touching x=0 and Z), top/bottom arcs
    on the overshoot lines refitted for tangency (table Y)"""
    g=Grid(); Z,J=DIM['Z'],DIM['J']
    cl=np.array([J,CAP/2]); cr=np.array([Z-J,CAP/2])
    dx=J-Z/2; h=TOP-CAP/2; Y=(J*J-dx*dx-h*h)/(2*(J-h))
    ct=np.array([Z/2,TOP-Y]); cb=np.array([Z/2,BOT+Y])
    out=pm(oval(ct,Y,cl,J,cb,Y,cr,J)); inn=pm(oval(ct,Y-A,cl,J-A,cb,Y-A,cr,J-A))
    return out&~inn,{'Y refit':round(Y,1),'Y table':round(DIM['Y'],1)}
