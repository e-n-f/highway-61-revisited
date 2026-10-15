"""Series D, page D-6 (figures 5 6 7 8 9 0)."""
import numpy as np
from shapely.geometry import Polygon
from shapely import affinity
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
from dimE4 import arc
from dimE6 import oval
R8={'A':'1 1/4','B':'5 3/8','C':'5 3/4','D':'2 11/16','E':'4 1/16','F':'1/8','G':'4 1/8','H':'3/8','I':'7/16','J':'3/16',
    'K':'7 1/8','L':'2 7/8','M':'2 5/16','N':'2 7/16','O':'4 7/8','P':'12','Q':'2 3/8','R':'2 7/16','S':'2 1/16',
    'T':'2 3/4','U':'4 1/2','V':'15/16','W':'5/16','X':'5 5/8','Y':'8','Z':'2 9/16'}
R12={'A':'1 7/8','B':'8 1/16','C':'8 5/8','D':'4 1/32','E':'6 3/32','F':'3/16','G':'6 3/16','H':'9/16','I':'21/32',
     'J':'9/32','K':'10 11/16','L':'4 5/16','M':'3 15/32','N':'3 21/32','O':'7 5/16','P':'18','Q':'3 9/16','R':'3 21/32',
     'S':'3 3/32','T':'4 1/8','U':'6 3/4','V':'1 13/32','W':'15/32','X':'8 7/16','Y':'12','Z':'3 27/32'}
R18={'A':'2 13/16','B':'12 3/32','C':'12 15/16','D':'6 1/16','E':'9 5/32','F':'9/32','G':'9 9/32','H':'27/32','I':'1',
     'J':'7/16','K':'16 1/32','L':'6 15/32','M':'5 3/16','N':'5 15/32','O':'10 31/32','P':'27','Q':'5 11/32','R':'5 1/2',
     'S':'4 5/8','T':'6 3/16','U':'10 1/8','V':'2 1/8','W':'11/16','X':'12 21/32','Y':'18','Z':'5 3/4'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['F']; TOP=CAP+OV; BOT=-OV
deg=lambda v: np.degrees(np.arctan2(v[1],v[0]))
def pm(pts): return Region(Polygon([tuple(map(float,p)) for p in pts]).buffer(0))
def fig5():
    """bowl D touching the right side, its top C above the bottom overshoot line; bottom-left arc E touching the overshoot
    line, internally tangent to D; stem's left edge from (J, cap-G) to (2J, cap), A wide; its foot runs to the notch on
    the bowl's inner arc, H lower; top bar to I short of B; lower terminal tip on the left edge, cut along the radius
    (the table's I above the overshoot line is about 4 units higher: over-determined, as in Series E)"""
    g=Grid(); B,C,D,E,G,H,I,J=[DIM[k] for k in 'BCDEGHIJ']
    cD=np.array([B-D,BOT+C-D])
    yE=BOT+E; cE=np.array([cD[0]-np.sqrt((E-D)**2-(yE-cD[1])**2),yE])
    yn=CAP-G-H; xn=cD[0]-np.sqrt((D-A)**2-(yn-cD[1])**2)        # notch on the bowl's inner arc
    p0=np.array([J,CAP-G]); p1=np.array([2*J,CAP]); u=(p1-p0)/np.linalg.norm(p1-p0)
    foot=lambda: ((g.Y-p0[1])*(xn-p0[0])-(g.X-p0[0])*(yn-p0[1]))>=0     # above the line from the stem's foot to the notch
    xt=0.0; yt=cE[1]-np.sqrt(E*E-(xt-cE[0])**2)
    a_DE=ang(cE,cD)                                                   # tangent point direction (from E through D)
    m=ring(g,cD,D,D-A)&g.sector(cD,a_DE,180)&(foot()|(g.X>=xn))
    m|=ring(g,cE,E,E-A)&g.sector(cE,ang(cE,(xt,yt)),a_DE)
    stem=strip(g,p0,u,A,(B,CAP))&(g.Y<=CAP)&foot()
    m|=stem
    m&=~(strip(g,p0,u,1e4,(-1e3,CAP))&~strip(g,p0,u,0.0,(-1e3,CAP))&(g.Y>=yn))  # nothing left of the stem's left edge
    m|=rect(g,2*J,B-I,CAP-A,CAP)
    return m,{'lower tip y':round(yt,1),'BOT+I':round(BOT+I,1),'notch':(round(xn,1),round(yn,1))}
def bowl6():
    """four-centre oval B wide, O tall from the bottom overshoot line: top/bottom arcs L, side arcs touching x=0 and B.
    The side radius M misses tangency with L by 0.5 units; refitted as the connecting arc (195.9 for 198.9)."""
    B,L,O=DIM['B'],DIM['L'],DIM['O']
    yc=BOT+O/2; h=L-O/2; M=(L*L-B*B/4-h*h)/(2*L-B)
    ct=np.array([B/2,BOT+O-L]); cb=np.array([B/2,BOT+L]); cl=np.array([M,yc]); cr=np.array([B-M,yc])
    return pm(oval(ct,L,cl,M,cb,L,cr,M)), pm(oval(ct,L-A,cl,M-A,cb,L-A,cr,M-A)), yc, M
def fig6():
    """tail: arc K centred on the bowl's centre line at x=K, so it leaves the bowl's left side tangentially; inner edge
    K-A; ends at the cap line, cut along the radius"""
    g=Grid(); K=DIM['K']
    out,inn,yc,M=bowl6()
    c=np.array([K,yc]); xt=c[0]-np.sqrt(K*K-(CAP-c[1])**2)
    tail=ring(g,c,K,K-A)&g.sector(c,ang(c,(xt,CAP)),180)
    return (out|tail)&~inn,{'side radius':round(M,1),'M':round(DIM['M'],1),'N vs O/2':(round(DIM['N'],1),round(DIM['O']/2,1))}
def fig9():
    m,info=fig6()
    return Region(affinity.rotate(m.geom,180,origin=(DIM['B']/2,CAP/2))),info
def fig7():
    """stem: ring A wide whose right edge (radius P-A) runs from (B, cap-A) to (B-Q, 0); centre to the right"""
    g=Grid(); B,P,Q=DIM['B'],DIM['P'],DIM['Q']
    r=P-A; p=np.array([B-Q,0.0]); q=np.array([B,CAP-A]); mid=(p+q)/2; ch=q-p; Lc=np.linalg.norm(ch)
    perp=np.array([ch[1],-ch[0]])/Lc; c=mid+perp*np.sqrt(r*r-(Lc/2)**2)
    if c[0]<mid[0]: c=mid-perp*np.sqrt(r*r-(Lc/2)**2)
    return (ring(g,c,P,r)&(g.Y>=0)&(g.Y<=CAP-A)&(g.X<=B))|rect(g,0,B,CAP-A,CAP),{'centre':c.round(1)}
def _ovals():
    B,R,S,T,W=[DIM[k] for k in 'BRSTW']
    xc=B/2
    # upper bowl: top arc R on the top overshoot line; side arcs S touching x=W and B-W, centred above R's centre
    ct=np.array([xc,TOP-R]); dx=xc-W-S; ys=ct[1]+np.sqrt((R-S)**2-dx*dx)
    up=lambda d: oval(ct,R-d,np.array([W+S,ys]),S-d,np.array([xc,2*ys-ct[1]]),R-d,np.array([B-W-S,ys]),S-d)
    # lower bowl: bottom arc T on the bottom overshoot line; side arcs R touching 0 and B, centred below T's centre
    cb=np.array([xc,BOT+T]); dx2=xc-R; ys2=cb[1]-np.sqrt((T-R)**2-dx2*dx2)
    lo=lambda d: oval(np.array([xc,2*ys2-cb[1]]),T-d,np.array([R,ys2]),R-d,cb,T-d,np.array([B-R,ys2]),R-d)
    return up,lo
def fig8():
    up,lo=_ovals()
    m=(pm(up(0))|pm(lo(0)))&~pm(up(A))&~pm(lo(A))
    from shapely.geometry import LineString
    U,V,B=DIM['U'],DIM['V'],DIM['B']
    bd=Polygon(up(0)).exterior.intersection(Polygon(lo(0)).exterior)
    pts=sorted([(p.x,p.y) for p in getattr(bd,'geoms',[bd])],key=lambda p:-p[0])
    return m,{'right notch':tuple(round(v,1) for v in pts[0]),'table (B-V, BOT+U)':(round(B-V,1),round(BOT+U,1))}
def fig0():
    """same dimensions as the letter O (width X, side arcs Y, top/bottom Z): built by the same construction"""
    import dimD3
    return dimD3.letter_O()
