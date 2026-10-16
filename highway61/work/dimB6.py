"""Series B, page B-6 (figures 5 6 7 8 9 0).  The scan's 18-inch row is partly illegible: dimensions from the 12-inch
row, checked against the 8-inch row."""
import numpy as np
from shapely import affinity
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
from dimE6 import oval
from dimD6 import pm
R8={'A':'1','B':'3 3/8','C':'1/8','D':'1/4','E':'3 1/4','F':'3 7/8','G':'1 5/8','H':'2 3/16','I':'4 3/8','J':'5 1/2',
    'K':'4 7/8','L':'7 3/8','M':'1 9/16','N':'2 7/16','O':'17 1/2','P':'7/8','Q':'4','R':'4 11/16','S':'1 3/8',
    'T':'1 1/2','U':'5/8','V':'16','W':'2 5/16','X':'3 5/8'}
R12={'A':'1 1/2','B':'5 1/16','C':'3/16','D':'3/8','E':'4 7/8','F':'5 13/16','G':'2 7/16','H':'3 9/32','I':'6 9/16',
     'J':'8 1/4','K':'7 5/16','L':'11 1/16','M':'2 11/32','N':'3 21/32','O':'26 1/4','P':'1 5/16','Q':'6','R':'7 1/32',
     'S':'2 1/16','T':'2 1/4','U':'15/16','V':'24','W':'3 15/32','X':'5 7/16'}
DIM=units(R12,12); CHECK=check({8:R8,12:R12})
A=DIM['A']; CAP=690.0; OV=DIM['C']; TOP=CAP+OV; BOT=-OV
def fig5():
    """bowl: top arc G (top J above the bottom overshoot line), right arc I centred E above that line and touching the
    right side, bottom arc H on the overshoot line; G and H internally tangent to I.  Stem's left edge from (C, cap-F) to
    (2C, cap), A wide, its foot running to the notch on G-A, C lower; top bar to D short of B; lower terminal tip on the
    left edge, cut along the radius"""
    g=Grid(); B,C,D,E,F,G,H,I,J=[DIM[k] for k in 'BCDEFGHIJ']
    cI=np.array([B-I,BOT+E])
    yG=BOT+J-G; cG=np.array([cI[0]+np.sqrt((I-G)**2-(yG-cI[1])**2),yG])
    yH=BOT+H; cH=np.array([cI[0]+np.sqrt((I-H)**2-(yH-cI[1])**2),yH])
    yn=CAP-F-C; xn=cG[0]-np.sqrt((G-A)**2-(yn-cG[1])**2)
    p0=np.array([C,CAP-F]); p1=np.array([2*C,CAP]); u=(p1-p0)/np.linalg.norm(p1-p0)
    foot=lambda: ((g.Y-p0[1])*(xn-p0[0])-(g.X-p0[0])*(yn-p0[1]))>=0
    xt=0.0; yt=cH[1]-np.sqrt(H*H-(xt-cH[0])**2)
    a_IG,a_IH=ang(cI,cG),ang(cI,cH)
    m=ring(g,cG,G,G-A)&g.sector(cG,a_IG,180)&(foot()|(g.X>=xn))
    m|=ring(g,cI,I,I-A)&g.sector(cI,a_IH,a_IG)
    m|=ring(g,cH,H,H-A)&g.sector(cH,ang(cH,(xt,yt)),a_IH)
    m|=strip(g,p0,u,A,(B,CAP))&(g.Y<=CAP)&foot()
    m&=~(strip(g,p0,u,1e4,(-1e3,CAP))&~strip(g,p0,u,0.0,(-1e3,CAP))&(g.Y>=yn))
    m|=rect(g,2*C,B-D,CAP-A,CAP)
    return m,{'G right x':round(cG[0]+G,1),'B':round(B,1),'tip y':round(yt,1),'notch':(round(xn,1),round(yn,1))}
def bowl6():
    """four-centre oval B wide, K tall from the bottom overshoot line (centre line N): side arcs L touching x=0 and B;
    top/bottom arcs refitted for tangency with L (table G)"""
    B,K,L=DIM['B'],DIM['K'],DIM['L']
    yc=BOT+K/2; G=(L*L-(L-B/2)**2-K*K/4)/(2*L-K)
    ct=np.array([B/2,BOT+K-G]); cb=np.array([B/2,BOT+G]); cl=np.array([L,yc]); cr=np.array([B-L,yc])
    return pm(oval(ct,G,cl,L,cb,G,cr,L)),pm(oval(ct,G-A,cl,L-A,cb,G-A,cr,L-A)),yc,G
def fig6():
    """tail: the bowl's left side arc L carried on up to the cap line (inner edge L-A), cut along the radius"""
    g=Grid(); L=DIM['L']
    out,inn,yc,G=bowl6()
    c=np.array([L,yc]); xt=c[0]-np.sqrt(L*L-(CAP-c[1])**2)
    tail=ring(g,c,L,L-A)&g.sector(c,ang(c,(xt,CAP)),180)
    return (out|tail)&~inn,{'G refit':round(G,1),'G table':round(DIM['G'],1),'N vs K/2':(round(DIM['N'],1),round(DIM['K']/2,1))}
def fig9():
    m,info=fig6()
    return Region(affinity.rotate(m.geom,180,origin=(DIM['B']/2,CAP/2))),info
def fig7():
    """stem: ring centred on the baseline whose left edge (radius O) passes through (P, 0); its right edge meets the top
    bar's underside at the right side"""
    g=Grid(); B,O,P=DIM['B'],DIM['O'],DIM['P']
    c=np.array([P+O,0.0])
    return (ring(g,c,O,O-A)&(g.Y>=0)&(g.Y<=CAP-A)&(g.X<=B))|rect(g,0,B,CAP-A,CAP),\
        {'right edge at cap-A':round(c[0]-np.sqrt((O-A)**2-(CAP-A)**2),1),'B':round(B,1)}
def fig8():
    """two four-centre ovals.  Upper: top arc S on the overshoot line, side arcs Q on the line W below it, touching C in
    from the sides, bottom arc T tangent to them.  Lower: bottom arc M on the overshoot line, side arcs R on the line N
    above it touching x=0 and B, top arc T tangent to them.  Counters concentric (A).  Notch: where the bowls cross."""
    B,C,I,M,N,Q,R,S,T,U,W=[DIM[k] for k in 'BCIMNQRSTUW']
    xc=B/2; yW=TOP-W; yN=BOT+N
    ctu=np.array([xc,TOP-S]); clu=np.array([C+Q,yW]); cru=np.array([B-C-Q,yW])
    cbu=np.array([xc,yW-np.sqrt((Q-T)**2-(C+Q-xc)**2)])
    cbl=np.array([xc,BOT+M]); cll=np.array([R,yN]); crl=np.array([B-R,yN])
    ctl=np.array([xc,yN+np.sqrt((R-T)**2-(R-xc)**2)])
    up=lambda d: pm(oval(ctu,S-d,clu,Q-d,cbu,T-d,cru,Q-d))
    lo=lambda d: pm(oval(ctl,T-d,cll,R-d,cbl,M-d,crl,R-d))
    m=(up(0)|lo(0))&~up(A)&~lo(A)
    bd=up(0).geom.exterior.intersection(lo(0).geom.exterior)
    pts=sorted([(p.x,p.y) for p in getattr(bd,'geoms',[bd])],key=lambda p:-p[0])
    return m,{'notch':tuple(round(v,1) for v in pts[0]),'table (B-U, BOT+I)':(round(B-U,1),round(BOT+I,1)),
              'upper tangency S/Q gap':round((Q-S)-np.linalg.norm(clu-ctu),1),'lower M/R gap':round((R-M)-np.linalg.norm(cll-cbl),1)}
def fig0():
    """same dimensions as the letter O (width X, side arcs V, top/bottom M): the same construction"""
    import dimB3
    return dimB3.letter_O()
