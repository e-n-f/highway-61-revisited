"""Series A, page A-6 (figures 5 6 7 8 9 0).  Constructions as in Series B (dimB6) with this page's letters."""
import numpy as np
from shapely import affinity
from vgeom import VGrid as Grid, Region
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
from dimE6 import oval
from dimD6 import pm
R8={'A':'3/4','B':'3','C':'1 3/8','D':'4 1/4','E':'1 15/16','F':'1/16','G':'1/8','H':'4','I':'5 3/8','J':'15',
    'K':'5 1/8','L':'4 7/8','M':'3 7/8','N':'1 7/16','O':'2 7/16','P':'36','Q':'1/4','R':'3 5/8','S':'6 3/4',
    'T':'2 1/8','U':'4 3/8','V':'1 5/16','W':'1 1/4','X':'9/16','Y':'1/2','Z':'18','AA':'3 1/4'}
R12={'A':'1 1/8','B':'4 1/2','C':'2 1/16','D':'6 3/8','E':'2 29/32','F':'3/32','G':'3/16','H':'6','I':'8 1/16',
     'J':'22 1/2','K':'7 11/16','L':'7 5/16','M':'5 13/16','N':'2 5/32','O':'3 21/32','P':'54','Q':'3/8',
     'R':'5 7/16','S':'10 1/8','T':'3 3/16','U':'6 9/16','V':'1 31/32','W':'1 7/8','X':'27/32','Y':'3/4','Z':'27',
     'AA':'4 7/8'}
R18={'A':'1 11/16','B':'6 3/4','C':'3 3/32','D':'9 9/16','E':'4 3/8','F':'1/8','G':'9/32','H':'9','I':'12 3/32',
     'J':'33 3/4','K':'11 17/32','L':'10 31/32','M':'8 23/32','N':'3 1/4','O':'5 1/2','P':'81','Q':'9/16',
     'R':'8 5/32','S':'15 3/16','T':'4 25/32','U':'9 27/32','V':'2 15/16','W':'2 13/16','X':'1 1/4','Y':'1 1/8',
     'Z':'40 1/2','AA':'7 5/16'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['G']; TOP=CAP+OV; BOT=-OV
I_FROM='BOT'; E_ON='overshoot'; D_FROM='BOT'      # heights on the curved figure from the overshoot line (best fit)
def fig5():
    """width B.  Bowl: right arc D centred B above the bottom overshoot line, touching the right side; top arc C (its
    top I above that line) and bottom arc E (on it), both internally tangent to D.  (Measured from the baseline, B and I
    fit the drawing worse: p95 8.2 against 3.3 units.)  Stem A wide, its
    left edge from (F, cap) to (0, cap-H), its foot at cap-H running to the counter's arc C-A; top bar A thick to
    G short of B.  Lower terminal: tip on the left edge, cut along E's radius."""
    g=Grid(); B,C,D,E,F,G,H,I=[DIM[k] for k in 'BCDEFGHI']
    cD=np.array([B-D,B if D_FROM=='baseline' else BOT+B])
    yC=(I if I_FROM=='baseline' else BOT+I)-C; cC=np.array([cD[0]+np.sqrt((D-C)**2-(yC-cD[1])**2),yC])
    yE=(BOT if E_ON=='overshoot' else 0.0)+E; cE=np.array([cD[0]+np.sqrt((D-E)**2-(yE-cD[1])**2),yE])
    yf=CAP-H; p0=np.array([0.0,yf]); p1=np.array([F,CAP]); u=(p1-p0)/np.linalg.norm(p1-p0)
    stem=strip(g,p0,u,A,(B,yf))&(g.Y>=yf)&(g.Y<=CAP)
    xt=0.0; yt=cE[1]-np.sqrt(E*E-(xt-cE[0])**2)
    a_DC,a_DE=ang(cD,cC),ang(cD,cE)
    m=ring(g,cC,C,C-A)&g.sector(cC,a_DC,270)&(g.Y>=yf)
    m|=ring(g,cD,D,D-A)&g.sector(cD,a_DE,a_DC)
    m|=ring(g,cE,E,E-A)&g.sector(cE,ang(cE,(xt,yt)),a_DE)
    m|=stem|rect(g,F,B-G,CAP-A,CAP)
    xn=cC[0]-np.sqrt(max((C-A)**2-(yf-cC[1])**2,0))
    return m,{'tip y':round(yt,1),'notch x on C-A at the foot':round(xn,1),'C left x':round(cC[0]-C,1)}
def bowl6():
    """four-centre oval B wide, L tall from the bottom overshoot line (centre line O): side arcs J touching x=0 and B;
    top/bottom arcs refitted for tangency with J (table N)"""
    B,L,J=DIM['B'],DIM['L'],DIM['J']
    yc=BOT+L/2; G=(J*J-(J-B/2)**2-L*L/4)/(2*J-L)
    ct=np.array([B/2,BOT+L-G]); cb=np.array([B/2,BOT+G]); cl=np.array([J,yc]); cr=np.array([B-J,yc])
    return pm(oval(ct,G,cl,J,cb,G,cr,J)),pm(oval(ct,G-A,cl,J-A,cb,G-A,cr,J-A)),yc,G
TIP=TOP
def fig6():
    """tail: arc K centred M above the baseline on the right, internally tangent to the bowl's left side arc J (inner
    edge K-A), reaching the top overshoot line, cut along the radius"""
    g=Grid(); B,J,K,M=DIM['B'],DIM['J'],DIM['K'],DIM['M']
    out,inn,yc,G=bowl6()
    cJ=np.array([J,yc]); cK=np.array([cJ[0]-np.sqrt((J-K)**2-(M-cJ[1])**2),M])
    yt=TIP; xt=cK[0]-np.sqrt(K*K-(yt-cK[1])**2)
    tail=ring(g,cK,K,K-A)&g.sector(cK,ang(cK,(xt,yt)),ang(cJ,cK)+30)    # past the tangency the tail lies in the bowl
    return (out|tail)&~inn,{'N refit':round(G,1),'N table':round(DIM['N'],1),'O vs L/2':(round(DIM['O'],1),round(DIM['L']/2,1)),
                            'tail tip x':round(xt,1),'K centre':cK.round(1)}
def fig9():
    m,info=fig6()
    return Region(affinity.rotate(m.geom,180,origin=(DIM['B']/2,CAP/2))),info
def fig7():
    """width B; top bar A thick; stem: a ring of radius P (outer, left edge) and P-A, centred far to the right, its left
    edge through (Q, 0) and its right edge reaching the top bar's underside at the right side (B, cap-A)"""
    g=Grid(); B,P,Q=DIM['B'],DIM['P'],DIM['Q']
    from scipy.optimize import fsolve
    c=np.array(fsolve(lambda v:[np.hypot(v[0]-Q,v[1])-P,np.hypot(v[0]-B,v[1]-(CAP-A))-(P-A)],[Q+P,0.0]))
    return (ring(g,c,P,P-A)&(g.Y>=0)&(g.Y<=CAP-A)&(g.X<=B))|rect(g,0,B,CAP-A,CAP),{'centre':c.round(1)}
def fig8():
    """two four-centre ovals.  Upper: top arc W on the overshoot line, side arcs R on the line T below it, touching F in
    from the sides, bottom arc W tangent to them.  Lower: bottom arc N on the overshoot line, side arcs S on the line T
    above it touching x=0 and B, top arc tangent to them.  Counters concentric (A; W-A is the table's Y).  The waist
    between the counters is dimensioned A thick: with the lower top arc W it comes out 36 units, so that arc is refitted
    (W -> about 118) to put the lower counter's top A below the upper one's bottom; this also brings the notch, where
    the bowls cross, from 8.5 above to 3.7 below the table's U above the bottom overshoot line."""
    from scipy.optimize import brentq
    B,F,N,R,S,T,U,W,Y=[DIM[k] for k in 'BFNRSTUWY']
    xc=B/2; yT=TOP-T; yB=BOT+T
    ctu=np.array([xc,TOP-W]); clu=np.array([F+R,yT]); cru=np.array([B-F-R,yT])
    cbu=np.array([xc,yT-np.sqrt((R-W)**2-(F+R-xc)**2)])
    cbl=np.array([xc,BOT+N]); cll=np.array([S,yB]); crl=np.array([B-S,yB])
    yct=lambda w: yB+np.sqrt((S-w)**2-(S-xc)**2)             # centre of a top arc w tangent to the sides S
    W2=brentq(lambda w: yct(w)+w-A-(cbu[1]-(W-A)-A),W,xc-1)    # lower counter top A below the upper counter bottom
    ctl=np.array([xc,yct(W2)])
    up=lambda d: pm(oval(ctu,W-d,clu,R-d,cbu,W-d,cru,R-d))
    lo=lambda d: pm(oval(ctl,W2-d,cll,S-d,cbl,N-d,crl,S-d))
    m=(up(0)|lo(0))&~up(A)&~lo(A)
    bd=up(0).geom.exterior.intersection(lo(0).geom.exterior)
    pts=sorted([(p.x,p.y) for p in getattr(bd,'geoms',[bd])],key=lambda p:-p[0])
    return m,{'notch':tuple(round(v,1) for v in pts[0]),'table notch height':round(BOT+U,1),'lower top arc':round(W2,1),'W-A vs Y':(round(W-A,1),round(Y,1)),
              'upper top W/R gap':round((R-W)-np.linalg.norm(clu-ctu),1),'lower N/S gap':round((S-N)-np.linalg.norm(cll-cbl),1)}
def fig0():
    """width AA, side arcs Z: the dimensions of the letter O (G, I on page A-3), so the same construction; the top and
    bottom arcs refit for tangency (table V)"""
    import dimA3
    return dimA3.letter_O()
