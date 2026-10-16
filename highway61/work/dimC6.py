"""Series C, page C-6 (figures 5 6 7 8 9 0)."""
import numpy as np
from shapely import affinity
from vgeom import VGrid as Grid, Region, fuse as vgeom_fuse
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
from dimE4 import arc
from dimE6 import oval
from dimD6 import pm
R8={'A':'1 1/8','B':'4 3/8','C':'5 1/2','D':'2 1/8','E':'2 1/2','F':'4 1/8','G':'3 1/8','H':'5/16','I':'3/16','J':'1/8',
    'K':'11/16','L':'2 3/16','M':'7','N':'5','O':'12 3/4','P':'1','Q':'4 1/2','R':'3 1/32','S':'2 1/16','T':'15/16',
    'U':'3/4','V':'7/8','W':'2 9/16','X':'4 5/8','Y':'2 1/32','Z':'10'}
R12={'A':'1 11/16','B':'6 9/16','C':'8 1/4','D':'3 3/16','E':'3 3/4','F':'6 3/16','G':'4 11/16','H':'15/32','I':'9/32',
     'J':'3/16','K':'1 1/32','L':'3 9/32','M':'10 1/2','N':'7 1/2','O':'19 1/8','P':'1 1/2','Q':'6 3/4','R':'4 35/64',
     'S':'3 3/32','T':'1 13/32','U':'1 1/8','V':'1 5/16','W':'3 27/32','X':'6 15/16','Y':'3 3/64','Z':'15'}
R18={'A':'2 17/32','B':'9 27/32','C':'12 3/8','D':'4 25/32','E':'5 5/8','F':'9 9/32','G':'7 1/32','H':'11/16','I':'7/16',
     'J':'9/32','K':'1 9/16','L':'4 15/16','M':'15 3/4','N':'11 1/4','O':'28 11/16','P':'2 1/4','Q':'10 1/8','R':'6 13/16',
     'S':'4 41/64','T':'2 7/64','U':'1 11/16','V':'1 31/32','W':'5 49/64','X':'10 13/32','Y':'4 9/16','Z':'22 1/2'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
def fig5():
    """bowl: top arc D (top C above the bottom overshoot line), right arc F centred near the left edge at height G, both
    touching the right side; bottom-left arc E on the overshoot line; D and E internally tangent to F.  Stem's left edge
    from (I, cap-F) to (2I, cap), A wide, its foot running to the notch on D-A, J lower; top bar to H short of B; lower
    terminal K above the overshoot line, cut along the radius"""
    g=Grid(); B,C,D,E,F,G,H,I,J,K=[DIM[k] for k in 'BCDEFGHIJK']
    cF=np.array([B-F,G])
    yD=BOT+C-D; cD=np.array([cF[0]+np.sqrt((F-D)**2-(yD-cF[1])**2),yD])
    yE=BOT+E; cE=np.array([cF[0]+np.sqrt((F-E)**2-(yE-cF[1])**2),yE])
    yn=CAP-F-J; xn=cD[0]-np.sqrt((D-A)**2-(yn-cD[1])**2)
    p0=np.array([I,CAP-F]); p1=np.array([2*I,CAP]); u=(p1-p0)/np.linalg.norm(p1-p0)
    foot=lambda: ((g.Y-p0[1])*(xn-p0[0])-(g.X-p0[0])*(yn-p0[1]))>=0
    yt=BOT+K; xt=cE[0]-np.sqrt(E*E-(yt-cE[1])**2)
    a_FD,a_FE=ang(cF,cD),ang(cF,cE)
    m=ring(g,cD,D,D-A)&g.sector(cD,a_FD,180)&(foot()|(g.X>=xn))
    m|=ring(g,cF,F,F-A)&g.sector(cF,a_FE,a_FD)
    m|=ring(g,cE,E,E-A)&g.sector(cE,ang(cE,(xt,yt)),a_FE)
    m|=strip(g,p0,u,A,(B,CAP))&(g.Y<=CAP)&foot()
    m&=~(strip(g,p0,u,1e4,(-1e3,CAP))&~strip(g,p0,u,0.0,(-1e3,CAP))&(g.Y>=yn))
    m|=rect(g,2*I,B-H,CAP-A,CAP)
    return m,{'right extremes D/F':(round(cD[0]+D,1),round(cF[0]+F,1)),'tip x':round(xt,1),'notch':(round(xn,1),round(yn,1))}
def bowl6():
    """four-centre oval B wide, N tall from the bottom overshoot line, centre line E above it: side arcs M (centred on
    the far side, touching x=0 and B); top/bottom arcs L refitted for tangency with M (the table's L misses by 1.5)"""
    B,E,M,N=DIM['B'],DIM['E'],DIM['M'],DIM['N']
    yc=BOT+N/2; L=(M*M-(M-B/2)**2-N*N/4)/(2*M-N)
    ct=np.array([B/2,BOT+N-L]); cb=np.array([B/2,BOT+L]); cl=np.array([M,yc]); cr=np.array([B-M,yc])
    return pm(oval(ct,L,cl,M,cb,L,cr,M)),pm(oval(ct,L-A,cl,M-A,cb,L-A,cr,M-A)),yc,L
def fig6():
    """tail: the bowl's left side arc M carried on up to the cap line (inner edge M-A), cut along the radius"""
    g=Grid(); M=DIM['M']
    out,inn,yc,L=bowl6()
    c=np.array([M,yc]); xt=c[0]-np.sqrt(M*M-(CAP-c[1])**2)
    tail=ring(g,c,M,M-A)&g.sector(c,ang(c,(xt,CAP)),180)
    return (out|tail)&~inn,{'L refit':round(L,1),'L table':round(DIM['L'],1),'E vs N/2':(round(DIM['E'],1),round(DIM['N']/2,1))}
def fig9():
    m,info=fig6()
    return Region(affinity.rotate(m.geom,180,origin=(DIM['B']/2,CAP/2))),info
def fig7():
    """stem: ring centred on the baseline whose left edge (radius O) passes through (P, 0); its right edge (O-A) meets
    the top bar's underside at the right side"""
    g=Grid(); B,O,P=DIM['B'],DIM['O'],DIM['P']
    c=np.array([P+O,0.0])
    return (ring(g,c,O,O-A)&(g.Y>=0)&(g.Y<=CAP-A)&(g.X<=B))|rect(g,0,B,CAP-A,CAP),\
        {'right edge at cap-A':round(c[0]-np.sqrt((O-A)**2-(CAP-A)**2),1),'B':round(B,1)}
def fig8():
    """upper bowl: circle S on the top overshoot line (J in from the sides) with side arcs W below its centre line,
    touching it at the extremes; counter: V at the same centre, carried straight down to a V semicircle A above the lower
    counter (the waist).  Lower bowl: side arcs R on the centre line E
    (touching x=0 and B), bottom arc S on the overshoot line (refitted for tangency with R); counter: R-A sides, S-A
    bottom, T on top tangent to both sides.  The notch is where W and R cross."""
    g=Grid(); B,D,E,J,Q,R,S,T,U,V,W=[DIM[k] for k in 'BDEJQRSTUVW']
    xc=B/2; cu=np.array([xc,TOP-S])
    cwl=np.array([J+W,cu[1]]); cwr=np.array([B-J-W,cu[1]])
    up=(g.disk(cu,S)&(g.Y>=cu[1]))|(g.disk(cwl,W)&g.disk(cwr,W)&(g.Y<=cu[1]))
    yc=BOT+E; crl=np.array([R,yc]); crr=np.array([B-R,yc])
    Sb=(R*R-(R-xc)**2-E*E)/(2*(R-E)); cb=np.array([xc,BOT+Sb])         # bottom arc refitted: tangent to R
    a_l,a_r=ang(crl,cb),ang(crr,cb)                                  # bottom arc between its tangent points with R
    lens=g.disk(crl,R)&g.disk(crr,R); bot=g.sector(cb,a_l,a_r)
    lo=(lens&~bot)|(g.disk(cb,Sb)&bot)
    yT=yc+np.sqrt((R-A-T)**2-(R-xc)**2); cT=np.array([xc,yT])        # counter top: T tangent to both R-A sides
    lens_in=g.disk(crl,R-A)&g.disk(crr,R-A); tp=g.sector(cT,ang(crr,cT),ang(crl,cT))
    lo_cnt=(lens_in&~bot&~tp)|(g.disk(cT,T)&tp)|(g.disk(cb,Sb-A)&bot)
    lo_cnt=Region(vgeom_fuse(lo_cnt.geom).buffer(1e-4,join_style='mitre'))   # weld the pieces before cutting the counter
    # upper counter: circle V on the bowl's centre, carried straight down to a second V semicircle whose bottom is
    # A above the lower counter's top (the waist is A thick, as dimensioned); the scan shows this stadium shape
    cv=np.array([xc,yT+T+A+V]); up_cnt=g.disk(cu,V)|g.disk(cv,V)|rect(g,xc-V,xc+V,cv[1],cu[1])
    m=(up|lo)&~up_cnt&~lo_cnt
    # notch on the right: crossing of the right-hand side arcs W and R (both centred left of the axis)
    from shapely.geometry import Point
    cW=Point(*cwr).buffer(W,quad_segs=512).exterior; cR=Point(*crr).buffer(R,quad_segs=512).exterior
    pts=[p for p in getattr(cW.intersection(cR),'geoms',[]) ]
    notch=max(((p.x,p.y) for p in pts),key=lambda p:p[0]) if pts else None
    return m,{'notch':tuple(round(v,1) for v in notch) if notch else None,'table (B-U, BOT+Q)':(round(B-U,1),round(BOT+Q,1)),
              'bottom S refit':round(Sb,1),'upper counter straight sides':round(cu[1]-(yT+T+A+V),1),'upper centre vs TOP-D':(round(cu[1],1),round(TOP-D,1))}
def fig0():
    """same dimensions as the letter O (width X, side arcs Z, top/bottom Y): the same construction"""
    import dimC3
    return dimC3.letter_O()
