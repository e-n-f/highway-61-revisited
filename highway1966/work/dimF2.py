"""Series F, page F-2 (G H I J K L).  The scan cuts off the 18-inch row: dimensions from the 12-inch row, checked
against the 8-inch row."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
R8={'A':'1 1/2','B':'7 3/8','C':'4 3/4','D':'3 7/8','E':'3 5/8','F':'3 3/8','G':'1 3/4','H':'1/8','I':'3 1/2',
    'J':'3 1/4','K':'2 5/8','L':'6 1/2','M':'6 3/4','N':'5/8','O':'2 3/4','P':'2 1/4','Q':'1 3/8','R':'4 1/2'}
R12={'A':'2 1/4','B':'11 1/16','C':'7 1/8','D':'5 13/16','E':'5 7/16','F':'5 1/16','G':'2 5/8','H':'3/16','I':'5 1/4',
     'J':'4 7/8','K':'3 15/16','L':'9 3/4','M':'10 1/8','N':'15/16','O':'4 1/8','P':'3 3/8','Q':'2 1/16','R':'6 3/4'}
DIM=units(R12,12); CHECK=check({8:R8,12:R12})
A=DIM['A']; CAP=690.0; OV=DIM['H']; TOP=CAP+OV; BOT=-OV
def letter_G():
    """left arc C centred at mid-height C from the left; upper arc D and lower arc E on the overshoot lines, internally
    tangent to C; upper terminal G below the top overshoot line, cut along the radius; bar F long, top R above the
    bottom overshoot line; jaw at the lower arc's right extreme"""
    g=Grid(); B,C,D,E,F,G,R=[DIM[k] for k in 'BCDEFGR']
    cC=np.array([C,CAP/2])
    yD=TOP-D; cD=np.array([cC[0]-np.sqrt((C-D)**2-(yD-cC[1])**2),yD])
    yE=BOT+E; cE=np.array([cC[0]-np.sqrt((C-E)**2-(yE-cC[1])**2),yE])
    yt=TOP-G; xt=cD[0]+np.sqrt(D*D-(yt-cD[1])**2)
    xr=cE[0]+E; yb=BOT+R
    m=ring(g,cC,C,C-A)&g.sector(cC,ang(cC,cD),ang(cC,cE))
    m|=ring(g,cD,D,D-A)&g.sector(cD,ang(cD,(xt,yt)),ang(cC,cD))
    m|=ring(g,cE,E,E-A)&g.sector(cE,ang(cC,cE),0)
    # jaw and bar; below the lower arc's centre the right edge follows that arc (the bar's foot lies below it)
    m|=(rect(g,xr-A,xr,cE[1],yb)|rect(g,xr-F,xr,yb-A,yb))&((g.Y>=cE[1])|g.disk(cE,E))
    return m,{'right x':round(xr,1),'B':round(B,1),'upper tip x':round(xt,1)}
def letter_H():
    g=Grid(); B,I=DIM['B'],DIM['I']
    return rect(g,0,A,0,CAP)|rect(g,B-A,B,0,CAP)|rect(g,0,B,I,I+A)
def letter_I():
    g=Grid(); return rect(g,0,A,0,CAP)
def letter_J():
    """stem at the right of width L; bowl J (= L/2) centred J above the bottom overshoot line; terminal K above the
    overshoot line, cut along the radius"""
    g=Grid(); J,K,L=DIM['J'],DIM['K'],DIM['L']
    c=np.array([L-J,BOT+J]); yt=BOT+K; xt=c[0]-np.sqrt(J*J-(yt-c[1])**2)
    return rect(g,L-A,L,c[1],CAP)|(ring(g,c,J,J-A)&g.sector(c,ang(c,(xt,yt)),0)),{'tip x':round(xt,1)}
def letter_K():
    """arm's lower edge from the stem at height P to the top-right corner (B-N, cap); crotch O below the cap on it;
    leg's right edge from the crotch to (B, 0)"""
    g=Grid(); B,N,O,P=DIM['B'],DIM['N'],DIM['O'],DIM['P']
    P1=np.array([A,P]); P2=np.array([B-N,CAP]); u=(P2-P1)/np.linalg.norm(P2-P1); nu=np.array([-u[1],u[0]])
    yc=CAP-O; crotch=P1+u*((yc-P1[1])/u[1])
    v=crotch-np.array([B,0.0]); v/=np.linalg.norm(v); nv=np.array([-v[1],v[0]])
    if nv[0]>0: nv=-nv
    arm=strip(g,P1,u,A,P1+nu)&(g.X>=A-1)&(g.Y<=CAP)
    leg=strip(g,np.array([B,0.0]),v,A,np.array([B,0.0])+nv)&(g.Y>=0)&(g.Y<=yc)
    return rect(g,0,A,0,CAP)|arm|leg,{'crotch':crotch.round(1)}
def letter_L():
    g=Grid(); M,Q=DIM['M'],DIM['Q']
    return rect(g,0,A,0,CAP)|rect(g,0,M,0,Q)
