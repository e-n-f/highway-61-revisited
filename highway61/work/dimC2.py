"""Series C, page C-2 (G H I J K L)."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
R8={'A':'1 1/8','B':'4 3/8','C':'2 1/32','D':'10','E':'2','F':'4','G':'1/8','H':'3 5/8','I':'2 1/8','J':'1',
    'K':'1 1/16','L':'4 3/16','M':'2 9/16','N':'2 13/16'}
R12={'A':'1 11/16','B':'6 9/16','C':'3 3/64','D':'15','E':'3','F':'6','G':'3/16','H':'5 7/16','I':'3 3/16','J':'1 1/2',
     'K':'1 19/32','L':'6 9/32','M':'3 27/32','N':'4 7/32'}
R18={'A':'2 17/32','B':'9 27/32','C':'4 9/16','D':'22 1/2','E':'4 1/2','F':'9','G':'9/32','H':'8 5/32','I':'4 25/32',
     'J':'2 1/4','K':'2 3/8','L':'9 7/16','M':'5 3/4','N':'6 5/16'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['G']; TOP=CAP+OV; BOT=-OV
def letter_G():
    """big arc D centred at mid-height (left side at x=0); arcs C tangent to it and to the overshoot lines; the upper one
    ends at its rightmost point with a horizontal cut (as the C); bar E long with its top F above the bottom overshoot line; the jaw
    drops from the bar into the lower arc at its right extreme"""
    g=Grid(); B,C,D,E,F=[DIM[k] for k in 'BCDEF']
    cD=np.array([D,CAP/2]); yu=TOP-C; yl=BOT+C
    xc=cD[0]-np.sqrt((D-C)**2-(yu-cD[1])**2); cu=np.array([xc,yu]); cl=np.array([xc,yl])
    xr=xc+C
    m=ring(g,cD,D,D-A)&g.sector(cD,ang(cD,cu),ang(cD,cl))
    m|=ring(g,cu,C,C-A)&g.sector(cu,0,ang(cD,cu))
    m|=ring(g,cl,C,C-A)&g.sector(cl,ang(cD,cl),0)
    yb=BOT+F                                          # bar top: F above the bottom overshoot line (as the scan shows)
    m|=rect(g,xr-A,xr,cl[1],yb)|rect(g,xr-E,xr,yb-A,yb)
    return m,{'right x':round(xr,1),'B':round(B,1)}
def letter_H():
    g=Grid(); B,H=DIM['B'],DIM['H']
    return rect(g,0,A,0,CAP)|rect(g,B-A,B,0,CAP)|rect(g,0,B,H,H+A)
def letter_I():
    g=Grid(); return rect(g,0,A,0,CAP)
def letter_J():
    """stem at the right of width F; bowl I centred I in from the right and I above the bottom overshoot line; terminal
    J above the overshoot line, cut along the radius"""
    g=Grid(); F,I,J=DIM['F'],DIM['I'],DIM['J']
    c=np.array([F-I,BOT+I]); yt=BOT+J; xt=c[0]-np.sqrt(I*I-(yt-c[1])**2)
    return rect(g,F-A,F,c[1],CAP)|(ring(g,c,I,I-A)&g.sector(c,ang(c,(xt,yt)),0)),{'tip x':round(xt,1)}
def letter_K():
    """arm's lower edge from the stem at height M to the top-right corner (L, cap); crotch N below the cap on it;
    leg's right edge from the crotch to (B, 0)"""
    g=Grid(); B,L,M,N=DIM['B'],DIM['L'],DIM['M'],DIM['N']
    P1=np.array([A,M]); P2=np.array([L,CAP]); u=(P2-P1)/np.linalg.norm(P2-P1); nu=np.array([-u[1],u[0]])
    yc=CAP-N; crotch=P1+u*((yc-P1[1])/u[1])
    v=crotch-np.array([B,0.0]); v/=np.linalg.norm(v); nv=np.array([-v[1],v[0]])
    if nv[0]>0: nv=-nv
    arm=strip(g,P1,u,A,P1+nu)&(g.X>=A-1)&(g.Y<=CAP)
    leg=strip(g,np.array([B,0.0]),v,A,np.array([B,0.0])+nv)&(g.Y>=0)&(g.Y<=yc)
    return rect(g,0,A,0,CAP)|arm|leg,{'crotch':crotch.round(1)}
def letter_L():
    g=Grid(); F,K=DIM['F'],DIM['K']
    return rect(g,0,A,0,CAP)|rect(g,0,F,0,K)
