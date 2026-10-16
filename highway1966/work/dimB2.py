"""Series B, page B-2 (G H I J K L).  Same constructions as Series C (dimC2)."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
R8={'A':'1','B':'3 3/8','C':'4','D':'1 9/16','E':'16','F':'1/8','G':'1 1/2','H':'3 3/4','I':'1 3/4','J':'3 1/2',
    'K':'3','L':'2 3/4','M':'15/16'}
R12={'A':'1 1/2','B':'5 1/16','C':'6','D':'2 11/32','E':'24','F':'3/16','G':'2 1/4','H':'5 5/8','I':'2 5/8','J':'5 1/4',
     'K':'4 1/2','L':'4 1/8','M':'1 13/32'}
R18={'A':'2 1/4','B':'7 19/32','C':'9','D':'3 1/2','E':'36','F':'9/32','G':'3 3/8','H':'8 7/16','I':'3 15/16',
     'J':'7 7/8','K':'6 3/4','L':'6 3/16','M':'2 1/8'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['F']; TOP=CAP+OV; BOT=-OV
BAR_FROM_BOT=True
def letter_G():
    """big arc E centred at mid-height (left side at x=0); arcs D tangent to it and to the overshoot lines; the upper one
    ends at its rightmost point with a horizontal cut; bar G long with its top C above the bottom overshoot line; the
    jaw drops from the bar into the lower arc at its right extreme"""
    g=Grid(); B,C,D,E,G=[DIM[k] for k in 'BCDEG']
    cE=np.array([E,CAP/2]); yu=TOP-D; yl=BOT+D
    xc=cE[0]-np.sqrt((E-D)**2-(yu-cE[1])**2); cu=np.array([xc,yu]); cl=np.array([xc,yl])
    xr=xc+D; yb=(BOT if BAR_FROM_BOT else 0)+C
    m=ring(g,cE,E,E-A)&g.sector(cE,ang(cE,cu),ang(cE,cl))
    m|=ring(g,cu,D,D-A)&g.sector(cu,0,ang(cE,cu))
    m|=ring(g,cl,D,D-A)&g.sector(cl,ang(cE,cl),0)
    m|=rect(g,xr-A,xr,cl[1],yb)|rect(g,xr-G,xr,yb-A,yb)
    return m,{'right x':round(xr,1),'B':round(B,1)}
def letter_H():
    g=Grid(); B,H=DIM['B'],DIM['H']
    return rect(g,0,A,0,CAP)|rect(g,B-A,B,0,CAP)|rect(g,0,B,H,H+A)
def letter_I():
    g=Grid(); return rect(g,0,A,0,CAP)
def letter_J():
    """stem at the right of width J; bowl I centred I in from the right and I above the bottom overshoot line; terminal
    G above the overshoot line, cut along the radius"""
    g=Grid(); G,I,J=DIM['G'],DIM['I'],DIM['J']
    c=np.array([J-I,BOT+I]); yt=BOT+G; xt=c[0]-np.sqrt(I*I-(yt-c[1])**2)
    return rect(g,J-A,J,c[1],CAP)|(ring(g,c,I,I-A)&g.sector(c,ang(c,(xt,yt)),0)),{'tip x':round(xt,1)}
def letter_K():
    """arm's lower edge from the stem at height K to the top-right corner (J-F, cap); crotch L below the cap on it;
    leg's right edge from the crotch to (J, 0)"""
    g=Grid(); F,J,K,L=DIM['F'],DIM['J'],DIM['K'],DIM['L']
    P1=np.array([A,K]); P2=np.array([J-F,CAP]); u=(P2-P1)/np.linalg.norm(P2-P1); nu=np.array([-u[1],u[0]])
    yc=CAP-L; crotch=P1+u*((yc-P1[1])/u[1])
    v=crotch-np.array([J,0.0]); v/=np.linalg.norm(v); nv=np.array([-v[1],v[0]])
    if nv[0]>0: nv=-nv
    arm=strip(g,P1,u,A,P1+nu)&(g.X>=A-1)&(g.Y<=CAP)
    leg=strip(g,np.array([J,0.0]),v,A,np.array([J,0.0])+nv)&(g.Y>=0)&(g.Y<=yc)
    return rect(g,0,A,0,CAP)|arm|leg,{'crotch':crotch.round(1)}
def letter_L():
    g=Grid(); K,M=DIM['K'],DIM['M']
    return rect(g,0,A,0,CAP)|rect(g,0,K,0,M)
