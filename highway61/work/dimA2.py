"""Series A, page A-2 (G H I J K L).  The constructions of Series B (dimB2) with this page's letters."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
R8={'A':'3/4','B':'3','C':'1 3/8','D':'3 7/8','E':'18','F':'1/8','G':'1 3/4','H':'3 3/4','I':'1 9/16','J':'2 7/8',
    'K':'3 1/4','L':'2 11/16','M':'11/16','N':'1'}
R12={'A':'1 1/8','B':'4 1/2','C':'2 1/16','D':'5 13/16','E':'27','F':'3/16','G':'2 5/8','H':'5 5/8','I':'2 11/32',
     'J':'4 5/16','K':'4 7/8','L':'4 1/32','M':'1 1/32','N':'1 1/2'}
R18={'A':'1 11/16','B':'6 3/4','C':'3 3/32','D':'8 23/32','E':'40 1/2','F':'9/32','G':'3 15/16','H':'8 7/16',
     'I':'3 33/64','J':'6 15/32','K':'7 5/16','L':'6 1/16','M':'1 9/16','N':'2 1/4'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['F']; TOP=CAP+OV; BOT=-OV
def letter_G():
    """big arc E centred at mid-height (left side at x=0); arcs C tangent to it and to the overshoot lines; the upper
    one runs straight on from its rightmost point to G below the overshoot line; bar C long with its top D above the
    bottom overshoot line; the jaw drops from the bar into the lower arc at its right extreme"""
    g=Grid(); B,C,D,E,G=[DIM[k] for k in 'BCDEG']
    cE=np.array([E,CAP/2]); yu=TOP-C; yl=BOT+C
    xc=cE[0]-np.sqrt((E-C)**2-(yu-cE[1])**2); cu=np.array([xc,yu]); cl=np.array([xc,yl])
    xr=xc+C; yb=BOT+D
    m=ring(g,cE,E,E-A)&g.sector(cE,ang(cE,cu),ang(cE,cl))
    m|=ring(g,cu,C,C-A)&g.sector(cu,0,ang(cE,cu))
    m|=ring(g,cl,C,C-A)&g.sector(cl,ang(cE,cl),0)
    m|=rect(g,xr-A,xr,TOP-G,yu)
    m|=rect(g,xr-A,xr,cl[1],yb)|rect(g,xr-C,xr,yb-A,yb)
    return m,{'right x':round(xr,1),'B':round(B,1)}
def letter_H():
    """stems B apart outside; bar A thick with its lower edge H above the baseline"""
    g=Grid(); B,H=DIM['B'],DIM['H']
    return rect(g,0,A,0,CAP)|rect(g,B-A,B,0,CAP)|rect(g,0,B,H,H+A)
def letter_I():
    g=Grid(); return rect(g,0,A,0,CAP)
def letter_J():
    """stem at the right of width B; bowl I centred I in from the right and I above the bottom overshoot line;
    terminal N above the overshoot line, cut along the radius"""
    g=Grid(); B,I,N=DIM['B'],DIM['I'],DIM['N']
    c=np.array([B-I,BOT+I]); yt=BOT+N; xt=c[0]-np.sqrt(I*I-(yt-c[1])**2)
    return rect(g,B-A,B,c[1],CAP)|(ring(g,c,I,I-A)&g.sector(c,ang(c,(xt,yt)),0)),{'tip x':round(xt,1)}
def letter_K():
    """arm's lower edge from the stem at height K to the top-right corner (J, cap); crotch B below the cap on it;
    leg's right edge from the crotch to (B, 0)"""
    g=Grid(); B,J,K=DIM['B'],DIM['J'],DIM['K']
    P1=np.array([A,K]); P2=np.array([J,CAP]); u=(P2-P1)/np.linalg.norm(P2-P1); nu=np.array([-u[1],u[0]])
    yc=CAP-B; crotch=P1+u*((yc-P1[1])/u[1])
    v=crotch-np.array([B,0.0]); v/=np.linalg.norm(v); nv=np.array([-v[1],v[0]])
    if nv[0]>0: nv=-nv
    arm=strip(g,P1,u,A,P1+nu)&(g.X>=A-1)&(g.Y<=CAP)
    leg=strip(g,np.array([B,0.0]),v,A,np.array([B,0.0])+nv)&(g.Y>=0)&(g.Y<=yc)
    return rect(g,0,A,0,CAP)|arm|leg,{'crotch':crotch.round(1)}
def letter_L():
    """foot L long and M thick"""
    g=Grid(); L,M=DIM['L'],DIM['M']
    return rect(g,0,A,0,CAP)|rect(g,0,L,0,M)
