"""Series D, page D-2 (G H I J K L)."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
R8={'A':'1 1/4','B':'5 3/8','C':'8','D':'2 9/16','E':'1 3/4','F':'1/8','G':'2 5/8','H':'3 1/2','I':'1 3/16','J':'1 7/16',
    'K':'5','L':'2 3/4','M':'3 1/8','N':'3/16','O':'5 1/2','P':'4 7/8'}
R12={'A':'1 7/8','B':'8 1/16','C':'12','D':'3 27/32','E':'2 5/8','F':'3/16','G':'3 15/16','H':'5 1/4','I':'1 25/32',
     'J':'2 5/32','K':'7 1/2','L':'4 1/8','M':'4 11/16','N':'9/32','O':'8 1/4','P':'7 5/16'}
R18={'A':'2 13/16','B':'12 3/32','C':'18','D':'5 3/4','E':'3 15/16','F':'9/32','G':'5 29/32','H':'7 7/8','I':'2 11/16',
     'J':'3 1/4','K':'11 1/4','L':'6 3/16','M':'7 1/32','N':'7/16','O':'12 3/8','P':'10 31/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['F']; TOP=CAP+OV; BOT=-OV
def letter_G():
    """big arc C (the letter height) centred at mid-height, C from the left; arcs D tangent to it and to the overshoot
    lines; upper terminal E below the top overshoot line, cut along the radius; bar top at mid-height, G long"""
    g=Grid(); B,C,D,E,G=[DIM[k] for k in 'BCDEG']
    cC=np.array([C,CAP/2]); yu=TOP-D; yl=BOT+D
    xc=cC[0]-np.sqrt((C-D)**2-(yu-cC[1])**2); cu=np.array([xc,yu]); cl=np.array([xc,yl])
    yt=TOP-E; xt=cu[0]+np.sqrt(D*D-(yt-cu[1])**2)
    xr=cl[0]+D                                # right side: the lower arc's right edge carried up (as Series E); B is nominal
    m=ring(g,cC,C,C-A)&g.sector(cC,ang(cC,cu),ang(cC,cl))
    m|=ring(g,cu,D,D-A)&g.sector(cu,ang(cu,(xt,yt)),ang(cC,cu))
    m|=ring(g,cl,D,D-A)&g.sector(cl,ang(cC,cl),0)&(g.X<=xr)
    m|=rect(g,xr-A,xr,cl[1],CAP/2)                                     # jaw: vertical down into the lower arc
    m|=rect(g,xr-G,xr,CAP/2-A,CAP/2)                                   # bar
    return m,{'upper tip x':round(xt,1),'lower arc right x':round(cl[0]+D,1),'B':round(B,1)}
def letter_H():
    g=Grid(); B,H=DIM['B'],DIM['H']
    return rect(g,0,A,0,CAP)|rect(g,B-A,B,0,CAP)|rect(g,0,B,H,H+A)
def letter_I():
    g=Grid(); return rect(g,0,A,0,CAP)
def letter_J():
    """stem A wide at the right; bowl G centred G in from the right and G above the bottom overshoot line;
    terminal J above the overshoot line, cut along the radius"""
    g=Grid(); G,J,K=DIM['G'],DIM['J'],DIM['K']
    c=np.array([K-G,BOT+G]); yt=BOT+J; xt=c[0]-np.sqrt(G*G-(yt-c[1])**2)
    m=rect(g,K-A,K,c[1],CAP)|(ring(g,c,G,G-A)&g.sector(c,ang(c,(xt,yt)),0))
    return m,{'tip x':round(xt,1)}
def letter_K():
    """arm's lower edge: stem at height L to the top-right corner (O-N, cap); crotch M below the cap on it;
    leg's right edge from the crotch to (O, 0); strokes A wide"""
    from dimE3 import strip
    g=Grid(); L,M,N,O=DIM['L'],DIM['M'],DIM['N'],DIM['O']
    P1=np.array([A,L]); P2=np.array([O-N,CAP]); u=(P2-P1)/np.linalg.norm(P2-P1); nu=np.array([-u[1],u[0]])
    yc=CAP-M; crotch=P1+u*((yc-P1[1])/u[1])
    v=crotch-np.array([O,0.0]); v/=np.linalg.norm(v); nv=np.array([-v[1],v[0]])
    if nv[0]>0: nv=-nv
    arm=strip(g,P1,u,A,P1+nu)&(g.X>=A-1)&(g.Y<=CAP)
    leg=strip(g,np.array([O,0.0]),v,A,np.array([O,0.0])+nv)&(g.Y>=0)&(g.Y<=yc)
    return rect(g,0,A,0,CAP)|arm|leg,{'crotch':crotch.round(1)}
def letter_L():
    g=Grid(); I,P=DIM['I'],DIM['P']
    return rect(g,0,A,0,CAP)|rect(g,0,P,0,I)
