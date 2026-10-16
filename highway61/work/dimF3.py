"""Series F, page F-3 (M N O P Q R).  The scan cuts off the 18-inch row: 12-inch row, checked against the 8-inch."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip, cross_solve
R8={'A':'1 1/2','B':'8 1/2','C':'4 3/4','D':'7 3/8','E':'2 1/2','F':'1 3/8','G':'8 1/4','H':'4 1/8','I':'3','J':'1/8',
    'K':'2 5/8','L':'2 1/8','M':'3/8'}
R12={'A':'2 1/4','B':'12 3/4','C':'7 1/8','D':'11 1/16','E':'3 3/4','F':'2 1/16','G':'12 3/8','H':'6 3/16','I':'4 1/2',
     'J':'3/16','K':'3 15/16','L':'3 3/16','M':'9/16'}
DIM=units(R12,12); CHECK=check({8:R8,12:R12})
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
def letter_M():
    g=Grid(); B,C=DIM['B'],DIM['C']
    stems=rect(g,0,A,0,CAP)|rect(g,B-A,B,0,CAP)
    v=np.array([B/2,0.0]); pr=np.array([B-A,C]); u=(pr-v)/np.linalg.norm(pr-v)
    right=strip(g,v,u,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X<=B-A+1)&(g.X>=B/2)
    ul=np.array([-u[0],u[1]]); left=strip(g,v,ul,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=B/2)
    return stems|right|left
def letter_N():
    g=Grid(); D,E,F=DIM['D'],DIM['E'],DIM['F']
    stems=rect(g,0,A,0,CAP)|rect(g,D-A,D,0,CAP)
    p0=np.array([A,CAP-E]); p1=np.array([D-F,0.0]); u=(p1-p0)/np.linalg.norm(p1-p0)
    diag=strip(g,p0,u,A,(D,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=D-A+1)
    below=((g.X-p0[0])*u[1]-(g.Y-p0[1])*u[0])>=0
    return (stems&~(below&(g.X>=D/2)))|diag
def letter_O():
    """a circle: width G = the overshoot-to-overshoot height, radius H"""
    g=Grid(); G,H=DIM['G'],DIM['H']
    return ring(g,np.array([G/2,CAP/2]),H,H-A),{'H vs G/2':(round(H,1),round(G/2,1))}
def letter_Q():
    """tail: rectangle F wide; its left edge from (G/2+J, BOT+K) to the bottom corner (G/2+L, BOT-M)"""
    g=Grid(); m,_=letter_O(); G,F,J,K,L,M=[DIM[k] for k in 'GFJKLM']
    TL=np.array([G/2+J,BOT+K]); BL=np.array([G/2+L,BOT-M]); u=(BL-TL)/np.linalg.norm(BL-TL)
    n=np.array([u[1],-u[0]])
    if n[0]<0: n=-n
    tail=strip(g,TL,u,F,TL+n)&strip(g,TL,n,np.linalg.norm(BL-TL),BL)
    return m|tail
def _bowl(g,right,bot):
    """semicircular bowl E (half its height) touching x=right; counter concentric"""
    E=DIM['E']; c=np.array([right-E,(CAP+bot)/2])
    outer=rect(g,0,c[0],bot,CAP)|(g.disk(c,E)&(g.X>=c[0]))
    inner=rect(g,A,c[0],bot+A,CAP-A)|(g.disk(c,E-A)&(g.X>=c[0]))
    return outer&~inner,c
def letter_P():
    g=Grid(); b,c=_bowl(g,DIM['D'],DIM['I'])
    return b|rect(g,0,A,0,CAP),{'E vs half bowl':(round(DIM['E'],1),round((CAP-DIM['I'])/2,1))}
def letter_R():
    """bowl as P, J in from the right; leg: left edge from (A+E, I), right edge to the corner (D, 0), A apart"""
    g=Grid(); D,E,I,J=[DIM[k] for k in 'DEIJ']
    b,c=_bowl(g,D-J,I)
    p=np.array([A+E,I]); q=np.array([D,0.0]); dx,dy=q-p
    th=cross_solve(p,q,A,-np.pi/2+0.01,np.arctan2(dy,dx)); u=np.array([np.cos(th),np.sin(th)])
    leg=strip(g,p,u,A,q)&(g.Y>=0)&(g.Y<=I+A)
    return b|leg|rect(g,0,A,0,CAP),{'leg angle':round(np.degrees(th),1)}
