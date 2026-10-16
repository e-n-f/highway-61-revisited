"""Series D, page D-3 (M N O P Q R)."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
R8={'A':'1 1/4','B':'6 3/16','C':'4 5/8','D':'2 11/16','E':'5 3/8','F':'1 1/8','G':'2 9/16','H':'8','I':'1/8','J':'3 1/2',
    'K':'2 1/4','L':'5 5/8','M':'2 1/4','N':'5/8','O':'3/8','P':'1 3/8'}
R12={'A':'1 7/8','B':'9 9/32','C':'6 15/16','D':'4 1/32','E':'8 1/16','F':'1 11/16','G':'3 27/32','H':'12','I':'3/16',
     'J':'5 1/4','K':'3 3/8','L':'8 7/16','M':'3 3/8','N':'15/16','O':'9/16','P':'2 1/16'}
R18={'A':'2 13/16','B':'13 15/16','C':'10 13/32','D':'6 1/16','E':'12 3/32','F':'2 17/32','G':'5 3/4','H':'18','I':'9/32',
     'J':'7 7/8','K':'5 1/16','L':'12 21/32','M':'5 1/16','N':'1 13/32','O':'27/32','P':'3 3/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['I']; TOP=CAP+OV; BOT=-OV
def letter_M():
    """diagonals meet at (B/2, 0); the right one's lower edge runs to the right stem at height C"""
    g=Grid(); B,C=DIM['B'],DIM['C']
    stems=(rect(g,0,A,0,CAP)|rect(g,B-A,B,0,CAP))
    v=np.array([B/2,0.0]); pr=np.array([B-A,C]); u=(pr-v)/np.linalg.norm(pr-v)
    right=strip(g,v,u,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X<=B-A+1)&(g.X>=B/2)
    ul=np.array([-u[0],u[1]]); left=strip(g,v,ul,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=B/2)
    return stems|right|left
def letter_N():
    """diagonal's lower edge from the left stem at D below the cap straight to (E-F, 0), trimming the right stem's foot"""
    g=Grid(); D,E,F=DIM['D'],DIM['E'],DIM['F']
    stems=rect(g,0,A,0,CAP)|rect(g,E-A,E,0,CAP)
    p0=np.array([A,CAP-D]); p1=np.array([E-F,0.0]); u=(p1-p0)/np.linalg.norm(p1-p0)
    diag=strip(g,p0,u,A,(E,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=E-A+1)
    below=((g.X-p0[0])*u[1]-(g.Y-p0[1])*u[0])>=0
    stems=stems&~(below&(g.X>=E/2))
    return stems|diag
def _O(g):
    """side arcs H centred at mid-height on the far side (width L); top and bottom arcs touch the overshoot lines.
    The table's G misses tangency with H by about 2 units, so G is refitted as a connecting arc (as in Series E)."""
    L,H=DIM['L'],DIM['H']
    cl=np.array([H,CAP/2]); cr=np.array([L-H,CAP/2])
    dx=H-L/2; h=TOP-CAP/2; G=(H*H-dx*dx-h*h)/(2*(H-h))
    cu=np.array([L/2,TOP-G]); cd=np.array([L/2,BOT+G])
    def area(rs,rt):
        lens=g.disk(cl,rs)&g.disk(cr,rs)
        top=g.sector(cu,ang(cr,cu),ang(cl,cu))&~g.disk(cu,rt)&(g.Y>cu[1])
        bot=g.sector(cd,ang(cl,cd),ang(cr,cd))&~g.disk(cd,rt)&(g.Y<cd[1])
        return lens&~top&~bot
    return area(H,G)&~area(H-A,G-A),{'G table':round(DIM['G'],1),'G tangent':round(G,1)}
def letter_O():
    g=Grid(); return _O(g)
def letter_Q():
    """tail: rectangle F wide; its left edge from (L/2-N, BOT+M) to the bottom corner (L/2+F, BOT-O)"""
    g=Grid(); m,info=_O(g); L,F,M,N,O=[DIM[k] for k in 'LFMNO']
    TL=np.array([L/2-N,BOT+M]); BL=np.array([L/2+F,BOT-O]); u=(BL-TL)/np.linalg.norm(BL-TL)
    n=np.array([u[1],-u[0]])
    if n[0]<0: n=-n                                  # to the right of the edge
    tail=strip(g,TL,u,F,TL+n)&strip(g,TL,n,np.linalg.norm(BL-TL),BL)
    return m|tail,info
def _bowl(g,E,K,J,I=0.0):
    c=np.array([E-I-K,(CAP+J)/2])
    outer=rect(g,0,c[0],J,CAP)|(g.disk(c,K)&(g.X>=c[0]))
    inner=rect(g,A,c[0],J+A,CAP-A)|(g.disk(c,K-A)&(g.X>=c[0]))
    return outer&~inner,c
def letter_P():
    g=Grid(); b,c=_bowl(g,DIM['E'],DIM['K'],DIM['J'])
    return b|rect(g,0,A,0,CAP),{'bowl half-height vs K':(round((CAP-DIM['J'])/2,1),round(DIM['K'],1))}
def letter_R():
    """bowl I in from the right; leg's right edge from (E-P, J) to the corner (E, 0), A wide"""
    g=Grid(); E,I,J,K,P=[DIM[k] for k in 'EIJKP']
    b,c=_bowl(g,E,K,J,I)
    p=np.array([E-P,J]); q=np.array([E,0.0]); u=(p-q)/np.linalg.norm(p-q)
    n=np.array([-u[1],u[0]])
    if n[0]>0: n=-n
    leg=strip(g,q,u,A,q+n)&(g.Y>=0)&(g.Y<=J+A)&~(g.disk(c,K-A)&(g.X>=c[0]))
    return b|leg|rect(g,0,A,0,CAP)
