"""Series B, page B-3 (M N O P Q R)."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
from dimC1 import _counter_end
R8={'A':'1','B':'3 7/8','C':'3 3/8','D':'3 1/4','E':'7/8','F':'3 5/8','G':'1 9/16','H':'16','I':'1/8','J':'3 1/2',
    'K':'2 1/8','L':'1 3/4','M':'3/4','N':'4 3/4','O':'3/8','P':'5/16','Q':'3 3/4','R':'13/16'}
R12={'A':'1 1/2','B':'5 13/16','C':'5 1/16','D':'4 7/8','E':'1 5/16','F':'5 7/16','G':'2 11/32','H':'24','I':'3/16',
     'J':'5 1/4','K':'3 3/16','L':'2 5/8','M':'1 1/8','N':'7 1/8','O':'9/16','P':'15/32','Q':'5 5/8','R':'1 7/32'}
R18={'A':'2 1/4','B':'8 23/32','C':'7 19/32','D':'7 5/16','E':'1 31/32','F':'8 5/32','G':'3 1/2','H':'36','I':'9/32',
     'J':'7 7/8','K':'4 25/32','L':'3 15/16','M':'1 11/16','N':'10 11/16','O':'27/32','P':'11/16','Q':'8 7/16',
     'R':'1 13/16'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['I']; TOP=CAP+OV; BOT=-OV
def letter_M():
    """diagonals' lower edges meet at (B/2, L); the right one runs to the right stem at height N"""
    g=Grid(); B,L,N=DIM['B'],DIM['L'],DIM['N']
    stems=rect(g,0,A,0,CAP)|rect(g,B-A,B,0,CAP)
    v=np.array([B/2,L]); pr=np.array([B-A,N]); u=(pr-v)/np.linalg.norm(pr-v)
    right=strip(g,v,u,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X<=B-A+1)&(g.X>=B/2)
    ul=np.array([-u[0],u[1]]); left=strip(g,v,ul,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=B/2)
    return stems|right|left
def letter_N():
    g=Grid(); C,D,E=DIM['C'],DIM['D'],DIM['E']
    stems=rect(g,0,A,0,CAP)|rect(g,C-A,C,0,CAP)
    p0=np.array([A,CAP-D]); p1=np.array([C-E,0.0]); u=(p1-p0)/np.linalg.norm(p1-p0)
    diag=strip(g,p0,u,A,(C,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=C-A+1)
    below=((g.X-p0[0])*u[1]-(g.Y-p0[1])*u[0])>=0
    return (stems&~(below&(g.X>=C/2)))|diag
def _O(g):
    """side arcs H centred at mid-height on the far side (width F); top/bottom arcs refitted for tangency (table G)"""
    W,R=DIM['F'],DIM['H']
    cl=np.array([R,CAP/2]); cr=np.array([W-R,CAP/2])
    dx=R-W/2; h=TOP-CAP/2; G=(R*R-dx*dx-h*h)/(2*(R-h))
    cu=np.array([W/2,TOP-G]); cd=np.array([W/2,BOT+G])
    def area(rs,rt):
        lens=g.disk(cl,rs)&g.disk(cr,rs)
        top=g.sector(cu,ang(cr,cu),ang(cl,cu))&~g.disk(cu,rt)&(g.Y>cu[1])
        bot=g.sector(cd,ang(cl,cd),ang(cr,cd))&~g.disk(cd,rt)&(g.Y<cd[1])
        return lens&~top&~bot
    return area(R,G)&~area(R-A,G-A),{'G table':round(DIM['G'],1),'G tangent':round(G,1)}
def letter_O():
    g=Grid(); return _O(g)
def letter_Q():
    """tail: rectangle M wide; its left edge from (F/2-P, BOT+G) to the bottom corner (F/2+A, BOT-O)"""
    g=Grid(); m,info=_O(g); W,G,M,O,P=[DIM[k] for k in 'FGMOP']
    TL=np.array([W/2-P,BOT+G]); BL=np.array([W/2+A,BOT-O]); u=(BL-TL)/np.linalg.norm(BL-TL)
    n=np.array([u[1],-u[0]])
    if n[0]<0: n=-n
    tail=strip(g,TL,u,M,TL+n)&strip(g,TL,n,np.linalg.norm(BL-TL),BL)
    return m|tail,info
def _bowl(g,top,bot):
    """flat-topped bowl C wide: side arc C (centred on the stem's left edge) on the centre line K below the cap,
    corners L; counter = the bowl offset by A"""
    C,K,L=DIM['C'],DIM['K'],DIM['L']
    c=np.array([C-C,CAP-K])
    out,_,_=_counter_end(g,0,top,bot,c,C,L)
    cnt,_,_=_counter_end(g,A,top-A,bot+A,c,C-A,L-A)
    return out,cnt,c
def letter_P():
    g=Grid(); out,cnt,c=_bowl(g,CAP,DIM['J'])
    return (out|rect(g,0,A,0,CAP))&~cnt
def letter_R():
    """bowl as P, bottom at Q; leg's right edge from the corner (C, 0) up to (C-R, Q)"""
    g=Grid(); C,Q,R=DIM['C'],DIM['Q'],DIM['R']
    out,cnt,c=_bowl(g,CAP,Q)
    p=np.array([C-R,Q]); q=np.array([C,0.0]); u=(p-q)/np.linalg.norm(p-q)
    n=np.array([-u[1],u[0]])
    if n[0]>0: n=-n
    leg=strip(g,q,u,A,q+n)&(g.Y>=0)&(g.Y<=Q+A)
    return (out|leg|rect(g,0,A,0,CAP))&~cnt
