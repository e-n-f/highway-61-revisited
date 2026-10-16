"""Series A, page A-3 (M N O P Q R).  The constructions of Series B (dimB3) with this page's letters."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
from dimC1 import _counter_end
R8={'A':'3/4','B':'3 1/2','C':'5 1/4','D':'3','E':'2 5/8','F':'5/8','G':'3 1/4','H':'1 3/8','I':'18','J':'1/8',
    'K':'2 1/4','L':'1 3/4','M':'1 1/4','N':'5/16','O':'3/8','P':'1'}
R12={'A':'1 1/8','B':'5 1/4','C':'7 7/8','D':'4 1/2','E':'3 15/16','F':'15/16','G':'4 7/8','H':'2 1/16','I':'27',
     'J':'3/16','K':'3 3/8','L':'2 5/8','M':'1 7/8','N':'15/32','O':'9/16','P':'1 1/2'}
R18={'A':'1 11/16','B':'7 7/8','C':'11 13/16','D':'6 3/4','E':'5 29/32','F':'1 13/32','G':'7 5/16','H':'3 3/32',
     'I':'40 1/2','J':'9/32','K':'5 1/16','L':'3 15/16','M':'2 13/16','N':'11/16','O':'27/32','P':'2 1/4'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
def letter_M():
    """width B; the diagonals' lower edges meet at (B/2, L) and reach the stems' inner edges at height C"""
    g=Grid(); B,C,L=DIM['B'],DIM['C'],DIM['L']
    stems=rect(g,0,A,0,CAP)|rect(g,B-A,B,0,CAP)
    v=np.array([B/2,L]); pr=np.array([B-A,C]); u=(pr-v)/np.linalg.norm(pr-v)
    right=strip(g,v,u,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X<=B-A+1)&(g.X>=B/2)
    ul=np.array([-u[0],u[1]]); left=strip(g,v,ul,A,(B/2,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=B/2)
    # where the diagonal's upper edge reaches the cap line (the drawing has it at the stem's inner corner, x = A)
    n=np.array([-u[1],u[0]]); n=n if n[1]>0 else -n; q=v+n*A; xcap=q[0]+u[0]*(CAP-q[1])/u[1]
    return stems|right|left,{'upper edge at cap: x from the right stem':round(B-xcap,1),'A':round(A,1)}
def letter_N():
    """width D; the diagonal's lower edge from the left stem's inner edge E below the cap to (D-F, 0)"""
    g=Grid(); D,E,F=DIM['D'],DIM['E'],DIM['F']
    stems=rect(g,0,A,0,CAP)|rect(g,D-A,D,0,CAP)
    p0=np.array([A,CAP-E]); p1=np.array([D-F,0.0]); u=(p1-p0)/np.linalg.norm(p1-p0)
    diag=strip(g,p0,u,A,(D,CAP))&(g.Y>=0)&(g.Y<=CAP)&(g.X>=A-1)&(g.X<=D-A+1)
    below=((g.X-p0[0])*u[1]-(g.Y-p0[1])*u[0])>=0
    return (stems&~(below&(g.X>=D/2)))|diag
def _O(g):
    """side arcs I centred at mid-height on the far side (width G); top/bottom arcs refitted for tangency (table H)"""
    W,R=DIM['G'],DIM['I']
    cl=np.array([R,CAP/2]); cr=np.array([W-R,CAP/2])
    dx=R-W/2; h=TOP-CAP/2; G=(R*R-dx*dx-h*h)/(2*(R-h))
    cu=np.array([W/2,TOP-G]); cd=np.array([W/2,BOT+G])
    def area(rs,rt):
        lens=g.disk(cl,rs)&g.disk(cr,rs)
        top=g.sector(cu,ang(cr,cu),ang(cl,cu))&~g.disk(cu,rt)&(g.Y>cu[1])
        bot=g.sector(cd,ang(cl,cd),ang(cr,cd))&~g.disk(cd,rt)&(g.Y<cd[1])
        return lens&~top&~bot
    return area(R,G)&~area(R-A,G-A),{'H table':round(DIM['H'],1),'H tangent':round(G,1)}
def letter_O():
    g=Grid(); return _O(g)
def letter_Q():
    """tail: rectangle F wide; its left edge from (G/2-N, BOT+M) to the bottom corner (G/2+A, BOT-O)"""
    g=Grid(); m,info=_O(g); W,F,M,N,O=[DIM[k] for k in 'GFMNO']
    TL=np.array([W/2-N,BOT+M]); BL=np.array([W/2+A,BOT-O]); u=(BL-TL)/np.linalg.norm(BL-TL)
    n=np.array([u[1],-u[0]])
    if n[0]<0: n=-n
    tail=strip(g,TL,u,F,TL+n)&strip(g,TL,n,np.linalg.norm(BL-TL),BL)
    return m|tail,info
def _bowl(g,top,bot):
    """flat-topped bowl D wide: side arc D centred on the stem's left edge on the centre line K below the cap,
    corners L; counter = the bowl offset by A"""
    D,K,L=DIM['D'],DIM['K'],DIM['L']
    c=np.array([0.0,CAP-K])
    out,_,_=_counter_end(g,0,top,bot,c,D,L)
    cnt,_,_=_counter_end(g,A,top-A,bot+A,c,D-A,L-A)
    return out,cnt,c
def letter_P():
    """bowl bottom B above the baseline"""
    g=Grid(); out,cnt,c=_bowl(g,CAP,DIM['B'])
    return (out|rect(g,0,A,0,CAP))&~cnt
def letter_R():
    """bowl as P; leg's right edge from the corner (D, 0) up to (D-P, B)"""
    g=Grid(); B,D,P=DIM['B'],DIM['D'],DIM['P']
    out,cnt,c=_bowl(g,CAP,B)
    p=np.array([D-P,B]); q=np.array([D,0.0]); u=(p-q)/np.linalg.norm(p-q)
    n=np.array([-u[1],u[0]])
    if n[0]>0: n=-n
    leg=strip(g,q,u,A,q+n)&(g.Y>=0)&(g.Y<=B+A)
    return (out|leg|rect(g,0,A,0,CAP))&~cnt
