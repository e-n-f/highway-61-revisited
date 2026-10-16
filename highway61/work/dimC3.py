"""Series C, page C-3 (M N O P Q R)."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
from dimE3 import strip
from dimC1 import _counter_end
R8={'A':'1 1/8','B':'5 3/16','C':'4 1/2','D':'4 3/8','E':'2 11/16','F':'7/8','G':'4 5/8','H':'2 1/32','I':'10','J':'1/8',
    'K':'2 1/4','L':'2 7/8','M':'1 7/8','N':'3 1/2','O':'3/8','P':'1 7/8','Q':'1','R':'7/16','S':'3 3/4','T':'2 1/8'}
R12={'A':'1 11/16','B':'7 25/32','C':'6 3/4','D':'6 9/16','E':'4 1/32','F':'1 5/16','G':'6 15/16','H':'3 3/64','I':'15',
     'J':'3/16','K':'3 3/8','L':'4 5/16','M':'2 13/16','N':'5 1/4','O':'9/16','P':'2 13/16','Q':'1 1/2','R':'21/32',
     'S':'5 5/8','T':'3 3/16'}
R18={'A':'2 17/32','B':'11 43/64','C':'10 1/8','D':'9 27/32','E':'6 1/16','F':'1 31/32','G':'10 13/32','H':'4 9/16',
     'I':'22 1/2','J':'9/32','K':'5 1/16','L':'6 15/32','M':'4 7/32','N':'7 7/8','O':'27/32','P':'4 7/32','Q':'2 1/4',
     'R':'1','S':'8 7/16','T':'4 25/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
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
def _O(g):
    """side arcs I centred at mid-height on the far side (width G); top/bottom arcs on the overshoot lines, refitted for
    tangency (the table's H misses by a little, as in Series D and E)"""
    W,R=DIM['G'],DIM['I']
    cl=np.array([R,CAP/2]); cr=np.array([W-R,CAP/2])
    dx=R-W/2; h=TOP-CAP/2; H=(R*R-dx*dx-h*h)/(2*(R-h))
    cu=np.array([W/2,TOP-H]); cd=np.array([W/2,BOT+H])
    def area(rs,rt):
        lens=g.disk(cl,rs)&g.disk(cr,rs)
        top=g.sector(cu,ang(cr,cu),ang(cl,cu))&~g.disk(cu,rt)&(g.Y>cu[1])
        bot=g.sector(cd,ang(cl,cd),ang(cr,cd))&~g.disk(cd,rt)&(g.Y<cd[1])
        return lens&~top&~bot
    return area(R,H)&~area(R-A,H-A),{'H table':round(DIM['H'],1),'H tangent':round(H,1)}
def letter_O():
    g=Grid(); return _O(g)
def letter_Q():
    """tail: rectangle F wide; its left edge from (G/2-R, BOT+P) to the bottom corner (G/2+Q, BOT-O)"""
    g=Grid(); m,info=_O(g); W,F,O,P,Q,R=[DIM[k] for k in 'GFOPQR']
    TL=np.array([W/2-R,BOT+P]); BL=np.array([W/2+Q,BOT-O]); u=(BL-TL)/np.linalg.norm(BL-TL)
    n=np.array([u[1],-u[0]])
    if n[0]<0: n=-n
    tail=strip(g,TL,u,F,TL+n)&strip(g,TL,n,np.linalg.norm(BL-TL),BL)
    return m|tail,info
def _bowl(g,D,top,bot,L,M):
    """flat-topped bowl (as the Series C B): side arc L centred on the bowl's mid-line touching x=D, corners M;
    counter = the bowl offset by A"""
    c=np.array([D-L,(top+bot)/2])
    out,_,_=_counter_end(g,0,top,bot,c,L,M)
    cnt,_,_=_counter_end(g,A,top-A,bot+A,c,L-A,M-A)
    return out&~cnt,c
def letter_P():
    g=Grid(); D,K,L,M,N=[DIM[k] for k in 'DKLMN']
    b,c=_bowl(g,D,CAP,N,L,M)
    return b|rect(g,0,A,0,CAP),{'bowl centre vs CAP-K':(round(c[1],1),round(CAP-K,1))}
def letter_R():
    """bowl as P, bottom at S; leg's right edge from the corner (D, 0) up to the bowl, A in from the right side"""
    g=Grid(); D,L,M,S,T=[DIM[k] for k in 'DLMST']
    b,c=_bowl(g,D,CAP,S,L,M)
    p=np.array([D-A,S]); q=np.array([D,0.0]); u=(p-q)/np.linalg.norm(p-q)
    n=np.array([-u[1],u[0]])
    if n[0]>0: n=-n
    leg=strip(g,q,u,A,q+n)&(g.Y>=0)&(g.Y<=S+A)
    cnt,_,_=_counter_end(g,A,CAP-A,S+A,c,L-A,M-A)
    return (b|leg|rect(g,0,A,0,CAP))&~cnt,{'bowl centre vs CAP-T':(round(c[1],1),round(CAP-T,1))}
