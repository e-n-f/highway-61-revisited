"""Series F, page F-1 (A B C D E F)."""
import numpy as np
from vgeom import VGrid as Grid
from dimtab import units, check
from dimD1 import rect, ring, ang
R8={'A':'1 1/2','B':'9','C':'1 3/4','D':'7 3/8','E':'2 1/8','F':'2 1/4','G':'1 1/4','H':'7/8','I':'1 1/8','J':'4 1/4',
    'K':'3/8','L':'4 1/8','M':'1/8','N':'4','O':'6 1/4','P':'6 3/8','Q':'3 1/2','R':'2 3/4','S':'2 7/8'}
R12={'A':'2 1/4','B':'13 1/2','C':'2 5/8','D':'11 1/16','E':'3 3/16','F':'3 3/8','G':'1 7/8','H':'1 5/16','I':'1 11/16',
     'J':'6 3/8','K':'9/16','L':'6 3/16','M':'3/16','N':'6','O':'9 3/8','P':'9 9/16','Q':'5 1/4','R':'4 1/8','S':'4 5/16'}
R18={'A':'3 3/8','B':'20 1/4','C':'3 15/16','D':'16 19/32','E':'4 25/32','F':'5 1/16','G':'2 13/16','H':'1 31/32',
     'I':'2 17/32','J':'9 9/16','K':'27/32','L':'9 9/32','M':'9/32','N':'9','O':'14 1/16','P':'14 11/32','Q':'7 7/8',
     'R':'6 3/16','S':'6 15/32'}
DIM=units(R18,18); CHECK=check({8:R8,12:R12,18:R18})
A=DIM['A']; CAP=690.0; OV=DIM['M']; TOP=CAP+OV; BOT=-OV
def letter_A():
    g=Grid(); B,C=DIM['B'],DIM['C']; xl,xr=(B-A)/2,(B+A)/2
    def edge(p0,p1): d=np.array(p1,float)-np.array(p0,float); return np.array(p0,float),d/np.linalg.norm(d)
    pL,dL=edge((0,0),(xl,CAP)); pR,dR=edge((B,0),(xr,CAP))
    side=lambda p,d: (g.X-p[0])*d[1]-(g.Y-p[1])*d[0]
    outer=(side(pL,dL)>=0)&(side(pR,dR)<=0)&(g.Y>=0)&(g.Y<=CAP)
    nL=np.array([dL[1],-dL[0]]); nR=np.array([-dR[1],dR[0]])
    inner=(side(pL+nL*A,dL)>=0)&(side(pR+nR*A,dR)<=0)
    return outer&~(inner&((g.Y<C)|(g.Y>C+A)))
def letter_B():
    """as the Series D B: outer bowls E (upper, K in from the right) and F (lower), tangent to the cap line and the
    baseline, crossing at the crotch; middle bar A thick centred on J.  Counter ends: arcs H (upper) and G (lower) level
    with the outer centres (tangent to the middle bar), then the outer arcs' offsets E-A, F-A towards the outer bars"""
    g=Grid(); D,E,F,G,H,I,J,K=[DIM[k] for k in 'DEFGHIJK']
    cu=np.array([D-K-E,CAP-E]); cl=np.array([D-F,F])
    # crotch: the outer circles' crossing on the right
    d=cu-cl; dist=np.linalg.norm(d); a=(F*F-E*E+dist*dist)/(2*dist); h=np.sqrt(F*F-a*a)
    base=cl+d/dist*a; nv=np.array([d[1],-d[0]])/dist; crotch=base+nv*h if (base+nv*h)[0]>base[0] else base-nv*h
    outer=rect(g,0,min(cu[0],cl[0]),0,CAP)|rect(g,0,cu[0],cu[1],CAP)|rect(g,0,cl[0],0,cl[1])
    outer|=(g.disk(cu,E)&(g.Y>=crotch[1]))|(g.disk(cl,F)&(g.Y<=crotch[1]))
    outer|=rect(g,0,crotch[0],cl[1],cu[1])
    ut,ub=CAP-A,J+A/2; lt,lb=J-A/2,A
    ku=np.array([cu[0]+E-A-H,cu[1]]); kl=np.array([cl[0]+F-A-G,cl[1]])
    cnt=rect(g,A,ku[0],ub,cu[1])|(g.disk(ku,H)&(g.X>=ku[0])&(g.Y<=cu[1])&(g.Y>=ub))
    cnt|=rect(g,A,cu[0],cu[1],ut)|(g.disk(cu,E-A)&(g.X>=cu[0])&(g.Y>=cu[1]))
    cnt|=rect(g,A,kl[0],cl[1],lt)|(g.disk(kl,G)&(g.X>=kl[0])&(g.Y>=cl[1])&(g.Y<=lt))
    cnt|=rect(g,A,cl[0],lb,cl[1])|(g.disk(cl,F-A)&(g.X>=cl[0])&(g.Y<=cl[1]))
    return outer&~cnt,{'crotch':crotch.round(1),'table (D-I, J)':(round(D-I,1),round(J,1)),
                       'H tangent to bar':(round(cu[1]-H,1),round(ub,1)),'G tangent to bar':(round(cl[1]+G,1),round(lt,1))}
def letter_C():
    """a ring of radius L (half the overshoot-to-overshoot height) centred L from the left at mid-height; terminals
    where it reaches x=D, cut along the radius"""
    g=Grid(); D,L=DIM['D'],DIM['L']
    c=np.array([L,CAP/2]); dy=np.sqrt(L*L-(D-L)**2)
    return ring(g,c,L,L-A)&g.sector(c,ang(c,(D,c[1]+dy)),ang(c,(D,c[1]-dy))),{'L vs (TOP-BOT)/2':(round(L,1),round((TOP-BOT)/2,1))}
def letter_D():
    """semicircular bowl N (= half the cap height) centred D-N from the left"""
    g=Grid(); D,N=DIM['D'],DIM['N']
    c=np.array([D-N,CAP/2])
    body=lambda d: rect(g,d,c[0],d,CAP-d)|(g.disk(c,N-d)&(g.X>=c[0]))
    return body(0)&~body(A),{'N vs CAP/2':(round(N,1),CAP/2)}
def letter_E():
    g=Grid(); O,P,Q,R=DIM['O'],DIM['P'],DIM['Q'],DIM['R']
    return rect(g,0,A,0,CAP)|rect(g,0,O,CAP-A,CAP)|rect(g,0,A+R,Q,Q+A)|rect(g,0,P,0,A)
def letter_F():
    g=Grid(); P,Q,S=DIM['P'],DIM['Q'],DIM['S']
    return rect(g,0,A,0,CAP)|rect(g,0,P,CAP-A,CAP)|rect(g,0,A+S,Q,Q+A)
