"""Series E, page E-6 (figures 5 6 7 8 9 0)."""
import numpy as np
from fractions import Fraction as Fr
from dim1961 import Grid, raster
from dimE3 import strip
def inch(s): return float(sum(Fr(p) for p in s.split()))
R8={'A':'1 3/8','B':'6 3/8','C':'8 1/8','D':'2 7/8','E':'5 3/4','F':'3 3/4','G':'3/8','H':'7/8','I':'5/8','J':'1/8','K':'1/4',
    'L':'4 7/8','M':'1 1/8','N':'2 1/4','O':'9 3/4','P':'4 7/16','Q':'4 3/4','R':'2 5/8','S':'1 63/64','T':'4 5/8','U':'6 5/8',
    'V':'3 1/8','W':'1 3/16','X':'17/32'}
DIM={k:inch(v)*690/8 for k,v in R8.items()}
A=DIM['A']; CAP=690.0; OV=DIM['J']; TOP=CAP+OV; BOT=-OV
ring=lambda g,c,ro,ri: g.disk(c,ro)&~g.disk(c,ri)
deg=lambda v: np.degrees(np.arctan2(v[1],v[0]))
def fig5():
    g=Grid(); B,C,D,E,F,G,H,I,K=[DIM[k] for k in 'BCDEFGHIK']
    cD=np.array([B-D,BOT+E-D])                    # bowl: touches the right edge, top at E above the overshoot line
    cE=np.array([cD[0],BOT]); cC=np.array([cD[0],CAP])   # big arcs: top-left from below, bottom-left from above (both tangent to D)
    m=ring(g,cD,D,D-A)&(g.X>=cD[0])
    yn=CAP-F-G; xn=cE[0]-np.sqrt((E-A)**2-(yn-cE[1])**2)          # where the stem's bottom edge meets the counter
    p0=np.array([K,CAP-F])
    above=((g.Y-p0[1])*(xn-p0[0])-(g.X-p0[0])*(yn-p0[1]))>=0       # above the line from the stem's foot to the notch
    m|=ring(g,cE,E,E-A)&(g.X<=cE[0])&(g.Y>=yn)&((g.X>=xn)|above)
    xt=0.0; yt=cC[1]-np.sqrt(C*C-(xt-cC[0])**2)                   # lower terminal: tip on the left edge, cut along the radius
    m|=ring(g,cC,C,C-A)&g.sector(cC,deg(np.array([xt,yt])-cC),-90)
    p0=np.array([K,CAP-F]); p1=np.array([2*K,CAP]); u=(p1-p0)/np.linalg.norm(p1-p0)
    stem=strip(g,p0,u,A,(B,CAP))&(g.Y<=CAP)
    stem&=((g.Y-p0[1])*(xn-p0[0])-(g.X-p0[0])*(yn-p0[1]))>=0        # above the stem bottom edge
    m|=stem
    m&=~(strip(g,p0,u,1e4,(-1e3,CAP)) & ~strip(g,p0,u,0.0,(-1e3,CAP)) & (g.Y>=yn))  # nothing left of the stem's left edge
    m|=(g.X>=2*K)&(g.X<=B-I)&(g.Y>=CAP-A)&(g.Y<=CAP)
    return g.bitmap(m),{'terminal y':round(yt,1),'H from overshoot':round(BOT+H,1),'notch':(round(xn,1),round(yn,1))}
def fig7():
    g=Grid(); B,O,R=DIM['B'],DIM['O'],DIM['R']
    p=np.array([B-R,0.0]); q=np.array([B,CAP-A]); mid=(p+q)/2; ch=q-p; L=np.linalg.norm(ch)
    perp=np.array([ch[1],-ch[0]])/L; c=mid+perp*np.sqrt(O*O-(L/2)**2)
    if c[0]<mid[0]: c=mid-perp*np.sqrt(O*O-(L/2)**2)            # centre to the right
    stem=ring(g,c,O+A,O)&(g.Y>=0)&(g.Y<=CAP-A)
    return g.bitmap(stem|((g.X>=0)&(g.X<=B)&(g.Y>=CAP-A)&(g.Y<=CAP))),{'centre':c.round(1)}

import cv2
from dimE4 import arc
def poly_mask(g,poly):
    H,W=g.X.shape; m=np.zeros((H,W),np.uint8)
    pts=np.round(np.column_stack([np.array(poly)[:,0]-g.x0, g.y1-np.array(poly)[:,1]])*16).astype(np.int32)
    cv2.fillPoly(m,[pts],1,shift=4); return m.astype(bool)
def oval(ct,rt,cl,rl,cb,rb,cr,rr):
    """four-centre oval: arcs joined where neighbouring circles touch (internal tangency)"""
    def tp(c1,r1,c2,r2):
        big,rbig,small=(c1,r1,c2) if r1>=r2 else (c2,r2,c1)
        d=small-big; return big+rbig*d/np.linalg.norm(d)
    Ttl,Tbl,Tbr,Ttr=tp(ct,rt,cl,rl),tp(cb,rb,cl,rl),tp(cb,rb,cr,rr),tp(ct,rt,cr,rr)
    a=lambda c,p: np.degrees(np.arctan2(p[1]-c[1],p[0]-c[0]))
    return np.vstack([arc(ct,rt,a(ct,Ttr),a(ct,Ttl),True),arc(cl,rl,a(cl,Ttl),a(cl,Tbl),True),
                      arc(cb,rb,a(cb,Tbl),a(cb,Tbr),True),arc(cr,rr,a(cr,Tbr),a(cr,Ttr),True)])
def bowl6(g):
    B,L,N=DIM['B'],DIM['L'],DIM['N']; V=DIM['V']
    xc=B-V; yt=BOT+L; ym=(yt+BOT)/2                     # bowl: top at L above the bottom overshoot line
    ct=np.array([xc,BOT]); cb=np.array([xc,yt]); cl=np.array([N,ym]); cr=np.array([B-N,ym])
    out=poly_mask(g,oval(ct,L,cl,N,cb,L,cr,N)); inn=poly_mask(g,oval(ct,L-A,cl,N-A,cb,L-A,cr,N-A))
    return out,inn
def fig6():
    g=Grid(); B,L,M=DIM['B'],DIM['L'],DIM['M']
    out,inn=bowl6(g)
    # tail: a connecting arc (radius found to fit, as the 1961 notes allow): it touches the left edge and the bowl's
    # end arc N at the same point, so it is centred level with the N centres, and its flat end reaches the cap line at B-M
    xt=B-M; yN=(BOT+L+BOT)/2; dy=CAP-yN
    r=(xt*xt+dy*dy)/(2*xt); c=np.array([r,yN])
    a0=np.degrees(np.arctan2(CAP-yN,xt-r))
    tail=ring(g,c,r,r-A)&g.sector(c,a0,181)&(g.X<=xt)
    return g.bitmap((out|tail)&~inn),{'tail centre':c.round(1),'tail radius':round(r,1),'L':round(L,1)}
def fig9():
    (bm,xl,base),info=fig6(); g=Grid()
    # the 9 is the 6 turned through 180 degrees about the centre of its box
    H,W=bm.shape; ys,xs=np.nonzero(bm); x=xs-xl; y=base-ys
    cx=DIM['B']/2; cy=CAP/2; X=(2*cx-x+xl).astype(int); Y=(base-(2*cy-y)).astype(int)
    out=np.zeros_like(bm); ok=(X>=0)&(X<W)&(Y>=0)&(Y<H); out[Y[ok],X[ok]]=1
    return (out,xl,base),info
def fig0():
    g=Grid(); U,B,V=DIM['U'],DIM['B'],DIM['V']
    cl=np.array([B,CAP/2]); cr=np.array([U-B,CAP/2]); ct=np.array([U/2,TOP-V]); cb=np.array([U/2,BOT+V])
    out=poly_mask(g,oval(ct,V,cl,B,cb,V,cr,B)); inn=poly_mask(g,oval(ct,V-A,cl,B-A,cb,V-A,cr,B-A))
    return g.bitmap(out&~inn)
def fig8():
    g=Grid(); B,K,P,Q,R,S,T,H,W,X=[DIM[k] for k in 'BKPQRSTHWX']
    xc=B/2
    def stack(r,cy,rs,xs,which):
        """oval with arcs of radius r at top and bottom (mirror images) and side arcs rs centred at x=xs and B-xs.
        which='top': cy is the centre of the top arc (below it); which='bottom': cy is the centre of the bottom arc (above it)"""
        d=np.sqrt((r-rs)**2-(xc-xs)**2); ys=cy+d if which=='top' else cy-d
        ctop=np.array([xc,cy]) if which=='top' else np.array([xc,2*ys-cy])
        cbot=np.array([xc,2*ys-cy]) if which=='top' else np.array([xc,cy])
        return oval(ctop,r,np.array([B-xs,ys]),rs,cbot,r,np.array([xs,ys]),rs), ys
    # upper bowl: top arc P (centre below), sides S, K in from each side
    cpt=TOP-P; up_out,yu=stack(P,cpt,S,B-K-S,'top')
    up_in,yui=stack(P-A,cpt,X,B-K-A-X,'top')
    # lower bowl: bottom arc Q (centre above), sides R
    cqb=BOT+Q; lo_out,yl=stack(Q,cqb,R,B-R,'bottom')
    lo_in,yli=stack(Q-A,cqb,W,B-A-W,'bottom')
    m=(poly_mask(g,up_out)|poly_mask(g,lo_out))&~poly_mask(g,up_in)&~poly_mask(g,lo_in)
    return g.bitmap(m),{'upper side centre y':round(yu,1),'lower side centre y':round(yl,1)}
