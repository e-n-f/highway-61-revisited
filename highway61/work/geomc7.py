"""e: top curve and bottom curve follow two different ellipse pairs (orange/red), plus the crossbar."""
import numpy as np,cv2
import geomc3 as g3, geomc5 as g5, geomc6 as g6
from geomc3 import Canvas,master,TM,edge_pts
from geomfit import fit_super
b,oy=master('e')
H=sorted([oy-L['c'][1] for L in g3.ML['e'] if abs(L['u'][1])<0.05]); YB,YT=H[0],H[-1]   # crossbar bottom/top (y-up)
def polar(lo,hi,ylo,yhi):
    """outermost / innermost outline points by angle around the letter's centre, within a band of heights"""
    ys,xs=np.nonzero(b); ctr=np.array([(xs.min()+xs.max())/2,(ys.min()+ys.max())/2])
    cs,hh=cv2.findContours(b.copy(),cv2.RETR_CCOMP,cv2.CHAIN_APPROX_NONE)
    Pall=np.vstack([c[:,0,:] for c in cs]).astype(float)
    d=Pall-ctr; ang=np.degrees(np.arctan2(-d[:,1],d[:,0]))%360; rad=np.hypot(d[:,0],d[:,1]); yu=oy-Pall[:,1]
    O=[];I=[]
    for a0 in range(lo,hi,3):
        k=(ang>=a0)&(ang<a0+3)&(yu>=ylo)&(yu<=yhi)
        if k.sum()<2: continue
        idx=np.nonzero(k)[0]; O.append(Pall[idx[np.argmax(rad[idx])]]); I.append(Pall[idx[np.argmin(rad[idx])]])
    return np.array(O),np.array(I)
# top: everything above the crossbar, all the way round to the left side
Ot,_It=polar(0,360,YT+5,700); po_t,_=fit_super(Ot,2.0)
cs,hh=cv2.findContours(b.copy(),cv2.RETR_CCOMP,cv2.CHAIN_APPROX_NONE)
eye=max([c[:,0,:] for i,c in enumerate(cs) if hh[0][i][3]!=-1],key=len).astype(float)
pi_t,_=fit_super(eye[(oy-eye[:,1])>YT+4],2.0)
# bottom: below the crossbar, away from the aperture on the right
Ob,_Ib=polar(150,345,-40,YB-5); po_b,_=fit_super(Ob,2.0)
xs_=np.nonzero(b.any(0))[0]; x0,x1=xs_[0],xs_[-1]; xm=(x0+x1)/2
def rows(P): return np.column_stack([P[:,0],oy-P[:,1]])
Ib=np.vstack([rows(edge_pts('e',(x0+40,xm+40),(60,YB-5),'top')),rows(edge_pts('e',(x0+60,xm),(80,YB-5),'right'))])
pi_b,_=fit_super(Ib,2.0)
g6.PROF['e_top']=g5.profile(po_t,pi_t,oy); g6.PROF['e_bot']=g5.profile(po_b,pi_b,oy)
def build(l,f,T):
    C=Canvas(f,T); m=C.mask()>0
    rt,gt=g6.ring(C,'e_top',side_theta=np.pi,junction=np.pi)
    rb,gb=g6.ring(C,'e_bot',side_theta=np.pi,junction=np.pi)
    k=T/TM; mid=(YB+YT)/2; th=(YT-YB)*k; yb,yt=mid-th/2,mid+th/2
    # line the bottom curve's left side up with the top curve's where they meet (bottom of the crossbar)
    r_=int(round(C.noy-yb)); lt=np.nonzero(rt[r_])[0]; lb=np.nonzero(rb[r_])[0]
    if len(lt) and len(lb): rb=np.roll(rb,int(lt[0]-lb[0]),axis=1)
    m|=rt&(C.yu>=yb); m|=rb&(C.yu<yb)        # the top curve runs down to the bottom of the crossbar
    # crossbar: from the eye's left side out to the top curve's outer edge
    outer=C.mask(); Q=gt['Q']; cv2.fillPoly(outer,[np.round(np.column_stack([Q[:,0],C.noy-Q[:,1]])*16).astype(np.int32)],1,shift=4)
    m|=(C.yu>=yb)&(C.yu<=yt)&(outer>0)
    # no bulge past the crossbar's end above it: the right side rises from the bar end

    # aperture: below the bar on the right, cut by the 1966 terminal (placed at the same point around the bottom curve)
    pr=g6.PROF['e_bot']
    TL=[L for L in g3.ML['e'] if 3<np.degrees(np.arctan2(abs(L['u'][1]),abs(L['u'][0])))<87]
    L=TL[0]; px,py=L['c'][0],oy-L['c'][1]
    thp=np.arctan2((py-pr['ci'][1])/pr['bi'],(px-pr['ci'][0])/pr['ai'])
    p=np.array([gb['cx']+gb['ai']*np.cos(thp),gb['cy']+gb['bi']*np.sin(thp)])
    slope=(-L['u'][1])/L['u'][0]
    above=C.yu>p[1]+slope*(C.xx-p[0])
    m&=~((C.xx>gb['cx'])&(C.yu<yb)&above)
    return m.astype(np.uint8),C.pad,C.noy
