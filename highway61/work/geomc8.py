"""g: q-style bowl + stem + a tail made from its own outer/inner ellipse pair."""
import numpy as np,cv2
import geomc3 as g3, geomc5 as g5, geomc6 as g6
from geomc3 import Canvas,master,runs_at,sc,TM,term_line
from geomfit import fit_super
b,oy=master('g'); ys,xs=np.nonzero(b); yu=oy-ys
s=runs_at('g',300)[-1]; XS=sc(s)
cs,hh=cv2.findContours(b.copy(),cv2.RETR_CCOMP,cv2.CHAIN_APPROX_NONE)
outer=max([c[:,0,:] for i,c in enumerate(cs) if hh[0][i][3]==-1],key=len).astype(float)
tl=term_line('g','bottom')
# tail outline points below the baseline and right of the terminal, split by radius around a rough centre
pts=outer[(oy-outer[:,1])<-15]
P_=np.column_stack([pts[:,0],oy-pts[:,1]])
uu=tl['u']; nn=np.array([-uu[1],uu[0]])
side=np.sign((np.array([XS,-150])-tl['p'])@nn)
pts=pts[((P_-tl['p'])@nn)*side>8]
ctr=np.array([XS-150, oy-(-40)])
d=pts-ctr; ang=np.degrees(np.arctan2(-d[:,1],d[:,0]))%360; rad=np.hypot(d[:,0],d[:,1])
O=[];I=[]
for a0 in range(0,360,4):
    k=(ang>=a0)&(ang<a0+4)
    if k.sum()<2: continue
    idx=np.nonzero(k)[0]
    if rad[idx].max()-rad[idx].min()<40: continue            # need both edges of the stroke in this direction
    O.append(pts[idx[np.argmax(rad[idx])]]); I.append(pts[idx[np.argmin(rad[idx])]])
po,eo=fit_super(np.array(O),2.0)
# inner ellipse: tangent to the outer at its upper-left (as in the drawing), scaled down; fit the scale to the tail's inner edge
from geomc3 import edge_pts
from scipy.optimize import least_squares
Pin=edge_pts('g',(tl['p'][0]+20,XS-80),(-260,-30),'top')         # (x, y-up) of the tail's upper/inner edge
ocx,ocy,oa,ob=po[0],oy-po[1],po[2],po[3]
def inner_of(sx,sy):
    ai,bi=oa*sx,ob*sy; return (ocx-oa)+ai,(ocy+ob)-bi,ai,bi
def res(p):
    cx,cy,a,bb=inner_of(*p); return (np.sqrt(((Pin[:,0]-cx)/a)**2+((Pin[:,1]-cy)/bb)**2)-1)*min(a,bb)
sx,sy=least_squares(res,[0.76,0.78],bounds=([0.4,0.4],[0.95,0.95])).x
icx,icy,ia,ib=inner_of(sx,sy)
pi=np.array([icx,oy-icy,ia,ib,2.0]); ei=(sx,sy)
g6.PROF['g_tail']=g5.profile(po,pi,oy)
def build(l,f,T):
    C=Canvas(f,T); m=C.mask()>0
    bw,_=g5.bowl(C,'q',dx=XS-sc(runs_at('q',300)[-1])); m|=bw&(C.yu>=-5)
    r,gg=g6.ring(C,'g_tail',side_theta=0.0,junction=0.0)
    keep=r&(C.yu<0)&(C.xx<=C.X(XS)+T/2)
    P0=np.array([C.X(tl['p'][0]),tl['p'][1]]); u=tl['u']; n=np.array([-u[1],u[0]])
    sd=(C.xx-P0[0])*n[0]+(C.yu-P0[1])*n[1]; ref=(C.X(XS)-P0[0])*n[0]+(-150-P0[1])*n[1]
    keep&=sd*np.sign(ref)>0
    m|=keep
    ocy=(gg['Q'][:,1].max()+gg['Q'][:,1].min())/2
    m|=C.stem(XS,min(ocy,0),525)
    return m.astype(np.uint8),C.pad,C.noy
