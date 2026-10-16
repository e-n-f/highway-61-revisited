"""t, f, a from the user's ellipse markups (converted from screenshot pixels to font units)."""
import numpy as np,cv2
import geomc3 as g3, geomc5 as g5
from geomc3 import Canvas,master,runs_at,sc,TM,term_line
def conv(X0,Y0,img_x,mst_x,img_y,mst_y):
    sx=(mst_x[1]-mst_x[0])/(img_x[1]-img_x[0]); sy=(mst_y[1]-mst_y[0])/(img_y[1]-img_y[0])
    def f(cx,cy,a,b): return (mst_x[0]+(cx-img_x[0])*sx, mst_y[0]+(cy-img_y[0])*sy, a*abs(sx), b*abs(sy))
    return f
# screenshot ellipses (centre x, centre y, half-width, half-height) in pixels
tc=conv(0,0,(400,645),(124,263),(110,1320),(702.6,3.6))
T_RED=tc(632.5,715,357.5,615); T_BLUE=tc(647.5,712.5,472.5,612.5); T_ORG=tc(750,993,118,112)
fc=conv(0,0,(470,700),(95,229),(110,1310),(697.3,4.3))
F_OUT=fc(737.5,282,267.5,180); F_IN=fc(789,386,89,74)
ac=conv(0,0,(880,1118),(336,469),(110,1025),(533.9,2.9))
A_TOP=ac(730,422.5,397,325.5); A_BIG=ac(730.5,754,387.5,436); A_TIN=ac(747.5,424,142.5,91)
def E(C,e,role):
    """ellipse mask in canvas C. role 'outer' keeps the pen side at T; 'inner' is a counter."""
    cx,cy,a,b=e; f=C.f; T=C.T; k=T/TM
    if role=='outer': a2=f*(a-TM/2)+T/2; b2=b
    else:             a2=max(3,f*(a+TM/2)-T/2); b2=b
    img=C.mask(); cv2.ellipse(img,(int(round(C.X(cx)*16)),int(round(C.R(cy)*16))),(int(round(a2*16)),int(round(b2*16))),0,0,360,1,-1,shift=4)
    return img>0, (C.X(cx),cy,a2,b2)
def thin_top(inner,outer,k):
    """move an inner ellipse's top so the stroke above it is k times as thick; its bottom stays"""
    cx,cy,a,b=inner; ot=outer[1]+outer[3]; it=cy+b; bot=cy-b
    it2=ot-(ot-it)*k; return (cx,(it2+bot)/2,a,(it2-bot)/2)
def thin_bottom(inner,outer,k):
    cx,cy,a,b=inner; ob=outer[1]-outer[3]; ib=cy-b; top=cy+b
    ib2=ob+(ib-ob)*k; return (cx,(top+ib2)/2,a,(top-ib2)/2)
def keep_side(C,tl,ref_xy):
    P=np.array([C.X(tl['p'][0]),tl['p'][1]]); u=tl['u']; n=np.array([-u[1],u[0]])
    sd=(C.xx-P[0])*n[0]+(C.yu-P[1])*n[1]; ref=(ref_xy[0]-P[0])*n[0]+(ref_xy[1]-P[1])*n[1]
    return sd*np.sign(ref)>0
def build(l,f,T):
    C=Canvas(f,T); m=C.mask()>0; k=T/TM; ASC=700
    if l=='t':
        s=max(runs_at('t',250),key=lambda r:r[1]-r[0]); xs=sc(s)
        red,er=E(C,T_RED,'outer'); blue,eb=E(C,T_BLUE,'outer'); org,eo=E(C,thin_bottom(T_ORG,T_BLUE,k),'inner')
        stem=C.stem(xs,-5,ASC,top=term_line('t','top'))
        tail=blue&~org&(C.xx>=C.X(xs))&(C.yu<=eo[1])
        tl=dict(term_line('t','bottom')); tl['p']=tl['p']-np.array([(TM-T)/3/max(C.f,1e-6),0])   # cut moves left as the pen thins
        tail&=keep_side(C,tl,(C.X(xs),eo[1]-eo[3]))
        body=stem|tail
        body&=red|(C.yu>er[1])|(C.xx>er[0])                 # the red ellipse rounds the stem's lower left
        body&=(C.yu>=0)
        br=max(runs_at('t',460),key=lambda r_:r_[1]-r_[0]); th=131*k
        m|=body|((C.yu>=525-th)&(C.yu<=525)&(C.xx>=C.X(br[0]))&(C.xx<=C.X(br[1])))
    elif l=='f':
        s=max(runs_at('f',250),key=lambda r:r[1]-r[0]); xs=sc(s)
        out,eo=E(C,F_OUT,'outer'); inn,ei=E(C,thin_top(F_IN,F_OUT,k),'inner')
        br=max(runs_at('f',455),key=lambda r_:r_[1]-r_[0]); th=141*k
        hook=out&~inn&(C.xx<=C.X(br[1]))&(C.xx>=C.X(xs)-T/2)&(C.yu>=525-th)
        m|=hook|C.stem(xs,0,eo[1])|((C.yu>=525-th)&(C.yu<=525)&(C.xx>=C.X(br[0]))&(C.xx<=C.X(br[1])))
    elif l=='a':
        s=runs_at('a',150)[-1]; xs=sc(s)
        top,et=E(C,A_TOP,'outer'); big,eb=E(C,thin_top(A_BIG,A_TOP,k),'outer'); tin,ei=E(C,thin_top(A_TIN,A_TOP,k),'inner')
        arc=top&~big&~tin&(C.xx<=C.X(xs)+T/2)
        arc&=keep_side(C,term_line('a','top'),(C.X(xs),et[1]+et[3]-10))
        bw,_=g5.bowl(C,'a',right_to=C.X(xs)-T/2+T*0.25)
        m|=arc|(bw&(C.xx<=C.X(xs)+T/2))|C.stem(xs,0,et[1])
        # stem top joins the arc: fill between stem and arc outline
        m|=top&(C.xx>=C.X(xs)-T/2)&(C.xx<=C.X(xs)+T/2)&(C.yu>=0)
    return m.astype(np.uint8),C.pad,C.noy
