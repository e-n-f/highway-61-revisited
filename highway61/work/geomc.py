import numpy as np,cv2,pickle
from geomfit import M,ML
G=pickle.load(open('geomfit.pkl','rb'))
TM=140.0
def ring_params(po,pi):
    """centre-line extremes and vertical thicknesses of a ring, in y-up units"""
    return dict(Lc=((po[0]-po[2])+(pi[0]-pi[2]))/2, Rc=((po[0]+po[2])+(pi[0]+pi[2]))/2,
                cy=po[1], bo=po[3], top=(pi[1]-pi[3])-(po[1]-po[3]), bot=(po[1]+po[3])-(pi[1]+pi[3]))
def draw_ellipse(img,cx,cy,a,b,val):
    cv2.ellipse(img,(int(round(cx*16)),int(round(cy*16))),(max(1,int(round(a*16))),max(1,int(round(b*16)))),0,0,360,val,-1,lineType=cv2.LINE_8,shift=4)
def ring(img,rp,f,T,pad,oy_m,noy,clip=None):
    """draw an even-stroke elliptical ring; rows: master row r -> noy-(oy_m-r)"""
    L=pad+rp['Lc']*f; R=pad+rp['Rc']*f
    cy=noy-(oy_m-rp['cy'])          # rp['cy'] is a master row
    k=T/TM
    ao=(R-L)/2+T/2; ai=(R-L)/2-T/2; cx=(L+R)/2
    bo=rp['bo']; top=rp['top']*k; bot=rp['bot']*k
    bi=bo-(top+bot)/2; cyi=cy+(bot-top)/2*-1*-1  # shift so top/bottom thicknesses match
    cyi=cy+(top-bot)/2
    o=np.zeros_like(img); draw_ellipse(o,cx,cy,ao,bo,1)
    i=np.zeros_like(img); draw_ellipse(i,cx,cyi,ai,bi,1)
    return o&(1-i), dict(cx=cx,cy=cy,ao=ao,ai=ai,bo=bo,bi=bi,cyi=cyi)
def term_lines(l):
    out=[]
    for L in ML[l]:
        u=L['u']; ang=np.degrees(np.arctan2(abs(u[1]),abs(u[0])))
        if 3<ang<87: out.append(L)
    return out
def halfplane(shape,c,n,side):
    H,W=shape; yy,xx=np.mgrid[0:H,0:W]
    return (((xx-c[0])*n[0]+(yy-c[1])*n[1])*side>0)
def build(l,f,T,pad=40):
    b,ox,oym,_=M[l]
    H=1300; noy=900.0
    W=int(600*max(f,1)+2*pad+200)
    img=np.zeros((H,W),np.uint8)
    X=lambda x: pad+x*f; Rw=lambda yup: noy-yup
    if l in 'oce':
        po,pi,_,_=G['o']; rp=ring_params(po,pi)
        # align the letter's own left edge: shift the o ring to this letter's position
        xs=np.nonzero(b.any(0))[0]; dx=xs[0]-(po[0]-po[2])
        rp['Lc']+=dx; rp['Rc']+=dx
        r,g=ring(img,rp,f,T,pad,oym,noy); img|=r
        if l=='c':
            for Lt in term_lines('c'):
                c=np.array([X(Lt['c'][0]),noy-(oym-Lt['c'][1])])
                cut=halfplane(img.shape,c,Lt['n'],-Lt['side'])
                right=np.zeros_like(img,bool); right[:,int(g['cx']):]=True
                half=np.zeros_like(img,bool)
                if c[1]<g['cy']: half[:int(g['cy']),:]=True      # upper terminal acts on the upper half
                else: half[int(g['cy']):,:]=True
                img[cut&right&half]=0
        if l=='e':
            hs=sorted([oym-L['c'][1] for L in ML['e'] if abs(L['u'][1])<0.05])
            yb,yt=hs[0],hs[-1]; mid=(yb+yt)/2; th=(yt-yb)*T/TM
            yb2,yt2=mid-th/2,mid+th/2
            outer=np.zeros_like(img); draw_ellipse(outer,g['cx'],g['cy'],g['ao'],g['bo'],1)
            bar=np.zeros_like(img); bar[int(round(Rw(yt2))):int(round(Rw(yb2)))+1,:]=1
            img|=bar&outer
            Lt=term_lines('e')[0]
            c=np.array([X(Lt['c'][0]),noy-(oym-Lt['c'][1])])
            cut=halfplane(img.shape,c,Lt['n'],-Lt['side'])
            ap=np.zeros_like(img,bool); ap[int(round(Rw(yb2))):,int(g['cx']):]=True
            img[cut&ap]=0
            # opening between the bar and the terminal: keep the counter open to the right
            inner=np.zeros_like(img); draw_ellipse(inner,g['cx'],g['cyi'],g['ai'],g['bi'],1)
            img[(ap)&(inner>0)]=0
            img[int(round(Rw(yb2))):int(round(Rw(g['cy']-g['bo']+1))), int(g['cx']):][
                (inner[int(round(Rw(yb2))):int(round(Rw(g['cy']-g['bo']+1))), int(g['cx']):]>0)]=0
    else:
        po,pi,_,se=G[l]; rp=ring_params(po,pi)
        r,g=ring(img,rp,f,T,pad,oym,noy); img|=r
        stems=[x for x in se]
        sc=(min(stems)+max(stems))/2 if l in 'bp' else (min(stems)+max(stems))/2
        # the stem: centre mapped, constant width T
        scx=X(sc); x0,x1=scx-T/2,scx+T/2
        up= l in 'bd'
        TL=[L for L in term_lines(l) if abs(L['u'][0])>0.5]; Lt=(min if up else max)(TL,key=lambda L:L['c'][1])
        y_top=690 if up else 525; y_bot=0 if up else -280
        st=np.zeros_like(img); st[int(round(Rw(y_top+60))):int(round(Rw(y_bot-60))), int(round(x0)):int(round(x1))]=1
        if up: st[int(round(Rw(-1))):,:]=0
        else:  st[:int(round(Rw(525))),:]=0
        c=np.array([X(Lt['c'][0]),noy-(oym-Lt['c'][1])])
        st[halfplane(img.shape,c,Lt['n'],-Lt['side'])]=0
        img|=st
        # nothing below the baseline (b,d) or above the x-height (p,q) except the stem
    return img,pad,noy
if __name__=='__main__':
    import sys
    from sheet import sheet
    items=[]
    for l in 'ocebdpq':
        img,pad,noy=build(l,1.0,140); items.append((img,pad,noy))
        b,ox,oy,_=M[l]; items.append((b,ox,oy))
    sheet(items,'geom_test.png',0.25,gap=30)
def compare(l,f=1.0,T=140):
    img,pad,noy=build(l,f,T); b,ox,oy,_=M[l]
    # align by left ink edge and baseline
    A=np.zeros((1300,1200),np.uint8); B=np.zeros_like(A)
    def place(dst,src,oyy):
        ys,xs=np.nonzero(src); Y=(ys-oyy+900).astype(int); X=xs-xs.min()+50
        ok=(Y>=0)&(Y<1300)&(X<1200); dst[Y[ok],X[ok]]=1
    place(A,img,noy); place(B,b,oy)
    return (A^B).sum()/B.sum()
