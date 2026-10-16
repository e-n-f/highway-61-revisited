import numpy as np,cv2,sys
import geomc7 as g6, geomc3 as g3
def ov(l):
    img,pad,noy=g6.build(l,1.0,140); b,oy=g3.master(l)
    H,W=1050,900; A=np.zeros((H,W),np.uint8); B=np.zeros_like(A)
    def place(dst,src,oyy,dx):
        ys,xs=np.nonzero(src); Y=(ys-oyy+720).astype(int); X=(xs+dx).astype(int)
        ok=(Y>=0)&(Y<H)&(X>=0)&(X<W); dst[Y[ok],X[ok]]=1
    place(B,b,oy,60); place(A,img,noy,60-pad)
    T_=np.full((H,W,3),255,np.uint8); T_[B>0]=(190,190,190)
    cs,_=cv2.findContours(A,cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE); cv2.drawContours(T_,cs,-1,(0,0,220),3)
    return T_,(A^B).sum()/max(1,B.sum())
tiles=[]
for l in sys.argv[1]:
    t,e=ov(l); print(l,'%.3f'%e); tiles.append(t)
cv2.imwrite('g7chk.png',cv2.resize(np.hstack(tiles),None,fx=0.4,fy=0.4))
