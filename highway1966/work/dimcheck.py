import numpy as np,cv2,pickle,sys
import dim1961 as d
CE=pickle.load(open('caps_E.pkl','rb'))
def overlay(ch,polys):
    img,xl,base=d.raster(polys)
    b,ox,oy,_=CE[ch]
    H,W=900,900; A_=np.zeros((H,W),np.uint8); B_=np.zeros_like(A_)
    def place(dst,src,oyy):
        ys,xs=np.nonzero(src); Y=(ys-oyy+760).astype(int); X=(xs-xs.min()+80).astype(int)
        ok=(Y>=0)&(Y<H)&(X<W); dst[Y[ok],X[ok]]=1
    place(A_,img,base); place(B_,b,oy)
    T=np.full((H,W,3),255,np.uint8); T[B_>0]=(185,185,185)
    cs,_=cv2.findContours(A_,cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE); cv2.drawContours(T,cs,-1,(0,0,220),3)
    e=(A_^B_).sum()/B_.sum()
    ysA=np.nonzero(A_.any(1))[0]; ysB=np.nonzero(B_.any(1))[0]; xsA=np.nonzero(A_.any(0))[0]; xsB=np.nonzero(B_.any(0))[0]
    info='w %d/%d h %d/%d'%(np.ptp(xsA)+1,np.ptp(xsB)+1,np.ptp(ysA)+1,np.ptp(ysB)+1)
    cv2.putText(T,'%s %.1f%% %s'%(ch,100*e,info),(10,860),cv2.FONT_HERSHEY_SIMPLEX,0.9,(0,0,0),2)
    return T,e,info
if __name__=='__main__':
    tiles=[]
    p,_=d.Z(); q,_=d.one()
    for ch,polys in [('Y',d.Y()),('Z',p),('1',q),('4',d.four())]:
        T,e,info=overlay(ch,polys); print(ch,'%.3f'%e,info); tiles.append(T)
    cv2.imwrite('dimchk.png',cv2.resize(np.hstack(tiles),None,fx=0.45,fy=0.45))
def overlay_bm(ch,bm,xl,base):
    b,ox,oy,_=CE[ch]
    H,W=900,900; A_=np.zeros((H,W),np.uint8); B_=np.zeros_like(A_)
    def place(dst,src,oyy):
        ys,xs=np.nonzero(src); Y=(ys-oyy+760).astype(int); X=(xs-xs.min()+80).astype(int)
        ok=(Y>=0)&(Y<H)&(X<W); dst[Y[ok],X[ok]]=1
    place(A_,bm,base); place(B_,b,oy)
    T=np.full((H,W,3),255,np.uint8); T[B_>0]=(185,185,185)
    cs,_=cv2.findContours(A_,cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE); cv2.drawContours(T,cs,-1,(0,0,220),3)
    e=(A_^B_).sum()/B_.sum(); cv2.putText(T,'%s %.1f%%'%(ch,100*e),(10,860),cv2.FONT_HERSHEY_SIMPLEX,1,(0,0,0),2)
    return T,e
