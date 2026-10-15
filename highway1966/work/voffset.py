import numpy as np,cv2,pickle
from fontTools.ttLib import TTFont
from fontTools.pens.basePen import BasePen
class Flat(BasePen):
    def __init__(s,gs): super().__init__(gs); s.cs=[]; s.cur=None
    def _moveTo(s,p): s.cur=[np.array(p,float)]
    def _lineTo(s,p): s.cur.append(np.array(p,float))
    def _qCurveToOne(s,a,b):
        p0=s.cur[-1]; a=np.array(a,float); b=np.array(b,float)
        L=np.hypot(*(a-p0))+np.hypot(*(b-a)); n=max(2,int(L/3))
        for t in np.linspace(0,1,n+1)[1:]: s.cur.append((1-t)**2*p0+2*(1-t)*t*a+t*t*b)
    def _curveToOne(s,a,b,c): raise NotImplementedError
    def _closePath(s):
        P=np.array(s.cur)
        if np.allclose(P[0],P[-1]): P=P[:-1]
        s.cs.append(P); s.cur=None
    _endPath=_closePath
def ymap(y,d):
    return np.where(y<0,y-d,np.where(y<=525,y*(525+2*d)/525-d,y+d))
def offset(P,d):
    n=len(P); out=np.zeros_like(P)
    # drop duplicate points
    keep=[0]
    for i in range(1,n):
        if np.hypot(*(P[i]-P[keep[-1]]))>0.05: keep.append(i)
    P=P[keep]; n=len(P)
    E=np.roll(P,-1,0)-P; L=np.hypot(E[:,0],E[:,1]); E=E/L[:,None]
    N=np.column_stack([E[:,1],-E[:,0]])        # right-hand normals (ink side in TrueType orientation)
    out=[]
    for i in range(n):
        n1=N[i-1]; n2=N[i]; c=1+np.dot(n1,n2)
        if c>=0.6: out.append(P[i]+d*(n1+n2)/c)       # miter
        else: out.append(P[i]+d*n1); out.append(P[i]+d*n2)   # bevel very sharp turns
    return np.array(out)
def thin_glyphs(fontpath,letters,d):
    f=TTFont(fontpath); gs=f.getGlyphSet(); res={}
    for l in letters:
        pen=Flat(gs); gs[l].draw(pen)
        polys=[]; origs=[]
        for P in pen.cs:
            Q=P.copy(); Q[:,1]=ymap(Q[:,1],d)
            origs.append(Q); polys.append(offset(Q,d))
        allp=np.vstack(polys); xmin,ymin=allp.min(0); xmax,ymax=allp.max(0)
        pad=30; W=int(xmax-xmin)+2*pad; H=int(ymax-ymin)+2*pad
        oy=ymax+pad   # row of y=0 is oy - 0
        img=np.zeros((H,W),np.uint8)
        pts=[np.round(np.column_stack([(p[:,0]-xmin+pad)*16,(oy-p[:,1])*16])).astype(np.int32) for p in polys]
        # even-odd fill: draw each contour with XOR
        for q in pts:
            m=np.zeros_like(img); cv2.fillPoly(m,[q],1,lineType=cv2.LINE_8,shift=4); img^=m
        org=np.zeros_like(img)
        for q in [np.round(np.column_stack([(p[:,0]-xmin+pad)*16,(oy-p[:,1])*16])).astype(np.int32) for p in origs]:
            m=np.zeros_like(img); cv2.fillPoly(m,[q],1,lineType=cv2.LINE_8,shift=4); org^=m
        ero=cv2.erode(org,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*d+1,2*d+1)))
        lim=cv2.dilate(ero,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(11,11)))
        img=img&lim
        res[l]=(img,pad,oy,'voffset')
    return res
if __name__=='__main__':
    d=9
    out=thin_glyphs('/mnt/user-data/outputs/Highway1966-E-Reference.ttf','abcdefghijklmnopqrstuvwxyz',d)
    pickle.dump(out,open('lcE_E.pkl','wb'))
    from src import stem_width
    for l,y in [('b',250),('n',250),('l',300),('o',250)]:
        b,ox,oy,_=out[l]; print(l,stem_width(b,oy,y/525))
    for l in 'hgx':
        b,ox,oy,_=out[l]; ys=np.nonzero(b.any(1))[0]; print(l,'top',round(oy-ys.min()),'bottom',round(oy-ys.max()))
