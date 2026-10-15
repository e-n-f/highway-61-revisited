import numpy as np,cv2,pickle
from geomfit import M,ML
G=pickle.load(open('geomfit.pkl','rb'))
def cl_params(po,pi):
    """centre-line of a ring (rows for y): left/right x, centre row, and outer half-height"""
    return dict(Lc=((po[0]-po[2])+(pi[0]-pi[2]))/2, Rc=((po[0]+po[2])+(pi[0]+pi[2]))/2, cy=po[1], bo=po[3])
def brush_ring(shape,cx,cy,ac,bc,T):
    img=np.zeros(shape,np.uint8)
    cv2.ellipse(img,(int(round(cx*16)),int(round(cy*16))),(int(round(ac*16)),int(round(bc*16))),0,0,360,1,
                thickness=max(1,int(round(T))),lineType=cv2.LINE_8,shift=4)
    return img
def fill_ellipse(shape,cx,cy,a,b):
    img=np.zeros(shape,np.uint8)
    cv2.ellipse(img,(int(round(cx*16)),int(round(cy*16))),(int(round(a*16)),int(round(b*16))),0,0,360,1,-1,shift=4)
    return img
def term_lines(l):
    out=[]
    for L in ML[l]:
        u=L['u']; ang=np.degrees(np.arctan2(abs(u[1]),abs(u[0])))
        if 3<ang<87 and abs(u[0])>0.5: out.append(L)
    return out
def build(l,f,T,pad=40):
    b,ox,oym,_=M[l]
    H=1300; noy=900.0; W=int(600*max(f,1)+2*pad+200); shape=(H,W)
    X=lambda x: pad+x*f; R=lambda yup: noy-yup
    yy,xx=np.mgrid[0:H,0:W]; yup=noy-yy
    if l in 'oce':
        po,pi,_,_=G['o']; c=cl_params(po,pi)
        xs=np.nonzero(b.any(0))[0]; dx=xs[0]-(po[0]-po[2])
        L,Rr=X(c['Lc']+dx),X(c['Rc']+dx); cx=(L+Rr)/2; ac=(Rr-L)/2
        cy=noy-(oym-c['cy']); bc=c['bo']-T/2
        img=brush_ring(shape,cx,cy,ac,bc,T)
        cyu=noy-cy
        if l in 'ce':
            # the e's lower terminal (kept at its 1966 angle) opens the bowl; the c uses it and its mirror image
            Lt=term_lines('e')[0]; eb,_,eoy,_=M['e']
            p=Lt['c']; u=Lt['u']
            exs=np.nonzero(eb.any(0))[0]; edx=xs[0]-exs[0]
            P=np.array([X(p[0]+edx), eoy-p[1]])          # point (x, yup)
            m=-u[1]/u[0]                                   # slope in y-up terms
            line=lambda x: P[1]+m*(x-P[0])
            img[(xx>cx)&(yup<cyu)&(yup>line(xx))]=0
            if l=='c':
                linem=lambda x: 2*cyu-line(x)
                img[(xx>cx)&(yup>cyu)&(yup<linem(xx))]=0
            if l=='e':
                hs=sorted([oym-L_['c'][1] for L_ in ML['e'] if abs(L_['u'][1])<0.05])
                yb,yt=hs[0],hs[-1]; mid=(yb+yt)/2; th=T*0.95
                outer=fill_ellipse(shape,cx,cy,ac+T/2,bc+T/2)
                img[(yup>=mid-th/2)&(yup<=mid+th/2)&(outer>0)]=1
                img[(yup<mid-th/2)&(yup>line(xx))&(xx>cx)]=0
        return img,pad,noy
    po,pi,_,se=G[l]; c=cl_params(po,pi)
    L,Rr=X(c['Lc']),X(c['Rc']); cx=(L+Rr)/2; ac=(Rr-L)/2
    cy=noy-(oym-c['cy']); bc=c['bo']-T/2
    img=brush_ring(shape,cx,cy,ac,bc,T)
    sc=(min(se)+max(se))/2; scx=X(sc); x0,x1=scx-T/2,scx+T/2
    up=l in 'bd'
    TL=term_lines(l); Lt=(min if up else max)(TL,key=lambda L_:L_['c'][1])
    p=Lt['c']; u=Lt['u']; P=np.array([X(p[0]), oym-p[1]]); m=-u[1]/u[0]
    line=lambda x: P[1]+m*(x-P[0])
    stem=(xx>=x0)&(xx<x1)
    if up: stem&=(yup>=0)&(yup<=line(xx))
    else:  stem&=(yup<=525)&(yup>=line(xx))
    img[stem]=1
    # the bowl may not poke past the stem's square end: clip it to the baseline / x-height near the stem
    near=(xx>=x0-T*1.2)&(xx<x1+T*1.2)
    if up: img[near&(yup<0)&(xx<x1) if l=='d' else near&(yup<0)&(xx>=x0)]=0
    else:  img[near&(yup>525)&(xx<x1) if l=='q' else near&(yup>525)&(xx>=x0)]=0
    return img,pad,noy
