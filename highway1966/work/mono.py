import pickle,numpy as np,cv2
from skimage.morphology import skeletonize
import build
from build import L, smooth
A=pickle.load(open('A.pkl','rb'))
S=525/115.5; SPEC_ASC=152/115.5*525
def width(b): xs=np.nonzero(b.any(0))[0]; return xs[-1]-xs[0]+1
def source(l):
    d=L[l]; x,y,w,h=d['box']; m=10
    c=d['img'][y-m:y+h+m, x-m:x+w+m].astype(np.float32)
    big=cv2.resize(c,None,fx=S,fy=S,interpolation=cv2.INTER_CUBIC)
    big=cv2.GaussianBlur(big,(0,0),0.8*S)
    b=(big>128).astype(np.uint8)
    n,lab,st,_=cv2.connectedComponentsWithStats(b)
    keep=[i for i in range(1,n) if st[i][4]>2000]; b=np.isin(lab,keep).astype(np.uint8)
    oy=(d['base']-(y-m))*S
    return b,oy
def mono(l,fx,T=100,TH=90,DK=1.0,HO=62.0,band=25):
    rx,ry=T/2,TH/2
    # 1. corner/terminal envelope: lighter squeezed version via existing synth
    build.FX[l]=fx; build.STEM_T=T; build.ERODE_R=17; build.DK=DK
    E,ox,newoy=build.synth2(l)
    # 2. skeleton of source letter
    b,oy=source(l)
    dt=cv2.distanceTransform(b,cv2.DIST_L2,5)
    sk=skeletonize(b>0)
    sk=(sk&(dt>=40)).astype(np.uint8)
    K=np.ones((3,3),np.uint8)
    for _ in range(45):
        nb=cv2.filter2D(sk,-1,K,borderType=cv2.BORDER_CONSTANT)-sk
        ends=(sk==1)&(nb<=1)
        if not ends.any(): break
        sk[ends]=0
    sk=sk.astype(bool)
    rr,cc=np.nonzero(sk)
    yup=oy-rr
    lo,hi=HO,525-HO
    k_asc=(680-(525-ry))/(SPEC_ASC-hi)
    G=np.where(yup<0,(ry-lo)+yup*DK,
      np.where(yup<lo,yup-lo+ry,
      np.where(yup<=hi,ry+(yup-lo)*(525-2*ry)/(hi-lo),(525-ry)+(yup-hi)*k_asc)))
    X=30+cc*fx; Y=newoy-G
    D=np.zeros_like(E)
    for x,y,d in zip(X,Y,dt[rr,cc]):
        ax=min(rx,d*fx); ay=min(ry,d)
        if ax<2 or ay<2: continue
        cv2.ellipse(D,(int(round(x*4)),int(round(y*4))),(int(ax*4),int(ay*4)),0,0,360,1,-1,lineType=cv2.LINE_8,shift=2)
    D=smooth(D,4)
    near=cv2.dilate(D,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*band+1,2*band+1)))
    Ed=cv2.dilate(E,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(25,25)))
    F=E|(D&Ed)
    return smooth(F,5),ox,newoy
if __name__=='__main__':
    import sys
    out={}
    for l in sys.argv[1]:
        Wstd=L[l]['box'][2]*S; Wn=width(A[l][0]); Wt=(Wstd+Wn)/2
        fx=0.92 if l in 'ijl' else (Wt-100)/(Wstd-137)
        for it in range(3):
            b,ox,oy=mono(l,fx,DK=(1+230/186.0)/2)
            if l in 'ijl': break
            got=width(b)
            if abs(got-Wt)<4: break
            fx*= (Wt-100)/max(1,(got-100))
        out[l]=(b,ox,oy,'mono'); print(l,'target',int(Wt),'got',width(b))
    pickle.dump(out,open('Amono_test.pkl','wb'))
