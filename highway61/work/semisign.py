import pickle,numpy as np,cv2
from build import L
A=pickle.load(open('A.pkl','rb')); A5=pickle.load(open('A5.pkl','rb'))
S=525/115.5; T=100; R=7
def width(b): xs=np.nonzero(b.any(0))[0]; return xs[-1]-xs[0]+1
out=dict(A5)
for l in 'aehorswy':
    b,ox,oy,_=A[l]
    Wn=width(b); Wt=(L[l]['box'][2]*S+Wn)/2
    fx=(Wt-T)/(Wn-2*R-T)
    dk=0.9 if l=='y' else 1.0
    def ymap(yup): return np.where(yup<0, yup*dk-R, np.where(yup<=525, yup*(525+2*R)/525-R, yup+R))
    H,W=b.shape; pad=40
    top=900; bot=-400; newH=top-bot; newoy=top
    tgt=newoy-np.arange(newH,dtype=np.float32)
    g=np.linspace(-450,950,3000); src_y=np.interp(tgt,ymap(g),g)
    src_r=(oy-src_y).astype(np.float32)
    newW=int(W*fx)+2*pad
    src_c=((np.arange(newW,dtype=np.float32)-pad)/fx).astype(np.float32)
    mx,my=np.meshgrid(src_c,src_r)
    f=cv2.GaussianBlur(b.astype(np.float32),(0,0),1.5)
    nb=(cv2.remap(f,mx,my.astype(np.float32),cv2.INTER_LINEAR,borderValue=0)>0.5).astype(np.uint8)
    se=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*round(R*fx)+1,2*R+1))
    nb=cv2.erode(nb,se)
    k=int(round((115-2*R)*fx-T))
    if k>=1: nb=cv2.erode(nb,np.ones((1,k+1),np.uint8))
    nb=(cv2.GaussianBlur(nb.astype(np.float32),(0,0),3)>0.5).astype(np.uint8)
    out[l]=(nb,pad,newoy,'semisign')
    print(l,'narrow',Wn,'target %d got %d fx %.2f'%(Wt,width(nb),fx))
pickle.dump(out,open('A6.pkl','wb'))
