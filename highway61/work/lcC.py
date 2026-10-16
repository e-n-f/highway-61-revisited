import pickle,numpy as np,cv2
A=pickle.load(open('A.pkl','rb'))
WR={**{c:0.81 for c in 'bdhijlmnpqrtuf'},**{c:0.82 for c in 'aceogs'},'k':0.79,'v':0.80,'w':0.87,'x':0.88,'y':0.78,'z':0.81}
T_OLD=115; R=5; T=T_OLD-2*R
def width(b): xs=np.nonzero(b.any(0))[0]; return xs[-1]-xs[0]+1
def ymap(yup):  # target y_up before erosion, given source y_up
    return np.where(yup<0, yup-R, np.where(yup<=525, yup*(525+2*R)/525-R, yup+R))
def transform(l):
    b,ox,oy,src=A[l]
    H,W=b.shape; Wd=width(b)
    wr=WR[l]
    fx=(wr*Wd-T)/(Wd-2*R-T) if l not in 'ijl' else 0.9
    fx=min(max(fx,0.5),1.0)
    pad=40
    # build inverse map for rows: target y_up range
    top=900; bot=-400; newH=top-bot; newoy=top
    tgt_y=newoy-np.arange(newH,dtype=np.float32)
    src_grid=np.linspace(-450,950,3000); mapped=ymap(src_grid)
    src_y=np.interp(tgt_y,mapped,src_grid)
    src_r=(oy-src_y).astype(np.float32)
    newW=int(W*fx)+2*pad
    src_c=((np.arange(newW,dtype=np.float32)-pad)/fx).astype(np.float32)
    mx,my=np.meshgrid(src_c,src_r)
    f=cv2.GaussianBlur(b.astype(np.float32),(0,0),1.5)
    out=cv2.remap(f,mx.astype(np.float32),my.astype(np.float32),cv2.INTER_LINEAR,borderValue=0)
    nb=(out>0.5).astype(np.uint8)
    # isotropic thinning (circular in unscaled space)
    se=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*max(1,round(R*fx))+1,2*R+1))
    nb=cv2.erode(nb,se)
    k=int(round(T*(1-fx)))
    if k>=1: nb=cv2.dilate(nb,np.ones((1,k+1),np.uint8))
    nb=(cv2.GaussianBlur(nb.astype(np.float32),(0,0),3)>0.5).astype(np.uint8)
    return nb,pad,newoy,'lcC',fx,Wd,width(nb)
A4={}
for l in 'abcdefghijklmnopqrstuvwxyz':
    nb,ox,oy,s,fx,w0,w1=transform(l); A4[l]=(nb,ox,oy,s)
    print(l,'fx %.2f'%fx,w0,w1,'%.3f'%(w1/w0))
pickle.dump(A4,open('A4.pkl','wb'))
