import pickle,numpy as np,cv2
def lighten(src,R):
    out={}
    for l,(b,ox,oy,s) in src.items():
        if not l.islower(): continue
        def ymap(yup): return np.where(yup<0, yup-R, np.where(yup<=525, yup*(525+2*R)/525-R, yup+R))
        H,W=b.shape; pad=30
        top=900; bot=-400; newH=top-bot; newoy=top
        tgt=newoy-np.arange(newH,dtype=np.float32)
        g=np.linspace(-450,950,3000); src_y=np.interp(tgt,ymap(g),g)
        src_r=(oy-src_y).astype(np.float32)
        src_c=(np.arange(W+2*pad,dtype=np.float32)-pad)
        mx,my=np.meshgrid(src_c,src_r)
        f=cv2.GaussianBlur(b.astype(np.float32),(0,0),1.0)
        o=cv2.remap(f,mx.astype(np.float32),my.astype(np.float32),cv2.INTER_LINEAR,borderValue=0)
        nb=(o>0.5).astype(np.uint8)
        nb=cv2.erode(nb,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*R+1,2*R+1)))
        nb=(cv2.GaussianBlur(nb.astype(np.float32),(0,0),1.2)>0.5).astype(np.uint8)
        out[l]=(nb,ox+pad,newoy,s)
    return out
if __name__=='__main__':
    A=pickle.load(open('A.pkl','rb')); A4=pickle.load(open('A4.pkl','rb'))
    pickle.dump(lighten(A,7),open('A_light.pkl','wb'))
    pickle.dump(lighten(A4,6),open('A4_light.pkl','wb'))
