import numpy as np,cv2,pickle,sys
from family import graph,prune,smooth_path,trim
from build import smooth
M=pickle.load(open('master_mk.pkl','rb'))     # 1966 lower case with marked straight edges
HO=62.0; STEM_M=140.0; BAR_M=125.0
def envelope(b,oy,f,R,T):
    # pre-stretch vertically so isotropic erosion R keeps baseline/x-height/ascender/descender
    def ymap(y): return np.where(y<0,y-R,np.where(y<=525,y*(525+2*R)/525-R,y+R))
    H,W=b.shape; pad=40
    top=950; bot=-450; newH=top-bot; newoy=top
    tgt=newoy-np.arange(newH,dtype=np.float32)
    g=np.linspace(-500,1000,4000); src_y=np.interp(tgt,ymap(g),g)
    src_r=(oy-src_y).astype(np.float32)
    newW=int(W*f)+2*pad
    src_c=((np.arange(newW,dtype=np.float32)-pad)/f).astype(np.float32)
    mx,my=np.meshgrid(src_c,src_r)
    o=cv2.remap(b.astype(np.float32),mx.astype(np.float32),my.astype(np.float32),cv2.INTER_LINEAR,borderValue=0)
    e=(o>0.5).astype(np.uint8)
    se=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*max(1,int(round(R*f)))+1,2*int(round(R))+1))
    e=cv2.erode(e,se)
    k=int(round(T-(STEM_M-2*R)*f))
    if k>=1: e=cv2.dilate(e,np.ones((1,k+1),np.uint8))
    elif k<=-1: e=cv2.erode(e,np.ones((1,-k+1),np.uint8))
    return e,pad,newoy
def glyph(l,f,T,hratio=0.9,band=4):
    b,ox,oy,_=M[l]
    rx=T/2; ry=T*hratio/2
    R=max(1.0,(BAR_M-(2*ry-3))/2)
    E,pad,newoy=envelope(b,oy,f,R,T-3)
    paths,nodes,dt=graph(b); paths=prune(paths,nodes,dt)
    lo,hi=HO,525-HO
    def G(y):
        return np.where(y<lo,y-lo+ry,np.where(y<=hi,ry+(y-lo)*(525-2*ry)/(hi-lo),y-hi+525-ry))
    D=np.zeros_like(E)
    for p in paths:
        P=smooth_path(p['pts'],p['closed'])
        if not p['closed']:
            if p['a'] and nodes[p['a']]['deg']==1: P=trim(P,rx,'a')
            if p['z'] and nodes[p['z']]['deg']==1: P=trim(P,rx,'z')
        seg=np.hypot(*np.diff(P,axis=0).T); s=np.r_[0,np.cumsum(seg)]
        if s[-1]<1: continue
        t=np.arange(0,s[-1],2.0); Rr=np.interp(t,s,P[:,0]); C=np.interp(t,s,P[:,1])
        dd=dt[np.clip(Rr.astype(int),0,b.shape[0]-1),np.clip(C.astype(int),0,b.shape[1]-1)]
        X=pad+C*f; Y=newoy-G(oy-Rr)
        for x,y,d0 in zip(X,Y,dd):
            ax=min(rx,d0*max(f,1)); ay=min(ry,d0)
            if ax<3 or ay<3: continue
            cv2.ellipse(D,(int(round(x*4)),int(round(y*4))),(int(ax*4),int(ay*4)),0,0,360,1,-1,shift=2)
    D=smooth(D,3)
    Ed=cv2.dilate(E,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*band+1,2*band+1)))
    F=E|(D&Ed)
    F=smooth(F,1.5)
    return F,pad,newoy
SER={'B':(0.654,91),'C':(0.78,102),'D':(0.89,110),'F':(1.10,132)}
if __name__=='__main__':
    import square
    s=sys.argv[1]; f,T=SER[s]
    out={}
    for l in 'abcdefghijklmnopqrstuvwxyz':
        F,pad,oy=glyph(l,f,T); out[l]=(F,pad,oy,'fam2-'+s)
    try: out=square.apply(out)
    except Exception as e: print('square skipped',e)
    pickle.dump(out,open('lcE_%s.pkl'%s,'wb'))
    from src import stem_width
    for l,y in [('n',250),('b',250),('o',250),('l',300)]:
        b,ox,oy,_=out[l]; print(s,l,stem_width(b,oy,y/525))
