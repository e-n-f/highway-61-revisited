import numpy as np,cv2,pickle
from skimage.morphology import skeletonize
import build
from build import L, smooth
S=525/115.5; SPEC_ASC=152/115.5*525; STEM_SPEC=137; HO=62.0
N8=[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
def source(l):
    d=L[l]; x,y,w,h=d['box']; m=10
    c=d['img'][y-m:y+h+m, x-m:x+w+m].astype(np.float32)
    big=cv2.resize(c,None,fx=S,fy=S,interpolation=cv2.INTER_CUBIC)
    big=cv2.GaussianBlur(big,(0,0),0.8*S)
    b=(big>128).astype(np.uint8)
    n,lab,st,_=cv2.connectedComponentsWithStats(b)
    b=np.isin(lab,[i for i in range(1,n) if st[i][4]>2000]).astype(np.uint8)
    b=smooth(b,4)
    return b,(d['base']-(y-m))*S
# ---------- skeleton graph ----------
def graph(b):
    dt=cv2.distanceTransform(b,cv2.DIST_L2,5)
    sk=skeletonize(b>0).astype(np.uint8)
    nb=cv2.filter2D(sk,-1,np.ones((3,3),np.uint8),borderType=cv2.BORDER_CONSTANT)-sk
    node=(sk==1)&(nb!=2)
    nn,nlab=cv2.connectedComponents(node.astype(np.uint8),connectivity=8)
    inner=(sk==1)&~node
    en,elab=cv2.connectedComponents(inner.astype(np.uint8),connectivity=8)
    H,W=b.shape
    def nbrs(p):
        r,c=p
        for dr,dc in N8:
            rr,cc=r+dr,c+dc
            if 0<=rr<H and 0<=cc<W: yield rr,cc
    paths=[]
    for e in range(1,en):
        pts=set(zip(*np.nonzero(elab==e)))
        ends=[p for p in pts if sum((q in pts) for q in nbrs(p))<=1]
        start=ends[0] if ends else next(iter(pts))
        order=[start]; seen={start}
        while True:
            nxt=[q for q in nbrs(order[-1]) if q in pts and q not in seen]
            if not nxt: break
            # prefer 4-neighbours
            nxt.sort(key=lambda q:abs(q[0]-order[-1][0])+abs(q[1]-order[-1][1]))
            order.append(nxt[0]); seen.add(nxt[0])
        # attach node clusters at both ends
        def node_at(p):
            for q in nbrs(p):
                if nlab[q]: return nlab[q]
            return 0
        a=node_at(order[0]); z=node_at(order[-1])
        closed=not ends
        paths.append(dict(pts=np.array(order,float),a=a,z=z,closed=closed))
    # node info
    nodes={}
    for k in range(1,nn):
        ys,xs=np.nonzero(nlab==k)
        deg=sum((p['a']==k)+(p['z']==k) for p in paths)
        nodes[k]=dict(c=np.array([ys.mean(),xs.mean()]),deg=deg)
    # add node centroids to path ends (so paths meet)
    for p in paths:
        if p['a']: p['pts']=np.vstack([nodes[p['a']]['c'],p['pts']])
        if p['z']: p['pts']=np.vstack([p['pts'],nodes[p['z']]['c']])
    return paths,nodes,dt
def plen(P): return np.sum(np.hypot(*np.diff(P,axis=0).T)) if len(P)>1 else 0
def prune(paths,nodes,dt):
    changed=True
    while changed:
        changed=False
        for k,nd in nodes.items():
            nd['deg']=sum((p['a']==k)+(p['z']==k) for p in paths)
        keep=[]
        for p in paths:
            ea=p['a'] and nodes[p['a']]['deg']==1; ez=p['z'] and nodes[p['z']]['deg']==1
            ja=p['a'] and nodes[p['a']]['deg']>=3; jz=p['z'] and nodes[p['z']]['deg']>=3
            if (ea and jz) or (ez and ja):
                j=p['z'] if ea else p['a']; c=nodes[j]['c'].astype(int)
                if plen(p['pts'])<1.25*dt[c[0],c[1]]: changed=True; continue
            keep.append(p)
        paths=keep
    for k,nd in nodes.items():
        nd['deg']=sum((p['a']==k)+(p['z']==k) for p in paths)
    return paths
def smooth_path(P,closed,w=7):
    if len(P)<5: return P
    Q=P.copy()
    for _ in range(3):
        if closed:
            R=sum(np.roll(Q,k,0) for k in range(-w,w+1))/(2*w+1)
        else:
            R=Q.copy()
            for i in range(1,len(Q)-1):
                k=min(w,i,len(Q)-1-i); R[i]=Q[i-k:i+k+1].mean(0)
        Q=R
    return Q
def trim(P,d,which):
    # remove length d from start ('a') or end ('z')
    if which=='z': P=P[::-1]
    acc=0; i=0
    while i<len(P)-2 and acc<d: acc+=np.hypot(*(P[i+1]-P[i])); i+=1
    P=P[i:]
    return P[::-1] if which=='z' else P
# ---------- one glyph at a width ----------
def glyph(l,f,stem,hratio=0.9,DK=1.1,band=12):
    rx=stem/2; ry=stem*hratio/2
    build.DOT=stem*1.03; build.FX[l]=f; build.STEM_T=stem-3
    build.ERODE_R=max(2,int(round((125-(2*ry-3))/2))); build.DK=DK
    E,ox,newoy=build.synth2(l)
    b,oy=source(l)
    paths,nodes,dt=graph(b)
    paths=prune(paths,nodes,dt)
    lo,hi=HO,525-HO
    k_asc=(680-(525-ry))/(SPEC_ASC-hi)
    def G(yup):
        return np.where(yup<0,(ry-lo)+yup*DK,np.where(yup<lo,yup-lo+ry,
               np.where(yup<=hi,ry+(yup-lo)*(525-2*ry)/(hi-lo),(525-ry)+(yup-hi)*k_asc)))
    D=np.zeros_like(E)
    for p in paths:
        P=smooth_path(p['pts'],p['closed'])
        if not p['closed']:
            if p['a'] and nodes[p['a']]['deg']==1: P=trim(P,rx,'a')
            if p['z'] and nodes[p['z']]['deg']==1: P=trim(P,rx,'z')
        # resample
        seg=np.hypot(*np.diff(P,axis=0).T); s=np.r_[0,np.cumsum(seg)]
        if s[-1]<1: continue
        t=np.arange(0,s[-1],2.0); R=np.interp(t,s,P[:,0]); C=np.interp(t,s,P[:,1])
        dd=dt[np.clip(R.astype(int),0,b.shape[0]-1),np.clip(C.astype(int),0,b.shape[1]-1)]
        X=30+C*f; Y=newoy-G(oy-R)
        for x,y,d0 in zip(X,Y,dd):
            ax=min(rx,d0*max(f,1)*1.0); ay=min(ry,d0)
            if ax<3 or ay<3: continue
            cv2.ellipse(D,(int(round(x*4)),int(round(y*4))),(int(ax*4),int(ay*4)),0,0,360,1,-1,shift=2)
    Ed=cv2.dilate(E,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*band+1,2*band+1)))
    F=E|(D&Ed)
    return smooth(F,5),ox,newoy
WIDTHS={'wide':1.15,'standard':1.0,'semi':0.89,'narrow':0.78,'condensed':0.66}
def stem_for(f): return 113*f**0.75
if __name__=='__main__':
    import sys
    letters=sys.argv[1]
    out={}
    for name,f in WIDTHS.items():
        for l in letters:
            out[(name,l)]=glyph(l,f,stem_for(f))
    pickle.dump(out,open('fam_test.pkl','wb'))
    print('ok')
