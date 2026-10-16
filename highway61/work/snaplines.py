import cv2,numpy as np,pickle
def unit(v): return v/np.hypot(*v)
def refine(b,A,B,band=14,minrun=8,maxang=6):
    """find boundary run along marked line A->B; return fitted (c,u,t0,t1,side) or None"""
    u=unit(B-A); n=np.array([-u[1],u[0]]); L=np.hypot(*(B-A))
    cs,_=cv2.findContours(b.copy(),cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)
    best=None; runs_out=[]
    for c in cs:
        P=c[:,0,:].astype(float)
        d=(P-A)@n; t=(P-A)@u
        ok=(np.abs(d)<=band)&(t>-40)&(t<L+40)
        if not ok.any(): continue
        # cyclic runs
        idx=np.nonzero(ok)[0]; N=len(P)
        if ok.all(): runs=[np.arange(N)]
        else:
            start=np.nonzero(~ok)[0][0]; order=np.roll(np.arange(N),-start); f=ok[order]
            runs=[]; cur=[]
            for i,v in zip(order,f):
                if v: cur.append(i)
                elif cur: runs.append(np.array(cur)); cur=[]
            if cur: runs.append(np.array(cur))
        for r in runs:
            if len(r)<minrun: continue
            Q=P[r]; m=Q.mean(0); _,_,vt=np.linalg.svd(Q-m); uu=vt[0]
            if uu@u<0: uu=-uu
            ang=np.degrees(np.arccos(min(1,abs(uu@u))))
            if ang>maxang: continue
            for _ in range(2):
                nn_=np.array([-uu[1],uu[0]]); res=np.abs((Q-m)@nn_)
                good=res<=5.0
                # longest contiguous good stretch
                bi=None; cur=0; st=0
                for i,g in enumerate(good):
                    if g:
                        if cur==0: st=i
                        cur+=1
                        if bi is None or cur>bi[1]-bi[0]: bi=(st,i+1)
                    else: cur=0
                if bi is None or bi[1]-bi[0]<minrun: break
                Q=Q[bi[0]:bi[1]]; m=Q.mean(0); _,_,vt=np.linalg.svd(Q-m); uu=vt[0]
                if uu@u<0: uu=-uu
            if len(Q)<minrun: continue
            tt=(Q-m)@uu
            runs_out.append((m,uu,tt.min(),tt.max(),len(Q)))
    res=[]
    for (m,uu,t0,t1,_) in runs_out:
        res.append(_finish(b,m,uu,t0,t1))
    return res
def _finish(b,m,uu,t0,t1):
    nn=np.array([-uu[1],uu[0]])
    # ink side
    H,W=b.shape; mid=m+uu*(t0+t1)/2
    def ink(p):
        x,y=int(round(p[0])),int(round(p[1]))
        return 0<=x<W and 0<=y<H and b[y,x]>0
    s=sum(ink(m+uu*t+nn*8) for t in np.linspace(t0,t1,9))-sum(ink(m+uu*t-nn*8) for t in np.linspace(t0,t1,9))
    side=1 if s>0 else -1
    return dict(c=m,u=uu,n=nn,t0=t0,t1=t1,side=side)
def intersect(a,b):
    M=np.array([a['u'],-b['u']]).T
    if abs(np.linalg.det(M))<1e-6: return None
    s,_=np.linalg.solve(M,b['c']-a['c']); return s
def apply(b,lines,band=5,ext=25):
    b=b.copy(); H,W=b.shape
    for i,L in enumerate(lines):
        t0,t1=L['t0'],L['t1']
        for j,K in enumerate(lines):
            if i==j: continue
            s=intersect(L,K)
            if s is None: continue
            if t1<s<=t1+ext: t1=s
            if t0-ext<=s<t0: t0=s
        L['T0'],L['T1']=t0,t1
        L['c0']=(t0!=L['t0']); L['c1']=(t1!=L['t1'])   # ends fixed at a corner
    yy,xx=np.mgrid[0:H,0:W]; P=np.stack([xx,yy],-1).astype(float)
    for L in lines:
        d=(P-L['c'])@L['n']; t=(P-L['c'])@L['u']
        zone=(np.abs(d)<=band)&(t>=L['T0'])&(t<=L['T1'])
        inkside=(d*L['side']>0)
        b[zone&inkside]=1; b[zone&~inkside]=0
    return b
def constrain(b,marks):
    lines=[]
    for A,B in marks:
        lines+=refine(b,np.asarray(A,float),np.asarray(B,float))
    nb=apply(b,lines)
    se=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(7,7))
    nb=cv2.morphologyEx(cv2.morphologyEx(nb,cv2.MORPH_OPEN,se),cv2.MORPH_CLOSE,se)
    nb=apply(nb,lines)   # re-assert the straight edges/corners after cleanup
    return nb,lines
if __name__=='__main__':
    M=pickle.load(open('hr_lc_sq.pkl','rb')); R=pickle.load(open('letterlines.pkl','rb'))
    out={}; LL={}
    for l,(b,ox,oy,t) in M.items():
        nb,lines=constrain(b,R.get(l,[])); out[l]=(nb,ox,oy,t); LL[l]=lines
        print(l,len(R.get(l,[])),'->',len(lines),'changed',int((nb!=b).sum()))
    pickle.dump(out,open('master_mk.pkl','wb')); pickle.dump(LL,open('master_lines.pkl','wb'))
