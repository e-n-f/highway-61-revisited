import pickle,numpy as np,sys
def edges(b,side):
    H,W=b.shape; e=np.full(H,np.nan)
    for r in range(H):
        xs=np.nonzero(b[r])[0]
        if len(xs): e[r]=xs[0] if side=='L' else xs[-1]
    return e
def scan(b,side,win=25,thr=3.0):
    e=edges(b,side); H=len(e); hits=[]
    for r in range(win,H-win):
        seg=e[r-win:r+win+1]
        if np.isnan(seg).any(): continue
        d=np.abs(np.diff(seg))
        if d.max()>6: continue          # edge jumps: bars/corners, skip
        t=np.arange(-win,win+1); A=np.vstack([t**2,t,np.ones_like(t)]).T
        coef,*_=np.linalg.lstsq(A,seg,rcond=None)
        res=seg[win]-coef[2]
        if abs(res)>thr: hits.append((r,res))
    # cluster
    out=[]
    for r,res in hits:
        if out and r-out[-1][1]<=3: out[-1][1]=r; out[-1][2]=max(out[-1][2],abs(res))
        else: out.append([r,r,abs(res)])
    return out
def fix(b,side,r0,r1,pad=12):
    b=b.copy(); e=edges(b,side)
    a=r0-pad; z=r1+pad
    for r in range(a+1,z):
        x=e[a]+(e[z]-e[a])*(r-a)/(z-a); xi=int(round(x))
        if side=='L':
            cur=int(e[r]); 
            if xi<cur: b[r,xi:cur]=1
            else: b[r,cur:xi]=0
        else:
            cur=int(e[r])
            if xi>cur: b[r,cur+1:xi+1]=1
            else: b[r,xi+1:cur+1]=0
    return b
if __name__=='__main__':
    for f in sys.argv[1:]:
        d=pickle.load(open(f,'rb'))
        for k in sorted(d):
            b,ox,oy,_=d[k]
            for side in 'LR':
                h=scan(b,side)
                if h: print(f[:6],k,side,[(int(oy-r1),int(oy-r0),round(a,1)) for r0,r1,a in h])
