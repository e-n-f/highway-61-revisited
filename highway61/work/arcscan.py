import numpy as np,cv2,pickle
from scipy.optimize import least_squares
import dim1961 as d
CE=pickle.load(open('caps_E.pkl','rb'))
def fitc(P):
    x,y=P[:,0],P[:,1]
    A_=np.column_stack([x,y,np.ones_like(x)]); b=x*x+y*y
    c=np.linalg.lstsq(A_,b,rcond=None)[0]; cx,cy=c[0]/2,c[1]/2; r=np.sqrt(c[2]+cx*cx+cy*cy)
    return cx,cy,r,np.abs(np.hypot(x-cx,y-cy)-r).max()
def scan(ch):
    b,ox,oy,_=CE[ch]; x0=np.nonzero(b.any(0))[0][0]
    cs,_=cv2.findContours(b.copy(),cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)
    out=[]
    for c in cs:
        P=c[:,0,:].astype(float); P=np.column_stack([P[:,0]-x0,oy-P[:,1]])
        n=len(P)
        if n<100: continue
        # greedy growth of arcs: extend while circle fit error < 2.5 units
        i=0
        while i<n:
            j=i+30
            if j>=n: break
            best=None
            while j<n:
                cx,cy,r,e=fitc(P[i:j])
                if e>2.5: break
                best=(i,j,cx,cy,r); j+=10
            if best and best[1]-best[0]>=60:
                i0,j0,cx,cy,r=best
                if r<2000: out.append((round(cx),round(cy),round(r),j0-i0, tuple(np.round(P[i0])), tuple(np.round(P[j0-1]))))
                i=j0
            else: i+=10
    return out
names={k:v for k,v in d.DIM.items() if k in 'GIKLMNPQR'}
for ch in '23':
    print('==',ch)
    for cx,cy,r,n,p0,p1 in scan(ch):
        near=min(list(names.items())+[(k+'-A',v-d.A) for k,v in names.items()]+[(k+'+A',v+d.A) for k,v in names.items()],key=lambda kv:abs(kv[1]-r))
        print('centre (%4d,%4d) r %4d  pts %4d  from %s to %s   nearest %s=%.0f'%(cx,cy,r,n,p0,p1,near[0],near[1]))
