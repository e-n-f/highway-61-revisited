"""Grid rectification for the 1966 lower-case sheets (lc1966/p33-35.jpg).

Every printed grid line is fitted individually (x = a + b*y for verticals, y = a + b*x for horizontals); a scan point is
mapped to grid coordinates by interpolating between the lines either side of it, so the letters are measured against
the squares they were drawn on, whatever the printing or scanning did to the pitch.  One square = 1 inch of an 8-inch
letter = 690/8 = 86.25 font units."""
import numpy as np, cv2
SQ=690/8
def _peaks(s,thr,sep=25):
    idx=[i for i in range(1,len(s)-1) if s[i]>thr and s[i]>=s[i-1] and s[i]>=s[i+1]]
    out=[]
    for i in idx:
        if out and i-out[-1]<sep:
            if s[i]>s[out[-1]]: out[-1]=i
        else: out.append(i)
    return out
def _fill(pos,pitch=59.0):
    pos=list(pos); out=[pos[0]]
    for p in pos[1:]:
        n=int(round((p-out[-1])/pitch))
        if n<1: continue
        for k in range(1,n): out.append(out[-1]+(p-out[-1])/(n-k+1))
        out.append(p)
    return out
class Page:
    def __init__(self,name,x0=150,x1=1960,y0=150,y1=2500):
        im=cv2.imread(f'lc1966/{name}.jpg',0); self.im=im
        d=(im<150).astype(np.uint8)
        blob=cv2.dilate(cv2.morphologyEx(d,cv2.MORPH_OPEN,np.ones((9,9),np.uint8)),np.ones((13,13),np.uint8))
        thin=d&(1-blob); self.thin=thin; self.ink=(im<128)&(blob>0)
        ys_,xs_=np.nonzero(thin)
        self.V=self._lines(xs_,ys_,x0,x1,y0,y1,len(im[0]))     # x = a + b*y
        self.H=self._lines(ys_,xs_,y0,y1,x0,x1,len(im))        # y = a + b*x
    def _lines(self,c,t,c0,c1,t0,t1,n):
        """lines c = a + b*t: best common shear, lattice of peaks of the sheared projection, then each line refitted"""
        k=(t>=t0)&(t<=t1); c,t=c[k].astype(float),t[k].astype(float)
        best=None
        for b in np.arange(-0.02,0.02,0.0005):
            h=np.bincount(np.clip((c-b*t).astype(int),0,n-1),minlength=n).astype(float); h[:c0]=0; h[c1:]=0
            sc=(h**2).sum()
            if best is None or sc>best[0]: best=(sc,b,h)
        _,b,h=best
        pk=_peaks(h,np.percentile(h[h>0],90)*0.3,sep=25); pk=self._lattice(pk,h)
        out=[]
        for p in pk:
            a=p
            for w in (8,4,3):
                sel=np.abs(c-(a+b*t))<w
                if sel.sum()<50: break
                A=np.vstack([np.ones(sel.sum()),t[sel]]).T; a,b2=np.linalg.lstsq(A,c[sel],rcond=None)[0]
                a=a+(b2-b)*0 ; aa,bb=a,b2
            out.append((a,b2) if sel.sum()>=50 else (p,b))
        return out
    @staticmethod
    def _lattice(p,s):
        """robust regular lattice through the peaks, then each line re-found at the projection's maximum nearby"""
        p=np.array(p,float); best=None
        for i in range(len(p)):
            for j in range(i+1,len(p)):
                n=round((p[j]-p[i])/59.0)
                if n<5: continue
                pitch=(p[j]-p[i])/n; k=np.round((p-p[i])/pitch); r=np.abs(p-(p[i]+k*pitch)); sc=(r<4).sum()
                if best is None or sc>best[0]: best=(sc,p[i],pitch)
        _,p0,pitch=best; k=np.round((p-p0)/pitch); inl=np.abs(p-(p0+k*pitch))<4
        pitch,p0=np.polyfit(k[inl],p[inl],1)
        out=[]
        for kk in range(int(k.min())-1,int(k.max())+2):
            e=p0+kk*pitch
            if e<5 or e>len(s)-5: continue
            lo=int(e-8); seg=s[lo:int(e+9)]
            if seg.max()<=0: out.append(e); continue
            m=lo+int(np.argmax(seg)); out.append(m if seg.max()>0.25*np.percentile(s[s>0],90) else e)
        return out
    @staticmethod
    def _regular(p):
        p=np.array(p,float); d=np.diff(p)
        keep=[p[0]]
        for v in p[1:]:
            if v-keep[-1]>40: keep.append(v)
        return _fill(keep)
    @staticmethod
    def _fit(t,c,w=9):
        """t: array indexed [line coord, along]; returns (a,b) with line coord = a + b*along"""
        c=int(round(c)); band=t[max(c-w,0):c+w+1]; ii,jj=np.nonzero(band)
        if len(ii)<50: return (float(c),0.0)
        ii=ii+max(c-w,0)
        for _ in range(3):
            A=np.vstack([np.ones_like(jj),jj]).T.astype(float); a,b=np.linalg.lstsq(A,ii.astype(float),rcond=None)[0]
            r=ii-(a+b*jj); k=np.abs(r)<2.5
            if k.all(): break
            ii,jj=ii[k],jj[k]
        return (a,b)
    def to_grid(self,x,y):
        """scan px -> (gx, gy) in squares; gx counts vertical lines from the first, gy horizontal lines from the first
        (downwards)"""
        x=np.asarray(x,float); y=np.asarray(y,float)
        vx=np.array([a+b*y for a,b in self.V])       # positions of each vertical line at height y
        hy=np.array([a+b*x for a,b in self.H])
        return _interp(x,vx),_interp(y,hy)
def _interp(v,lines):
    """fractional index of v among monotone line positions (per point), linear extrapolation outside"""
    lines=np.atleast_2d(lines.T) if lines.ndim>1 else lines[None]
    v=np.atleast_1d(v); out=np.empty(v.shape)
    for k,(vv,L) in enumerate(zip(v.ravel(),lines.reshape(-1,lines.shape[-1]))):
        i=np.clip(np.searchsorted(L,vv)-1,0,len(L)-2); out.ravel()[k]=i+(vv-L[i])/(L[i+1]-L[i])
    return out
if __name__=='__main__':
    for n in ['p33','p34','p35']:
        P=Page(n); print(n,len(P.V),len(P.H),np.diff([a for a,b in P.V]).round(1).tolist())
        print('   H',np.diff([a for a,b in P.H]).round(1).tolist())
        print('   tilt V',np.round([b for a,b in P.V],4).tolist()[:6],' H',np.round([b for a,b in P.H],4).tolist()[:6])
