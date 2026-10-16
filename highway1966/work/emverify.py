"""Overlay the widened Series E capitals and figures (Highway61Revisited-EM) on the 1983 E(M) tracings (caps_EM.pkl):
python3 emverify.py -> outputs/verify-EM-caps.png (grey = 1983, red = widened 1961 E; numbers: p95 / median distance)"""
import pickle, numpy as np, cv2
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen
from scipy.spatial import cKDTree
D=pickle.load(open('caps_EM.pkl','rb'))
F=TTFont('../outputs/Highway61Revisited-EM.otf'); gs=F.getGlyphSet(); cm=F.getBestCmap()
def rings(name):
    rp=RecordingPen(); gs[name].draw(rp); out=[]; cur=[]
    for op,a in rp.value:
        if op=='moveTo': cur=[np.array(a[0])]
        elif op=='lineTo': cur.append(np.array(a[0]))
        elif op=='curveTo':
            p0=cur[-1]; p1,p2,p3=map(np.array,a)
            for t in np.linspace(0,1,17)[1:]: cur.append((1-t)**3*p0+3*(1-t)**2*t*p1+3*(1-t)*t*t*p2+t**3*p3)
        elif op=='closePath': out.append(np.array(cur))
    return out
tiles=[]; stats=[]
for ch in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789':
    b,ox,oy,t=D[ch]; ys,xs=np.nonzero(b)
    R=rings(cm[ord(ch)]); P=np.vstack(R)
    # align horizontally on the bounding-box centre (the tracings have their own side bearings)
    dx=(xs.min()+xs.max())/2-(P[:,0].min()+P[:,0].max())/2
    cs,_=cv2.findContours((b>0).astype(np.uint8),cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)
    S=np.vstack([c[:,0,:] for c in cs if len(c)>20]).astype(float); S=np.c_[S[:,0],oy-S[:,1]]
    def dense(r,step=2.0):
        Q=np.vstack([r,r[:1]]); out=[]
        for p,q in zip(Q[:-1],Q[1:]):
            n=max(1,int(np.hypot(*(q-p))/step)); out.append(p+(q-p)*np.linspace(0,1,n,endpoint=False)[:,None])
        return np.vstack(out)
    M=np.vstack([dense(np.c_[r[:,0]+dx,r[:,1]]) for r in R])
    a,_=cKDTree(M).query(S); c,_=cKDTree(S).query(M); dd=np.r_[a,c]; p95=np.percentile(dd,95); stats.append(p95)
    s=0.35; H=int(820*s); W=int(780*s); im=np.full((H+40,W,3),255,np.uint8)
    T=lambda Q: np.round(np.c_[(Q[:,0]+40)*s,(760-Q[:,1])*s]*16).astype(np.int32)
    sm=cv2.resize((b>0).astype(np.uint8)*255,None,fx=s,fy=s,interpolation=cv2.INTER_AREA)
    y0=int(round((760-oy)*s)); x0=int(round(40*s))
    hh,ww=sm.shape; sub=im[y0:y0+hh,x0:x0+ww]; hh=min(hh,sub.shape[0]); ww=min(ww,sub.shape[1])
    im[y0:y0+hh,x0:x0+ww]=np.minimum(im[y0:y0+hh,x0:x0+ww],(255-sm[:hh,:ww]*0.45)[...,None].astype(np.uint8))
    for r in R: cv2.polylines(im,[T(np.c_[r[:,0]+dx,r[:,1]])],True,(0,0,230),1,cv2.LINE_AA,shift=4)
    cv2.putText(im,'%s  %.0f / %.1f'%(ch,p95,np.median(dd)),(5,H+28),0,0.6,(0,0,0),1,cv2.LINE_AA)
    tiles.append(im)
rows=[np.hstack(tiles[i:i+9]) for i in range(0,36,9)]
cv2.imwrite('../outputs/verify-EM-caps.png',np.vstack(rows))
print('mean p95 %.1f, median of p95 %.1f'%(np.mean(stats),np.median(stats)))
