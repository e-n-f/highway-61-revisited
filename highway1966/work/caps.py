import numpy as np,cv2
from src import *
CAP=1000
def spec_cap(l):
    o=np.load('p12.npy')
    box=dict(B=[645,250,207,303],D=[304,667,205,304],F=[988,667,185,303])[l]
    b,ox,oy,S=raster(o,box,box[1]+box[3],box[3]*CAP/525.0*525/CAP if False else box[3]*525/CAP,thr=128)
    return b,oy
def sign_cap(src,idx):
    r,res,t=sign_letters(src)
    g=res[idx[0]][idx[1]]
    b,ox,oy,S=raster(r,g,g[1]+g[3],g[3]*525/CAP)
    return b,oy
def metrics(b,oy):
    ys,xs=np.nonzero(b); x0=xs.min(); H=oy-ys.min()
    W=xs.max()-x0+1
    def hrun(fr):
        row=b[int(oy-H*fr)]; xs_=np.nonzero(row)[0]; 
        runs=[];s=None
        for x in range(row.size):
            if row[x] and s is None: s=x
            if not row[x] and s is not None: runs.append((s-x0,x-s)); s=None
        return runs
    def vrun(fx):
        col=b[:,int(x0+W*fx)]; runs=[];s=None
        for y in range(col.size):
            if col[y] and s is None: s=y
            if not col[y] and s is not None: runs.append((round(oy-y),y-s)); s=None
        return runs
    return dict(W=W,H=H,stem=hrun(0.5)[0][1] if hrun(0.5) else None,
                r25=hrun(0.25),r75=hrun(0.75),v55=vrun(0.55),v80=vrun(0.8))
