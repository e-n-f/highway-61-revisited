import numpy as np,pickle,sys
import geomc4 as g4, geomc3 as g3
def err(l):
    img,pad,noy=g4.build(l,1.0,140); b,oy=g3.master(l)
    H,W=1100,950; A=np.zeros((H,W),np.uint8); B=np.zeros_like(A)
    def place(dst,src,oyy,dx):
        ys,xs=np.nonzero(src); Y=(ys-oyy+720).astype(int); X=(xs+dx).astype(int)
        ok=(Y>=0)&(Y<H)&(X>=0)&(X<W); dst[Y[ok],X[ok]]=1
    place(B,b,oy,60); place(A,img,noy,60-pad); return (A^B).sum()/B.sum()
tags={'r':['r'],'n':['n1'],'h':['h1'],'m':['m1','m2'],'t':['t'],'j':['j'],'a':['a'],'f':['f'],'c':['c']}
for l,tl in tags.items():
    for tg in tl:
        best=(err(l),0,0)
        for step in (16,6,2):
            cx,cy=best[1],best[2]
            for ddx in (-2,-1,0,1,2):
                for ddy in (-2,-1,0,1,2):
                    g4.OFF[tg]=(cx+ddx*step,cy+ddy*step); e=err(l)
                    if e<best[0]: best=(e,cx+ddx*step,cy+ddy*step)
        g4.OFF[tg]=(best[1],best[2]); print(l,tg,'offset',best[1],best[2],'err %.3f'%best[0],flush=True)
pickle.dump(g4.OFF,open('ring_off.pkl','wb'))
