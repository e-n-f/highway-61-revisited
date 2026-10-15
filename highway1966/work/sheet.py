import numpy as np,cv2
def sheet(items,fn,scale=0.5,gap=40):
    base=850; H=1250; tot=sum(b.shape[1] for b,_,_ in items)+gap*len(items)+40
    S=np.zeros((H,tot),np.uint8); x=20
    for b,ox,oy in items:
        h,w=b.shape; top=int(round(base-oy))
        for r in range(h):
            yy=top+r
            if 0<=yy<H: S[yy,x:x+w]|=(b[r]*255).astype(np.uint8)
        x+=w+gap
    S=255-S
    for yy in (base,base-525,base-680,base+230): cv2.line(S,(0,yy),(tot,yy),160,2)
    cv2.imwrite(fn,cv2.resize(S,None,fx=scale,fy=scale,interpolation=cv2.INTER_AREA))
