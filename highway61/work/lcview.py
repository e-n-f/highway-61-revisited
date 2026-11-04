"""Letter outlines from the 1966 lower-case sheets in font units (baseline 0, x from the grid line left of the letter),
and a viewer with the grid drawn in."""
import numpy as np, sys, imgops
from lcscan import letters, SQ
DESC=set('gjpqy')
def units(ch):
    """contours in font units; origin: the grid line at the baseline and the grid line left of the letter"""
    r=letters()[ch]; cs=r['contours']; P=np.vstack(cs)
    y0=np.round(P[:,1].min()) if ch not in DESC else None
    if ch in DESC:   # baseline: 6 squares under the x-height line, found from the top of the bowl/stem
        top=P[:,1].max(); y0=np.round(top)-(8 if ch=='j' else 6)
        if ch=='j': y0=np.round(top)-8
    x0=np.floor(P[:,0].min()+0.2)
    return [np.c_[(c[:,0]-x0)*SQ,(c[:,1]-y0)*SQ] for c in cs]
def view(ch,extra=None,scale=1.0,out=None):
    cs=units(ch); P=np.vstack(cs); x1,y1=P.max(0)+60; x0,y0=P.min(0)-60
    W=int((x1-x0)*scale); H=int((y1-y0)*scale)
    im=np.full((H,W,3),255,np.uint8)
    T=lambda p: np.c_[(p[:,0]-x0)*scale,(y1-p[:,1])*scale]
    pts=[T(c) for c in cs]
    imgops.fill_polys(im,pts,(200,200,200))
    for k in range(int(np.floor(x0/SQ*4)),int(np.ceil(x1/SQ*4))+1):
        x=int((k*SQ/4-x0)*scale); imgops.line(im,(x,0),(x,H),(0,0,0) if k%4==0 else (255,220,220))
    for k in range(int(np.floor(y0/SQ*4)),int(np.ceil(y1/SQ*4))+1):
        y=int((y1-k*SQ/4)*scale); imgops.line(im,(0,y),(W,y),(0,0,0) if k%4==0 else (255,220,220))
    imgops.polylines(im,pts,(0,0,0))
    if extra is not None:
        for c,col in extra:                       # (colours given blue, green, red, as OpenCV took them)
            imgops.polylines(im,[T(np.asarray(c))],col[::-1],closed=False)
    if out: imgops.imwrite(out,im)
    return im
if __name__=='__main__':
    import tempfile; S=tempfile.gettempdir()+'/'
    for ch in sys.argv[1]:
        cs=units(ch); P=np.vstack(cs); print(ch,P.min(0).round(1),P.max(0).round(1))
        view(ch,scale=float(sys.argv[2]) if len(sys.argv)>2 else 1.0,out=S+f'lc_{ch}.png')
