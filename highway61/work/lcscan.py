"""Letter outlines from the 1966 lower-case sheets, in grid squares (lcgrid) and font units.

Ink is opened with a small disk to drop the grid lines that touch the letters, then contoured at sub-pixel precision
on the grey image.  Coordinates: gx to the right, gy upwards, both in squares from the page's first grid lines; font
units = squares x SQ, with the baseline and a reference vertical line chosen per letter."""
import numpy as np, cv2, pickle, os
from skimage import measure
from lcgrid import Page, SQ
PAGES={'p33':'abcdefghi','p34':'jklmnopqr','p35':'stuvwxyz'}
_cache={}
def page(n):
    if n not in _cache: _cache[n]=Page(n)
    return _cache[n]
def blobs(n):
    P=page(n); im=P.im
    d=(im<128).astype(np.uint8)
    k=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(9,9))
    L=cv2.morphologyEx(d,cv2.MORPH_OPEN,k)
    nl,lab,st,_=cv2.connectedComponentsWithStats(L)
    out=[]
    for i in range(1,nl):
        x,y,w,h,a=st[i]
        if a<3000 or y<250 or y+h>2450 or w<40 or h<40: continue
        out.append((x,y,w,h,lab==i))
    return out
def contours(n,mask):
    """sub-pixel contours of one letter: grey image, with everything away from the opened mask set to white"""
    P=page(n); g=P.im.astype(float)
    near=cv2.dilate(mask.astype(np.uint8),np.ones((5,5),np.uint8))>0
    inner=cv2.erode(mask.astype(np.uint8),np.ones((5,5),np.uint8))>0
    g=np.where(near,g,255.0); g=np.where(inner,0.0,g)
    g=cv2.GaussianBlur(g,(0,0),0.8)
    cs=measure.find_contours(g,128)
    out=[]
    for c in cs:
        if len(c)<40: continue
        y,x=c[:,0],c[:,1]; gx,gy=P.to_grid(x,y); out.append(np.c_[gx,-gy])
    return out
def letters():
    fn='lc1966/contours.pkl'
    if os.path.exists(fn): return pickle.load(open(fn,'rb'))
    res={}
    for n,chars in PAGES.items():
        bl=blobs(n)
        # group blobs into letters (i and j dots): sort by row band then x
        bl.sort(key=lambda b:b[3]*b[2],reverse=True)
        groups=[]
        for b in bl:
            for gp in groups:
                x,y,w,h,_=gp[0]
                if b[2]<130 and b[3]<130 and abs((b[0]+b[2]/2)-(x+w/2))<w/2+10 and (b[1]+b[3]<y+5 or y+h<b[1]+5) and min(abs(b[1]+b[3]-y),abs(y+h-b[1]))<120:
                    gp.append(b); break
            else: groups.append([b])
        groups.sort(key=lambda gp:(int(max(b[1]+b[3] for b in gp)//600),min(b[0] for b in gp)))
        print(n,len(groups),[ (min(b[0] for b in gp),min(b[1] for b in gp)) for gp in groups])
        for ch,gp in zip(chars,groups):
            cs=[]
            for b in gp: cs+=contours(n,b[4])
            res[ch]=dict(page=n,contours=cs)
    pickle.dump(res,open(fn,'wb')); return res
if __name__=='__main__':
    R=letters()
    for ch,r in R.items():
        allp=np.vstack(r['contours']); print(ch,r['page'],len(r['contours']),allp.min(0).round(2),allp.max(0).round(2))
