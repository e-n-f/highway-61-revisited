import numpy as np,cv2,sys
from PIL import Image,ImageDraw,ImageFont
from hires import page,boxes,LAY
font='/mnt/user-data/outputs/Highway1966-E-Reference.ttf'
want=sys.argv[1]; out=sys.argv[2]
pages=[p for p in LAY if any(l in ''.join(LAY[p]) for l in want)]
tiles={}
for p in pages:
    lc=p in ('0031','0032','0033'); S=(525/457.75) if lc else (690/607.0)
    F=ImageFont.truetype(font,int(round(1000/S)))
    o,ang,k=page(p); bx=boxes(o,LAY[p])
    for l in ''.join(LAY[p]):
        if l not in want: continue
        x,y,w,h=bx[l]; m=40
        crop=o[y-m:y+h+m,x-m:x+w+m].copy()
        crop[:, :m-3]=0; crop[:, m+w+3:]=0; crop[:m-3,:]=0; crop[m+h+3:,:]=0
        ink=(crop>128).astype(np.uint8); H,W=ink.shape
        im=Image.new('L',(W+1600,H+1600),0); ImageDraw.Draw(im).text((400,200),l,font=F,fill=255)
        g=(np.array(im)>128).astype(np.uint8)
        ys,xs=np.nonzero(ink); ys2,xs2=np.nonzero(g)
        dx=(xs.min()+xs.max())/2-(xs2.min()+xs2.max())/2; dy=ys.max()-ys2.max() if l not in 'gjpqyQ' else ys.min()-ys2.min()
        g=cv2.warpAffine(g,np.float32([[1,0,dx],[0,1,dy]]),(W,H))
        T=np.full((H,W,3),255,np.uint8); T[ink>0]=(185,185,185)
        cs,_=cv2.findContours(g,cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE); cv2.drawContours(T,cs,-1,(0,0,220),2)
        tiles[l]=T
row=[tiles[l] for l in want if l in tiles]
Hm=max(t.shape[0] for t in row)
img=np.hstack([np.pad(t,((0,Hm-t.shape[0]),(0,20),(0,0)),constant_values=255) for t in row])
cv2.imwrite(out,img)
