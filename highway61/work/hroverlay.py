import pickle,numpy as np,cv2,sys
from PIL import Image,ImageDraw,ImageFont
from hires import page,boxes,LAY
font='/mnt/user-data/outputs/Highway1966-E-Reference.ttf'
def overlay(pages,target,refpx,out):
    S=target/refpx; size=int(round(1000/S))
    F=ImageFont.truetype(font,size)
    tiles=[];errs=[]
    for p in pages:
        o,ang,k=page(p); bx=boxes(o,LAY[p])
        for r in LAY[p]:
            for l in r:
                x,y,w,h=bx[l]; m=40
                crop=o[y-m:y+h+m,x-m:x+w+m].copy()
                crop[:, :m-3]=0; crop[:, m+w+3:]=0; crop[:m-3,:]=0; crop[m+h+3:,:]=0
                ink=(crop>128).astype(np.uint8); H,W=ink.shape
                im=Image.new('L',(W+1600,H+1600),0); ImageDraw.Draw(im).text((400,200),l,font=F,fill=255)
                g=(np.array(im)>128).astype(np.uint8)
                ys,xs=np.nonzero(ink); ys2,xs2=np.nonzero(g)
                dx=(xs.min()+xs.max())/2-(xs2.min()+xs2.max())/2; dy=ys.max()-ys2.max() if l not in 'gjpqyQ' else (ys.min()-ys2.min())
                g=cv2.warpAffine(g,np.float32([[1,0,dx],[0,1,dy]]),(W,H))
                e=(ink^g).sum()/ink.sum(); errs.append((l,round(float(e),3)))
                T=np.full((H,W,3),255,np.uint8); T[ink>0]=(185,185,185)
                cs,_=cv2.findContours(g,cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE); cv2.drawContours(T,cs,-1,(0,0,220),2)
                ys,xs=np.nonzero(ink|g); tiles.append(T[max(0,ys.min()-15):ys.max()+15,max(0,xs.min()-15):xs.max()+15])
    Hm=max(t.shape[0] for t in tiles); rows=[]
    per=9
    for i in range(0,len(tiles),per):
        rows.append(np.hstack([np.pad(t,((0,Hm-t.shape[0]),(0,15),(0,0)),constant_values=255) for t in tiles[i:i+per]]))
    Wm=max(r.shape[1] for r in rows)
    img=np.vstack([np.pad(r,((0,10),(0,Wm-r.shape[1]),(0,0)),constant_values=255) for r in rows])
    cv2.imwrite(out,cv2.resize(img,None,fx=1600/img.shape[1],fy=1600/img.shape[1],interpolation=cv2.INTER_AREA))
    print(errs)
if __name__=="__main__" and sys.argv[1]=="lc": overlay(['0031','0032','0033'],525,457.75,'/mnt/user-data/outputs/hires-lowercase-vs-scan.png')
elif __name__=="__main__": overlay(['0016','0017','0018','0019'],690,607.0,'/mnt/user-data/outputs/hires-capitals-vs-scan.png')
