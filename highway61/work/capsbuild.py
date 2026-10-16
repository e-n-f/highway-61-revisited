import cv2,numpy as np,pickle
from src import rot
CAPH=690
PAGES={12:['ABC','DEF','GHI'],13:['JKL','MNO','PQR'],14:['STU','VWX','YZ1'],15:['234','567','890']}
def page(p):
    g=cv2.imread('/mnt/user-data/uploads/Image063017043007_%04d.jpg'%p,0)
    b=(g<170).astype(float)[200:1450,200:1150]
    ang=max((rot(b,a).sum(1).var(),a) for a in np.arange(-2,2.01,0.05))[1]
    inv=255-g
    o=cv2.morphologyEx(inv,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(6,6)))
    return rot(o.astype(float),ang),ang
def extract(pages=None):
    out={}
    for p,rows in (pages or PAGES).items():
        o,ang=page(p)
        bi=(o>128).astype(np.uint8); bi[:180]=0; bi[1460:]=0; bi[:,:190]=0
        n,lab,st,cen=cv2.connectedComponentsWithStats(bi)
        cs=[s[:4].tolist() for s in st[1:] if s[4]>2500 and s[3]>150]
        cs.sort(key=lambda s:s[1]+s[3]/2)
        bands=[]
        for s in cs:
            cy=s[1]+s[3]/2
            for bd in bands:
                if abs(bd['cy']-cy)<200: bd['c'].append(s); break
            else: bands.append({'cy':cy,'c':[s]})
        print(p,'angle',round(ang,2),[len(b['c']) for b in bands])
        for bd,letters in zip(bands,rows):
            g=[]
            for s in sorted(bd['c'],key=lambda s:s[0]):
                if g and s[0]<g[-1][0]+g[-1][2]-5:
                    a=g[-1]; x0=min(a[0],s[0]);y0=min(a[1],s[1]);x1=max(a[0]+a[2],s[0]+s[2]);y1=max(a[1]+a[3],s[1]+s[3]); g[-1]=[x0,y0,x1-x0,y1-y0]
                else: g.append(list(s))
            assert len(g)==len(letters),(p,letters,g)
            heights=[s[3] for s,l in zip(g,letters) if l not in 'Q']
            base=np.median([s[1]+s[3] for s,l in zip(g,letters) if l not in 'QCGOSJ0235689'] or [s[1]+s[3] for s in g])
            capp=np.median([s[3] for s,l in zip(g,letters) if l in 'BDEFHIKLMNPRTUVWXYZ14']) if any(l in 'BDEFHIKLMNPRTUVWXYZ14' for l in letters) else None
            for s,l in zip(g,letters): out[l]=dict(img=o,box=s,base=base,cap=capp)
    return out
def rasterize(d,cap_px):
    x,y,w,h=d['box']; m=10
    S=CAPH/cap_px
    c=d['img'][y-m:y+h+m,x-m:x+w+m].astype(np.float32)
    big=cv2.resize(c,None,fx=S,fy=S,interpolation=cv2.INTER_CUBIC)
    big=cv2.GaussianBlur(big,(0,0),0.6*S)
    b=(big>128).astype(np.uint8)
    n,lab,st,cen=cv2.connectedComponentsWithStats(b)
    lo,hi=m*S+3,(m+w)*S-3
    keep=[i for i in range(1,n) if st[i][4]>0.05*st[1:,4].max() and st[i][0]<hi and st[i][0]+st[i][2]>lo
          and not (st[i][0]+st[i][2]<lo+15*S or st[i][0]>hi-15*S)]
    b=np.isin(lab,keep).astype(np.uint8)
    # keep holes: fill by taking the main component only (holes remain as 0)
    f=cv2.GaussianBlur(b.astype(np.float32),(0,0),3); b=(f>0.5).astype(np.uint8)
    oy=(d['base']-(y-m))*S
    return b,10*S,oy
if __name__=='__main__':
    E=extract()
    caps=[d['cap'] for d in E.values() if d['cap']]
    CP=float(np.median(caps)); print('cap px',CP)
    A2={}
    for l,d in E.items():
        b,ox,oy=rasterize(d,CP); A2[l]=(b,ox,oy,'seriesD')
        ys,xs=np.nonzero(b); print(l,xs.max()-xs.min()+1,round(oy-ys.min()),round(oy-ys.max()))
    pickle.dump(A2,open('A2.pkl','wb'))
