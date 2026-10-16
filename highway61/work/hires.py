import cv2,numpy as np,pickle,sys
from src import rot
LAY={'0016':['ABC','DEF','GHI'],'0017':['JKL','MNO','PQR'],'0018':['STU','VWX','YZ1'],'0019':['234','567','890'],
     '0031':['abc','def','ghi'],'0032':['jkl','mno','pqr'],'0033':['stu','vwx','yz']}
def page(p):
    g=cv2.imread('hr/%s.png'%p,0).astype(np.float32)
    sm=cv2.resize(g,(1275,int(g.shape[0]*1275/g.shape[1])),interpolation=cv2.INTER_AREA)
    b=(sm<170).astype(float)[200:1450,200:1150]
    ang=max((rot(b,a).sum(1).var(),a) for a in np.arange(-2,2.01,0.05))[1]
    k=g.shape[1]/1275
    o=rot(255-g,ang)
    H,W=o.shape
    dark=(o>100).astype(np.float32)
    rowf=dark[:,int(W*.1):int(W*.95)].mean(1); colf=dark[int(H*.1):int(H*.9),:].mean(0)
    D=6
    def lines(prof):
        L=[]
        for i in range(D+2,len(prof)-D-2):
            if prof[i]>=prof[i-1] and prof[i]>=prof[i+1] and prof[i]-max(prof[i-D-2],prof[i+D+2])>0.15: L.append(i)
        return L
    o2=o.copy()
    for i in lines(rowf):
        a,b=i-3,i+4
        top=o[a-4,:]; bot=o[b+3,:]; n=b-a
        for j in range(n): o2[a+j,:]=top+(bot-top)*(j+4)/(n+7)
    o3=o2.copy()
    for i in lines(colf):
        a,b=i-3,i+4
        lf=o2[:,a-4]; rt=o2[:,b+3]; n=b-a
        for j in range(n): o3[:,a+j]=lf+(rt-lf)*(j+4)/(n+7)
    O=cv2.morphologyEx(o,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(11,11)))
    Id=cv2.dilate(o3,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(9,9)))
    return np.minimum(O,Id),ang,k
def boxes(o,rows):
    H,W=o.shape
    bi=(o>128).astype(np.uint8)
    bi[:int(H*0.085)]=0; bi[int(H*0.90):]=0; bi[:,:int(W*0.075)]=0; bi[:,int(W*0.985):]=0
    n,lab,st,_=cv2.connectedComponentsWithStats(bi)
    cs=[list(s[:4])+[s[4]] for s in st[1:] if s[4]>4000]
    big=[c for c in cs if c[4]>40000]; small=[c for c in cs if c[4]<=40000]
    # merge dots into the letter below them
    for d in small:
        cx=d[0]+d[2]/2
        best=min(big,key=lambda c:abs(c[0]+c[2]/2-cx)+abs(c[1]-(d[1]+d[3])))
        x0=min(best[0],d[0]);y0=min(best[1],d[1]);x1=max(best[0]+best[2],d[0]+d[2]);y1=max(best[1]+best[3],d[1]+d[3])
        best[:4]=[x0,y0,x1-x0,y1-y0]
    big.sort(key=lambda s:s[1]+s[3]/2)
    bands=[]
    for s in big:
        cy=s[1]+s[3]-150  # use bottom-ish
        for bd in bands:
            if abs(bd['cy']-cy)<500: bd['c'].append(s); break
        else: bands.append({'cy':cy,'c':[s]})
    bands.sort(key=lambda b:b['cy'])
    out={}
    assert len(bands)==len(rows),[len(b['c']) for b in bands]
    for bd,letters in zip(bands,rows):
        g=sorted(bd['c'],key=lambda s:s[0]); assert len(g)==len(letters),(letters,g)
        for s,l in zip(g,letters): out[l]=s[:4]
    return out
def run(pages,DESC,ROUND,heightset,target):
    E={}
    for p in pages:
        o,ang,k=page(p); bx=boxes(o,LAY[p]); print(p,'angle',ang,len(bx))
        rowsof={}
        for r in LAY[p]:
            base=np.median([bx[l][1]+bx[l][3] for l in r if l not in DESC and l not in ROUND] or [bx[l][1]+bx[l][3] for l in r if l not in DESC])
            for l in r: E[l]=dict(img=o,box=bx[l],base=base)
    hs=[E[l]['base']-E[l]['box'][1] for l in heightset if l in E]
    hp=float(np.median(hs)); print('ref height px',hp)
    S=target/hp
    out={}
    for l,d in E.items():
        x,y,w,h=d['box']; m=30
        c=d['img'][y-m:y+h+m, x-m:x+w+m]
        big=cv2.resize(c,None,fx=S,fy=S,interpolation=cv2.INTER_CUBIC)
        big=cv2.GaussianBlur(big,(0,0),1.2)
        b=(big>128).astype(np.uint8)
        n,lab,st,_=cv2.connectedComponentsWithStats(b)
        lo,hi=m*S+4,(m+w)*S-4
        keep=[i for i in range(1,n) if st[i][4]>200 and st[i][0]<hi and st[i][0]+st[i][2]>lo and st[i][1]<(m+h)*S and st[i][1]+st[i][3]>m*S]
        b=np.isin(lab,keep).astype(np.uint8)
        out[l]=(b,m*S,(d['base']-(y-m))*S,'hires')
    return out,S
if __name__=='__main__':
    which=sys.argv[1]
    if which=='caps':
        out,S=run(['0016','0017','0018','0019'],'Q','CGOSJ0235689','BDEFHIKLMNPRTUVWXYZ14',690)
        pickle.dump(out,open('hr_caps_E.pkl','wb'))
    else:
        out,S=run(['0031','0032','0033'],'gjpqy','aceos','mnruvwxz',525)
        pickle.dump(out,open('hr_lc.pkl','wb'))
    for l in sorted(out):
        b,ox,oy,_=out[l]; ys,xs=np.nonzero(b); print(l,xs.max()-xs.min()+1,round(oy-ys.min()),round(oy-ys.max()),end=' | ')
