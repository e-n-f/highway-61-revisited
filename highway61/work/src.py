import cv2,numpy as np
U='/mnt/user-data/uploads/'
def sign(name):
    if name=='bf':
        im=cv2.imread(U+'Screenshot_2026-10-08_at_2_12_28_PM.png'); box=(150,400,1540,660)
    else:
        im=cv2.imread(U+'Screenshot_2026-09-30_at_10_49_11_PM_copy.png'); box=(40,70,600,320)
    x0,y0,x1,y1=box
    c=im[y0:y1,x0:x1].astype(float)
    r=c[...,2]-0.5*c[...,1]
    r=(r-r.min())/(r.max()-r.min())*255
    return r  # bright = ink
def spec(page):
    g=cv2.imread(U+'Image063017043007_00%d.jpg'%page,0).astype(np.uint8)
    inv=255-g
    k=cv2.getStructuringElement(cv2.MORPH_RECT,(4,4))
    o=cv2.morphologyEx(inv,cv2.MORPH_OPEN,k)
    return o.astype(float)

def rot(img,deg):
    h,w=img.shape
    M=cv2.getRotationMatrix2D((w/2,h/2),deg,1)
    return cv2.warpAffine(img,M,(w,h),flags=cv2.INTER_CUBIC,borderValue=0)
def grid_angle(page):
    g=cv2.imread(U+'Image063017043007_00%d.jpg'%page,0)
    b=(g<170).astype(float)[60:780,80:620]
    best=None
    for a in np.arange(-2,2.01,0.05):
        r=rot(b,a); s=r.sum(1); v=s.var()
        if best is None or v>best[0]: best=(v,a)
    return best[1]

LAYOUT={31:['abc','def','ghi'],32:['jkl','mno','pqr'],33:['stu','vwx','yz']}
DESC=set('gjpqy'); XH=set('acemnorsuvwxz')
def spec_letters():
    out={}
    for p,rows in LAYOUT.items():
        a=grid_angle(p)
        img=rot(spec(p),a)
        b=(img>128).astype(np.uint8)
        b[:70,:]=0; b[770:,:]=0; b[:,:60]=0  # kill header/footer/binding
        n,lab,st,cen=cv2.connectedComponentsWithStats(b)
        comps=[s for s in st[1:].tolist() if s[4]>300]
        # group by row using y-center
        comps.sort(key=lambda s:s[1]+s[3]/2)
        bands=[]
        for s in comps:
            cy=s[1]+s[3]/2
            for bd in bands:
                if abs(bd['cy']-cy)<110: bd['c'].append(s); break
            else: bands.append({'cy':cy,'c':[s]})
        assert len(bands)==len(rows),(p,len(bands))
        for bd,letters in zip(bands,rows):
            cs=sorted(bd['c'],key=lambda s:s[0])
            # merge dots: components overlapping in x
            groups=[]
            for s in cs:
                if groups and s[0]<groups[-1][0]+groups[-1][2]-5 or (groups and s[0]+s[2]<=groups[-1][0]+groups[-1][2]):
                    g=groups[-1]; x0=min(g[0],s[0]);y0=min(g[1],s[1]);x1=max(g[0]+g[2],s[0]+s[2]);y1=max(g[1]+g[3],s[1]+s[3])
                    groups[-1]=[x0,y0,x1-x0,y1-y0]
                else: groups.append(s[:4])
            assert len(groups)==len(letters),(p,letters,groups)
            base=np.median([g[1]+g[3] for g,l in zip(groups,letters) if l not in DESC])
            for g,l in zip(groups,letters): out[l]=dict(img=img,box=g,base=base,page=p)
        # x-height per page
    return out

def sign_letters(name):
    r=sign(name)
    t,_=cv2.threshold(r.astype(np.uint8),0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    b=(r>t).astype(np.uint8)
    n,lab,st,cen=cv2.connectedComponentsWithStats(b)
    comps=[s[:4] for s in st[1:].tolist() if s[4]>150 and s[3]<200 and s[2]<150]
    if name=='bf':
        comps=[c for c in comps if c[3]>60]
        comps.sort(key=lambda s:s[0])
        names='Bayshor eFwy'.replace(' ','')
        # names: B a y s h o r e F w y
        names=['B','a','y','s','h','o','r','e','F','w','y2']
        lines=[comps]
    else:
        comps=[c for c in comps if c[0]>20]
        l1=sorted([c for c in comps if c[1]<130],key=lambda s:s[0])
        l2=sorted([c for c in comps if c[1]>=130],key=lambda s:s[0])
        lines=[l1,l2]
        names=None
    # merge dots
    res=[]
    for line in lines:
        groups=[]
        for s in line:
            if groups and s[0]<groups[-1][0]+groups[-1][2]-3:
                g=groups[-1]; x0=min(g[0],s[0]);y0=min(g[1],s[1]);x1=max(g[0]+g[2],s[0]+s[2]);y1=max(g[1]+g[3],s[1]+s[3])
                groups[-1]=[x0,y0,x1-x0,y1-y0]
            else: groups.append(list(s))
        res.append(groups)
    return r,res,t

XH_UNITS=525
def raster(gray, box, base, xh, thr=None, margin=6, sigma=0.5):
    x,y,w,h=box
    x0=max(0,x-margin); y0=max(0,y-margin); x1=x+w+margin; y1=y+h+margin
    c=gray[y0:y1,x0:x1].astype(np.float32)
    S=XH_UNITS/xh
    big=cv2.resize(c,None,fx=S,fy=S,interpolation=cv2.INTER_CUBIC)
    if sigma>0:
        big=cv2.GaussianBlur(big,(0,0),sigma*S)
    if thr is None:
        ink=np.percentile(c,97); bg=np.percentile(c,25)
        thr=(ink+bg)/2
    bi=(big>thr).astype(np.uint8)
    # keep components overlapping the original box region (drop neighbors' fragments)
    n,lab,st,cen=cv2.connectedComponentsWithStats(bi)
    keep=np.zeros_like(bi)
    bx0=(x-x0)*S; bx1=(x+w-x0)*S
    for i in range(1,n):
        cx=st[i][0]+st[i][2]/2
        if st[i][4]>200 and bx0-5<cx<bx1+5: keep[lab==i]=1
    # origin: left at box left, baseline
    ox=(x-x0)*S; oy=(base-y0)*S
    return keep, ox, oy, S
def stem_width(b, oy, frac=0.5):
    row=b[int(oy-XH_UNITS*frac)]
    runs=[];cur=0
    for v in row:
        if v: cur+=1
        elif cur: runs.append(cur); cur=0
    if cur: runs.append(cur)
    return runs

def vruns(b, xfrac=0.5):
    ys,xs=np.nonzero(b)
    x=int(xs.min()+(xs.max()-xs.min())*xfrac)
    col=b[:,x]; runs=[];cur=0;start=0
    for i,v in enumerate(col):
        if v and not cur: start=i
        if v: cur+=1
        elif cur: runs.append((start,cur)); cur=0
    if cur: runs.append((start,cur))
    return [r[1] for r in runs]

STEM_T=115; STEM_S=137; VTHIN=20
def synth(d, fx, xh=115.5, hstem_target=STEM_T, verode=VTHIN):
    # high-res smooth grayscale, then transform
    x,y,w,h=d['box']; m=8
    gray=d['img']; base=d['base']
    c=gray[y-m:y+h+m, x-m:x+w+m].astype(np.float32)
    S=XH_UNITS/xh
    big=cv2.resize(c,None,fx=S,fy=S,interpolation=cv2.INTER_CUBIC)
    big=cv2.GaussianBlur(big,(0,0),0.5*S)
    oy=(base-(y-m))*S
    # vertical pre-scale: y_up' = a*y_up - verode/2 ; with y_up = oy - row
    half=verode/2
    a=(XH_UNITS+verode)/XH_UNITS
    H,W=big.shape
    newoy=oy*a+half+20
    newH=int(H*a+half*2+40)
    # map: row' -> y_up' = newoy-row'; y_up=(y_up'+half)/a ; row = oy - y_up
    rows=np.arange(newH,dtype=np.float32)
    src_r=oy-((newoy-rows)+half)/a
    newW=int(W*fx)+60
    cols=np.arange(newW,dtype=np.float32)
    src_c=(cols-30)/fx
    mx,my=np.meshgrid(src_c,src_r)
    out=cv2.remap(big,mx.astype(np.float32),my.astype(np.float32),cv2.INTER_LINEAR,borderValue=0)
    b=(out>128).astype(np.uint8)
    # vertical erosion
    if verode>0:
        b=cv2.erode(b,np.ones((int(verode)+1,1),np.uint8))
    dl=hstem_target-STEM_S*fx
    k=int(round(abs(dl)))
    if k>=1:
        se=np.ones((1,k+1),np.uint8)
        b=cv2.dilate(b,se) if dl>0 else cv2.erode(b,se)
    return b, 30, newoy
def width(b):
    xs=np.nonzero(b.any(0))[0]; return xs[-1]-xs[0]+1
