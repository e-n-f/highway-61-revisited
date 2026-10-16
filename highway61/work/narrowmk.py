import numpy as np,pickle,cv2,sys,string
from family2 import glyph,M
from snaplines import refine,apply
from lighten import lighten
ML=pickle.load(open('master_lines.pkl','rb'))
def unit(v): return v/np.hypot(*v)
def constrain_narrow(l,F,pad,noy,f):
    oym=M[l][2]; lines=[]; pairs=[]
    cs,_=cv2.findContours(F.copy(),cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)
    P=np.vstack([c[:,0,:] for c in cs]).astype(float) if cs else np.zeros((0,2))
    H,W=F.shape
    def mp(p): return np.array([pad+p[0]*f, noy-(oym-p[1])])
    def ink(p):
        x,y=int(round(p[0])),int(round(p[1])); return 0<=x<W and 0<=y<H and F[y,x]>0
    for L in ML.get(l,[]):
        u=L['u']; ang=np.degrees(np.arctan2(abs(u[1]),abs(u[0])))
        axis=ang<3 or ang>87
        structural=(l in 'kvwxyz') and 35<=ang<=85
        p0=mp(L['c']+u*L['t0']); p1=mp(L['c']+u*L['t1']); c=(p0+p1)/2
        u2=unit(p1-p0) if (structural or axis) else u.copy()
        if u2@u<0: u2=-u2
        n2=np.array([-u2[1],u2[0]])
        if len(P)==0: continue
        d=(P-c)@n2; t=(P-c)@u2
        half=max(abs((p1-c)@u2),10)
        sel=(np.abs(d)<25)&(np.abs(t)<=half)
        if sel.sum()<5: continue
        off=np.median(d[sel]); c2=c+n2*off
        if structural or axis:
            t0,t1=-half,half
        else:
            # extent of the cut across the stroke: walk along the line just inside the ink
            side=L['side']; ts=np.arange(-300,301,2.0)
            inside=np.array([ink(c2+u2*tt+n2*side*6) for tt in ts])
            k0=np.argmin(np.abs(ts))
            if not inside[k0]:
                continue
            a=k0
            while a>0 and inside[a-1]: a-=1
            b=k0
            while b<len(ts)-1 and inside[b+1]: b+=1
            t0,t1=ts[a],ts[b]
        lines.append(dict(c=c2,u=u2,n=n2,t0=t0,t1=t1,side=L['side'],term=not (structural or axis))); pairs.append((L,lines[-1]))
    nb=apply(F,lines)
    se=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(7,7))
    nb=cv2.morphologyEx(cv2.morphologyEx(nb,cv2.MORPH_OPEN,se),cv2.MORPH_CLOSE,se)
    nb=apply(nb,lines)
    k=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(15,15))
    nb=(nb&cv2.dilate(F,k))|cv2.erode(F,k)      # never move an edge more than ~7 units
    for L in lines:
        if L['term']: nb=square_terminal(nb,L)
    nb=apply(nb,[L for L in lines if L['term']])
    sq=square_corners(nb,F,lines)
    k2=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(25,25))
    nb=sq&cv2.dilate(nb,k2)     # corner fixes may add at most ~12 units of ink
    nb=carve_notches(l,nb,noy,pairs)
    if l=='z': nb=rebuild_z(nb,lines,noy,pad,f,T_CUR)
    if l in 'tf': nb=rebuild_bar(nb,lines,noy)
    nb=special_fix(l,nb,lines,noy)
    nb=cv2.morphologyEx(nb,cv2.MORPH_OPEN,np.ones((5,5),np.uint8))   # drop hairline spikes
    return nb,len(lines)
T_CUR=0
SER={'B':(0.542,82),'C':(0.78,102),'D':(0.89,110),'F':(1.10,132)}
# ---------- z rebuilt from its straight edges ----------
def rebuild_z(nb,lines,noy,pad,f,T):
    """z as three simple shapes: two bars and a diagonal whose edges run
    into the bars' inner corners (top-right and bottom-left), as in the drawing."""
    b,ox,oym,_=M['z']; L=ML['z']
    H=[x for x in L if abs(x['u'][1])<0.05]; V=[x for x in L if abs(x['u'][0])<0.05]
    hs=sorted(set(round(oym-x['c'][1]) for x in H)); mids=[h for h in hs if 40<h<485]
    if len(mids)<2 or len(V)<2: return nb
    lbar_m,ubar_m=min(mids),max(mids)
    xlm=min(x['c'][0] for x in V); xrm=max(x['c'][0] for x in V)
    k=T/140.0
    ubar=525-(525-ubar_m)*k; lbar=lbar_m*k
    xl=pad+xlm*f; xr=pad+xrm*f
    R=lambda y: noy-y
    out=np.zeros_like(nb)
    def poly(pts): cv2.fillPoly(out,[np.round(np.array(pts)*16).astype(np.int32)],1,shift=4)
    poly([[xl,R(525)],[xr,R(525)],[xr,R(ubar)],[xl,R(ubar)]])
    poly([[xl,R(lbar)],[xr,R(lbar)],[xr,R(0)],[xl,R(0)]])
    dl=0.25*(525-ubar)                         # the diagonal starts a little into each bar
    yt,ybm=ubar+dl,lbar-dl
    dy=yt-ybm; h=T
    for _ in range(30):                      # horizontal thickness giving a perpendicular stroke of T
        run=(xr-h)-xl; th=np.arctan2(dy,run); h=T/np.sin(th)
    poly([[xl,R(ybm)],[xl+h,R(ybm)],[xr,R(yt)],[xr-h,R(yt)]])
    # cut the bars back flush with the diagonal where they would stick out past it
    H_,W_=out.shape; yy,xx=np.mgrid[0:H_,0:W_]; yu=noy-yy
    xr_edge=xl+h+(yu-ybm)*(xr-(xl+h))/(yt-ybm)          # diagonal's right edge
    xl_edge=xl+(yu-ybm)*((xr-h)-xl)/(yt-ybm)            # diagonal's left edge
    out[(yu>=ubar-1)&(yu<yt)&(xx>xr_edge)]=0
    out[(yu<=lbar+1)&(yu>ybm)&(xx<xl_edge)]=0
    return out

# ---------- t and f crossbars rebuilt square ----------
def rebuild_bar(nb,lines,noy):
    def yup(p): return noy-p[1]
    Hs=[L for L in lines if abs(L['u'][1])<0.05 and 300<yup(L['c'])<600]
    Vs=[L for L in lines if abs(L['u'][0])<0.05]
    if not Hs: return nb
    ys=sorted(yup(L['c']) for L in Hs)
    yb=ys[0]; yt=ys[-1] if ys[-1]-ys[0]>=40 else 525.0
    def yr(L):
        a=yup(L['c']+L['u']*L['t0']); b=yup(L['c']+L['u']*L['t1']); return min(a,b),max(a,b)
    R=lambda y:int(round(noy-y))
    # stem edges: scan the ink just below the bar (more reliable than whichever marks survived)
    row=nb[R(yb-80)]; xs=np.nonzero(row)[0]
    if len(xs)==0: return nb
    sl,sr=xs[0],xs[-1]
    # bar ends from the ink across the middle of the bar
    L_,R_=[],[]
    for y in np.arange(yb+15,yt-15,2):
        xs2=np.nonzero(nb[R(y)])[0]
        if len(xs2): L_.append(xs2[0]); R_.append(xs2[-1])
    if not L_: return nb
    xl=int(round(np.median(L_))); xr=int(round(np.median(R_)))
    out=nb.copy()
    # clear everything beside the stem in a band around the bar, then lay in an exact rectangle
    for c0,c1 in [(0,int(sl)),(int(sr)+1,out.shape[1])]:
        out[R(yt+30):R(yb-30)+1, c0:c1]=0
    out[R(yt):R(yb)+1, xl:xr+1]=1
    return out

# ---------- square outer stem corners (a foot, g top) and clean the k arm ----------
def special_fix(l,nb,lines,noy):
    H,W=nb.shape; R=lambda y:int(round(noy-y))
    Vs=[L for L in lines if abs(L['u'][0])<0.05]
    if l in 'ag' and Vs:
        xr=max(L['c'][0] for L in Vs); xi=int(round(xr))
        if l=='a':
            nb[R(150):R(-40), xi+1:]=0            # nothing to the right of the stem near the foot
            nb[R(0)+1:R(-40), xi-120:xi+1]=0       # nothing below the baseline under the stem
            nb[R(25):R(0)+1, xi-25:xi+1]=1         # square foot corner
        else:
            nb[R(560):R(380), xi+1:]=0
            nb[R(560):R(525), xi-120:xi+1]=0
            nb[R(525):R(500), xi-25:xi+1]=1
    return nb

# ---------- corner squaring ----------
def square_corners(nb,F,lines,r=30,ext=75):
    H,W=nb.shape; yy,xx=np.mgrid[0:H,0:W]; P=np.stack([xx,yy],-1).astype(float)
    out=nb.copy()
    for i in range(len(lines)):
        for j in range(i+1,len(lines)):
            A,B=lines[i],lines[j]
            if abs(A['u']@B['u'])>np.cos(np.radians(15)): continue
            Mx=np.array([A['u'],-B['u']]).T
            sa,sb=np.linalg.solve(Mx,B['c']-A['c'])
            def near_end(L,s): return (L['t0']-ext<=s<=L['t0']+ext) or (L['t1']-ext<=s<=L['t1']+ext)
            if not (near_end(A,sa) and near_end(B,sb)): continue
            X=A['c']+A['u']*sa
            x0,x1=int(max(0,X[0]-r)),int(min(W,X[0]+r+1)); y0,y1=int(max(0,X[1]-r)),int(min(H,X[1]+r+1))
            if x1<=x0 or y1<=y0: continue
            Q=P[y0:y1,x0:x1]; disk=np.hypot(Q[...,0]-X[0],Q[...,1]-X[1])<=r
            ia=((Q-A['c'])@A['n'])*A['side']>0; ib=((Q-B['c'])@B['n'])*B['side']>0
            small=np.hypot(Q[...,0]-X[0],Q[...,1]-X[1])<=12
            q1=small&ia&~ib; q2=small&~ia&ib
            if q1.sum()<5 or q2.sum()<5: continue
            f1=F[y0:y1,x0:x1][q1].mean(); f2=F[y0:y1,x0:x1][q2].mean()
            if f1>0.7 and f2>0.7: concave=True
            elif f1<0.3 and f2<0.3: concave=False
            else: continue
            ink=(~(~ia&~ib)) if concave else (ia&ib)
            sub=out[y0:y1,x0:x1]; sub[disk&ink]=1; sub[disk&~ink]=0
    return out

# ---------- notches and square terminal ends ----------
def master_notches(l):
    b,ox,oy,_=M[l]
    se=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(25,25))
    gap=(cv2.morphologyEx(b,cv2.MORPH_CLOSE,se)&(1-b)).astype(np.uint8)
    n,lab,st,cen=cv2.connectedComponentsWithStats(gap)
    out=[]
    for i in range(1,n):
        if not (30<st[i][4]<4000): continue
        yup=oy-cen[i][1]
        if not (yup<70 or 455<yup<600 or -260<yup<-150): continue
        ys,xs=np.nonzero(lab==i); out.append((xs,ys))
    return out
def carve_notches(l,nb,noy,pairs):
    if l not in 'abdghmnpqru': return nb
    b,ox,oy,_=M[l]
    for xs,ys in master_notches(l):
        cy=ys.mean(); cx=xs.mean(); best=None
        for Lm,Ln in pairs:
            u=Lm['u']
            if abs(u[1])<0.99: continue          # vertical edges only
            xl=Lm['c'][0]+u[0]*(cy-Lm['c'][1])/u[1]
            yr=sorted([Lm['c'][1]+u[1]*Lm['t0'],Lm['c'][1]+u[1]*Lm['t1']])
            dist=abs(cx-xl)
            if best is None or dist<best[0]: best=(dist,xl,Ln)
        if best is None or best[0]>60: continue
        _,xl,Ln=best
        xn=Ln['c'][0]
        X=np.round(xn+(xs-xl)).astype(int); Y=np.round(noy-(oy-ys)).astype(int)
        ok=(X>=0)&(X<nb.shape[1])&(Y>=0)&(Y<nb.shape[0])
        nb[Y[ok],X[ok]]=0
        # keep the notch open to the outside, as in the drawing
        yup=oy-cy
        for x in np.unique(X[ok]):
            col=Y[ok][X[ok]==x]
            if yup<70:   nb[col.max():int(round(noy+12)),x]=0          # down through the baseline
            elif yup>400: nb[int(round(noy-537)):col.min()+1,x]=0      # up through the x-height
    return nb
def square_terminal(nb,L):
    """extend the stroke sides straight to a terminal cut"""
    H,W=nb.shape
    def ink(p):
        x,y=int(round(p[0])),int(round(p[1])); return 0<=x<W and 0<=y<H and nb[y,x]>0
    c,u,n,side=L['c'],L['u'],L['n'],L['side']
    ts=np.arange(-300,301,1.0)
    def extent(s):
        ins=np.array([ink(c+u*t+n*side*s) for t in ts]); k=np.argmin(np.abs(ts-(L['t0']+L['t1'])/2))
        if not ins[k]: return None
        a=k
        while a>0 and ins[a-1]: a-=1
        bb=k
        while bb<len(ts)-1 and ins[bb+1]: bb+=1
        return ts[a],ts[bb]
    S=np.arange(18,40,3.0); E=[(s,extent(s)) for s in S]; E=[(s,e) for s,e in E if e]
    if len(E)<4: return nb
    ss=np.array([s for s,_ in E]); A=np.array([e[0] for _,e in E]); B=np.array([e[1] for _,e in E])
    ka=np.polyfit(ss,A,1); kb=np.polyfit(ss,B,1)
    poly=[c+u*np.polyval(ka,0)+n*side*0, c+u*np.polyval(kb,0), c+u*np.polyval(kb,22)+n*side*22, c+u*np.polyval(ka,22)+n*side*22]
    pts=np.round(np.array(poly)*16).astype(np.int32)
    cv2.fillPoly(nb,[pts],1,shift=4)
    return nb

if __name__=='__main__':
    s=sys.argv[1]; f,T=SER[s]; T_CUR=T; out={}; info=[]
    for l in string.ascii_lowercase:
        ff=0.66 if (s=='B' and l=='x') else f
        F,pad,noy=glyph(l,ff,T)
        nb,n=constrain_narrow(l,F,pad,noy,ff); out[l]=(nb,pad,noy,'mk-'+s); info.append('%s%d'%(l,n))
    out.update(lighten({'w':out['w']},4))
    if s=='B': out.update(lighten({'x':out['x']},4))
    d=pickle.load(open('lcE_%s.pkl'%s,'rb')); d.update(out); pickle.dump(d,open('lcE_%s.pkl'%s,'wb'))
    print(s,' '.join(info))

