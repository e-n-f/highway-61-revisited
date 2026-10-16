"""Brush-geometry construction of the remaining lowercase from 1966 measurements."""
import numpy as np,cv2,pickle
from scipy.optimize import least_squares
from geomfit import M,ML
TM=140.0
def unit(v): return v/np.hypot(*v)
# ---------- measurement helpers (master frame: x, yup) ----------
def master(l):
    b,ox,oy,_=M[l]; return b,oy
def edge_pts(l,xr,yr,kind):
    b,oy=master(l); H,W=b.shape
    up=np.zeros_like(b); up[1:]=b[:-1]; dn=np.zeros_like(b); dn[:-1]=b[1:]
    lf=np.zeros_like(b); lf[:,1:]=b[:,:-1]; rt=np.zeros_like(b); rt[:,:-1]=b[:,1:]
    m={'top':(b>0)&(up==0),'bottom':(b>0)&(dn==0),'left':(b>0)&(lf==0),'right':(b>0)&(rt==0)}[kind]
    ys,xs=np.nonzero(m); yu=oy-ys
    k=(xs>=xr[0])&(xs<=xr[1])&(yu>=yr[0])&(yu<=yr[1])
    return np.column_stack([xs[k],yu[k]]).astype(float)
def fit_ellipse(P,guess=None):
    x,y=P[:,0],P[:,1]
    p0=guess or [x.mean(),y.min() if True else y.mean(),(x.max()-x.min())/2+1,(y.max()-y.min())+1]
    def res(p):
        cx,cy,a,b=p; return (np.sqrt(((x-cx)/a)**2+((y-cy)/b)**2)-1)*min(a,b)
    s=least_squares(res,p0,bounds=([-1e4,-1e4,5,5],[1e4,1e4,1e4,1e4]))
    return s.x,np.abs(s.fun).max()
def ring_fit(l,xr,yr,outer_kind,inner_kind,guess=None):
    po,eo=fit_ellipse(edge_pts(l,xr,yr,outer_kind),guess)
    pi,ei=fit_ellipse(edge_pts(l,xr,yr,inner_kind),[po[0],po[1],po[2]*0.5,po[3]*0.5])
    return dict(cx=(po[0]+pi[0])/2,cy=(po[1]+pi[1])/2,a=(po[2]+pi[2])/2,b=(po[3]+pi[3])/2,
                otop=po[1]+po[3],obot=po[1]-po[3],err=(eo,ei))
def stems(l):
    """vertical marked edges -> list of (left,right) edge pairs (master x)"""
    xs=sorted(L['c'][0] for L in ML[l] if abs(L['u'][1])>0.99)
    xs=[x for i,x in enumerate(xs) if i==0 or x-xs[i-1]>8]
    pairs=[]
    i=0
    while i+1<len(xs):
        if 60<xs[i+1]-xs[i]<200: pairs.append((xs[i],xs[i+1])); i+=2
        else: i+=1
    return pairs
def term_line(l,which):
    """marked slanted terminal: which='top'|'bottom'|'left'|'right' (most extreme)"""
    b,oy=master(l)
    TL=[L for L in ML[l] if 3<np.degrees(np.arctan2(abs(L['u'][1]),abs(L['u'][0])))<87 and not
        (l in 'kvwxyz' and np.degrees(np.arctan2(abs(L['u'][1]),abs(L['u'][0])))>=35)]
    if not TL: return None
    key={'top':lambda L:-(oy-L['c'][1]),'bottom':lambda L:oy-L['c'][1],'left':lambda L:L['c'][0],'right':lambda L:-L['c'][0]}[which]
    L=min(TL,key=key); return dict(p=np.array([L['c'][0],oy-L['c'][1]]),u=np.array([L['u'][0],-L['u'][1]]))
# ---------- drawing in the output frame ----------
class Canvas:
    def __init__(s,f,T,pad=40):
        s.f,s.T,s.pad=f,T,pad; s.H=1300; s.noy=900.0; s.W=int(650*max(f,1)+2*pad+200)
        s.img=np.zeros((s.H,s.W),np.uint8)
        yy,xx=np.mgrid[0:s.H,0:s.W]; s.xx=xx; s.yu=s.noy-yy
    def X(s,x): return s.pad+x*s.f
    def R(s,y): return s.noy-y
    def mask(s): return np.zeros((s.H,s.W),np.uint8)
    def stem(s,xc,y0,y1,top=None,bot=None):
        """vertical stroke of width T centred on master x xc, from y0 to y1 (y-up); optional slanted cuts"""
        x=s.X(xc); m=(s.xx>=x-s.T/2)&(s.xx<x+s.T/2)&(s.yu>=y0)&(s.yu<=y1)
        if top is not None: m&=s.below(top)
        if bot is not None: m&=~s.below(bot)
        return m
    def below(s,L):
        """points below a (terminal) line given in master coords; the angle is kept, the point is mapped"""
        P=np.array([s.X(L['p'][0]),L['p'][1]]); u=L['u']
        if abs(u[0])<1e-6: return s.xx<P[0]
        mm=u[1]/u[0]; return s.yu<=P[1]+mm*(s.xx-P[0])
    def right_of(s,L):
        P=np.array([s.X(L['p'][0]),L['p'][1]]); u=L['u']
        if abs(u[1])<1e-6: return s.yu<P[1]
        return s.xx>=P[0]+(s.yu-P[1])*u[0]/u[1]
    def ring(s,rf,k=1.0):
        """brush ellipse: centre-line from the fit, x mapped, outer top/bottom kept; k scales the pen"""
        L=s.X(rf['cx']-rf['a']); R=s.X(rf['cx']+rf['a']); cx=(L+R)/2; ac=(R-L)/2
        T=s.T*k
        bc=rf['b']+( (TM*k-T)/2 )
        img=s.mask()
        cv2.ellipse(img,(int(round(cx*16)),int(round(s.R(rf['cy'])*16))),(int(round(ac*16)),int(round(bc*16))),0,0,360,1,
                    thickness=max(1,int(round(T))),shift=4)
        return img>0, dict(cx=cx,cy=rf['cy'],ac=ac,bc=bc)
    def diag(s,p0,p1,y0,y1):
        """straight stroke of perpendicular width T along the master centre-line p0->p1 (x mapped), cut flat at y0,y1"""
        a=np.array([s.X(p0[0]),p0[1]]); b=np.array([s.X(p1[0]),p1[1]]); d=unit(b-a); n=np.array([-d[1],d[0]])
        a2=a-d*1500; b2=b+d*1500; h=s.T/2
        poly=np.array([a2+n*h,b2+n*h,b2-n*h,a2-n*h]); poly[:,1]=s.noy-poly[:,1]
        img=s.mask(); cv2.fillPoly(img,[np.round(poly*16).astype(np.int32)],1,shift=4)
        return (img>0)&(s.yu>=y0)&(s.yu<=y1)
    def diag_taper(s,p0,p1,y0,y1,wb,wt):
        """straight stroke whose perpendicular width runs from wb*T at y0 to wt*T at y1"""
        a=np.array([s.X(p0[0]),p0[1]]); b=np.array([s.X(p1[0]),p1[1]]); d=unit(b-a); n=np.array([-d[1],d[0]])
        at=lambda y: a+d*((y-a[1])/d[1])
        c0=at(y0-2); c1=at(y1+2)
        poly=np.array([c0+n*wb*s.T/2,c1+n*wt*s.T/2,c1-n*wt*s.T/2,c0-n*wb*s.T/2]); poly[:,1]=s.noy-poly[:,1]
        img=s.mask(); cv2.fillPoly(img,[np.round(poly*16).astype(np.int32)],1,shift=4)
        return (img>0)&(s.yu>=y0)&(s.yu<=y1)
    def diag_sides(s,p0,p1,y0,y1,bl,br,tl,tr):
        """stroke with separate left/right half-widths (fractions of T/2) at the bottom (y0) and top (y1)"""
        a=np.array([s.X(p0[0]),p0[1]]); b=np.array([s.X(p1[0]),p1[1]]); d=unit(b-a); n=np.array([-d[1],d[0]])
        if n[0]>0: n=-n                                   # n points left
        at=lambda y: a+d*((y-a[1])/d[1])
        c0=at(y0-2); c1=at(y1+2); h=s.T/2
        poly=np.array([c0+n*bl*h,c1+n*tl*h,c1-n*tr*h,c0-n*br*h]); poly[:,1]=s.noy-poly[:,1]
        img=s.mask(); cv2.fillPoly(img,[np.round(poly*16).astype(np.int32)],1,shift=4)
        return (img>0)&(s.yu>=y0)&(s.yu<=y1)
    def dot(s,l):
        b,oy=master(l); n,lab,st,cen=cv2.connectedComponentsWithStats(b)
        main=np.argmax(st[1:,4])+1
        for i in range(1,n):
            if i==main or st[i][4]<500: continue
            r=(st[i][2]+st[i][3])/4*(s.T/TM); cx=s.X(cen[i][0]); cy=s.R(oy-cen[i][1])
            img=s.mask(); cv2.circle(img,(int(round(cx*16)),int(round(cy*16))),int(round(r*16)),1,-1,shift=4); return img>0
        return s.mask()>0
def runs_at(l,y):
    b,oy=master(l); row=b[int(round(oy-y))]; xs=np.nonzero(row)[0]
    if len(xs)==0: return []
    sp=np.where(np.diff(xs)!=1)[0]
    return [(a,z) for a,z in zip(np.r_[xs[0],xs[sp+1]],np.r_[xs[sp],xs[-1]])]
def sc(r): return (r[0]+r[1])/2
def ink_x(l):
    b,oy=master(l); xs=np.nonzero(b.any(0))[0]; return xs[0],xs[-1]
def build(l,f,T):
    C=Canvas(f,T); m=C.mask()>0
    ASC=700; DESC=-300
    if l in 'nhm':
        st=runs_at(l,150)
        m|=C.stem(sc(st[0]),0,525 if l!='h' else ASC,top=term_line(l,'top') if l=='h' else None)
        arches=[(st[0],st[1])] if l!='m' else [(st[0],st[1]),(st[1],st[2])]
        for a,bb in arches:
            c0,c1=sc(a),sc(bb)
            P=edge_pts(l,(a[1]+6,bb[1]),(300,600),'top')
            # centre-line runs from stem centre to stem centre; fit only the height
            x,y=P[:,0],P[:,1]; ac=(c1-c0)/2; cx=(c0+c1)/2
            def res(p):
                cy,b=p; return (np.sqrt(((x-cx)/(ac+70))**2+((y-cy)/(b+65))**2)-1)*100
            cy,bb_=least_squares(res,[300,170]).x
            otop=cy+bb_+65
            bc=max(bb_+(TM-T)/2, ac*f)           # keep the counter round, not flat
            rf=dict(cx=cx,cy=otop-T/2-bc,a=ac,b=bc-(TM-T)/2)
            r,g=C.ring(rf); m|=r&(C.yu>=g['cy'])&(C.xx>=C.X(c0))
            xr=g['cx']+g['ac']
            m|=(C.xx>=xr-T/2)&(C.xx<xr+T/2)&(C.yu>=0)&(C.yu<=g['cy'])
    elif l=='u':
        st=runs_at(l,450); c0,c1=sc(st[0]),sc(st[1])
        P=edge_pts(l,(st[0][0],st[1][0]-6),(-30,250),'bottom'); x,y=P[:,0],P[:,1]; ac=(c1-c0)/2; cx=(c0+c1)/2
        def res(p):
            cy,b=p; return (np.sqrt(((x-cx)/(ac+70))**2+((y-cy)/(b+65))**2)-1)*100
        cy,bb_=least_squares(res,[200,170]).x
        obot=cy-bb_-65
        bc=max(bb_+(TM-T)/2, ac*f)
        rf=dict(cx=cx,cy=obot+T/2+bc,a=ac,b=bc-(TM-T)/2)
        r,g=C.ring(rf); m|=r&(C.yu<=g['cy'])&(C.xx<=C.X(c1))
        m|=C.stem(c0,g['cy'],525)
        m|=C.stem(c1,0,525)
    elif l=='r':
        st=runs_at(l,150); x0,x1=ink_x(l); c0=sc(st[0])
        rf=fit_anchored(edge_pts(l,(st[0][1]+6,x1),(380,600),'top'),c0,+1,[150,330,150])
        r,g=C.ring(rf); tl=term_line(l,'right')
        m|=r&(C.yu>=g['cy'])&(C.xx>=C.X(c0))&~C.right_of(tl)
        m|=C.stem(c0,0,525)
    elif l=='a':
        st=runs_at(l,150); x0,x1=ink_x(l); s=st[-1]; c1=sc(s)
        m|=C.stem(c1,0,525)
        b_,oy_=master(l)
        cs,hh=cv2.findContours(b_.copy(),cv2.RETR_CCOMP,cv2.CHAIN_APPROX_NONE)
        hole=max([c[:,0,:] for i,c in enumerate(cs) if hh[0][i][3]!=-1],key=len)
        Hp=np.column_stack([hole[:,0],oy_-hole[:,1]]).astype(float); Hp=Hp[Hp[:,0]<s[0]-6]
        Po=np.vstack([edge_pts(l,(x0,s[0]-6),(-30,200),'bottom'),edge_pts(l,(x0,s[0]-6),(-30,300),'left'),edge_pts(l,(x0+30,s[0]-6),(220,380),'top')])
        po,_=fit_ellipse(Po,[x0+150,130,150,140]); pi,_=fit_ellipse(Hp,[po[0],po[1],po[2]*0.5,po[3]*0.5])
        rb=dict(cx=(po[0]+pi[0])/2,cy=(po[1]+pi[1])/2,a=(po[2]+pi[2])/2,b=(po[3]+pi[3])/2)
        r,g=C.ring(rb,0.88); m|=r&(C.xx<=C.X(c1))
        Pa=edge_pts(l,(x0,s[1]),(400,600),'top')
        ra=fit_anchored(Pa,c1,-1,[150,330,150])
        r2,g2=C.ring(ra); tl=term_line(l,'top')
        m|=r2&(C.yu>=g2['cy'])&(C.xx<=C.X(c1))&C.right_of(tl)
    elif l in 'gj':
        st=runs_at(l,250) if l=='g' else runs_at(l,300); s=st[-1]; x0,x1=ink_x(l)
        if l=='g':
            b_,oy_=master(l)
            cs,hh=cv2.findContours(b_.copy(),cv2.RETR_CCOMP,cv2.CHAIN_APPROX_NONE)
            hole=max([c[:,0,:] for i,c in enumerate(cs) if hh[0][i][3]!=-1],key=len)
            Hp=np.column_stack([hole[:,0],oy_-hole[:,1]]).astype(float); Hp=Hp[Hp[:,0]<s[0]-12]
            Po=np.vstack([edge_pts(l,(x0,s[0]-12),(0,560),k) for k in ('top','left','bottom')])
            Po=Po[(Po[:,1]>20)|(Po[:,0]<x0+200)]
            po,_=fit_ellipse(Po,[x0+200,262,200,262]); pi,_=fit_ellipse(Hp,[po[0],po[1],po[2]*0.45,po[3]*0.5])
            rb=dict(cx=(po[0]+pi[0])/2,cy=(po[1]+pi[1])/2,a=(po[2]+pi[2])/2,b=(po[3]+pi[3])/2)
            rr,gg=C.ring(rb); m|=rr&(C.xx<=C.X(sc(s)))
        b_,oy_=master(l); ys_,xs_=np.nonzero(b_); yu_=oy_-ys_
        tail=yu_<-40; xl_=xs_[tail].min(); ye_=yu_.min()
        a_=sc(s)-(xl_+40); cxh=xl_+40; cyh=ye_+65+0.9*a_
        rt=dict(cx=cxh,cy=cyh,a=a_,b=cyh-(ye_+65))
        r,g=C.ring(rt); tl=term_line(l,'bottom') if l=='g' else None
        keep=r&(C.yu<=g['cy'])&(C.xx<=C.X(sc(s)))&(C.xx>=C.X(xl_))
        if l=='g' and tl is not None: keep&=C.right_of(tl)
        if l=='j':   # j's hook ends at a vertical cut
            keep&=C.xx>=C.X(x0)
        m|=keep
        m|=C.stem(sc(s),g['cy']-1,525)
        if l=='j': m|=C.dot('j')
    elif l in 'tf':
        st=runs_at(l,250); s=max(st,key=lambda r:r[1]-r[0]); x0,x1=ink_x(l)
        if l=='t':
            rt=hook(l,s,'bottom-right')
            r,g=C.ring(rt); tl=term_line(l,'bottom')
            m|=r&(C.yu<=g['cy'])&(C.xx>=C.X(sc(s)))&~C.right_of(tl) if tl is not None else r&(C.yu<=g['cy'])&(C.xx>=C.X(sc(s)))
            m|=C.stem(sc(s),g['cy']-1,ASC,top=term_line(l,'top'))
            br=max(runs_at(l,460),key=lambda r:r[1]-r[0]); th=131*T/TM
            m|=(C.yu>=527-th)&(C.yu<=527)&(C.xx>=C.X(br[0]))&(C.xx<=C.X(br[1]))
        else:
            rt=hook(l,s,'top-right')
            r,g=C.ring(rt); tl=term_line(l,'right')
            keep=r&(C.yu>=g['cy'])&(C.xx>=C.X(sc(s)))
            if tl is not None: keep&=~C.right_of(tl)
            else: keep&=C.xx<=C.X(x1)
            m|=keep
            m|=C.stem(sc(s),0,g['cy']+1)
            br=max(runs_at(l,455),key=lambda r:r[1]-r[0]); th=141*T/TM
            m|=(C.yu>=525-th)&(C.yu<=525)&(C.xx>=C.X(br[0]))&(C.xx<=C.X(br[1]))
    elif l in 'il':
        st=runs_at(l,150)[0]
        if l=='i': m|=C.stem(sc(st),0,525); m|=C.dot('i')
        else: m|=C.stem(sc(st),0,ASC,top=term_line(l,'top'))
    elif l in 'kvwxy':
        D=[L for L in ML[l] if 35<=np.degrees(np.arctan2(abs(L['u'][1]),abs(L['u'][0])))<87]
        b,oy=master(l)
        # pair parallel edges into strokes (the band between them must be ink)
        H_,W_=b.shape
        def inkfrac(A,B):
            c=(A['c']+B['c'])/2; u=A['u']; ts=np.linspace(-120,120,25); ok=0; n=0
            for t in ts:
                p=c+u*t; x,y=int(round(p[0])),int(round(p[1]))
                if 0<=x<W_ and 0<=y<H_: n+=1; ok+=b[y,x]>0
            return ok/max(1,n)
        cand=[]
        for i,A in enumerate(D):
            for j,B in enumerate(D):
                if j<=i or abs(abs(A['u']@B['u'])-1)>0.003: continue
                dist=abs((B['c']-A['c'])@A['n'])
                if 60<dist<200 and inkfrac(A,B)>0.85: cand.append((abs(dist-140),i,j))
        cand.sort(); used=set(); strokes=[]
        for dist,i,j in cand:
            if i in used or j in used: continue
            used|={i,j}; A,B=D[i],D[j]; c=(A['c']+B['c'])/2; u=A['u']
            p0=np.array([c[0],oy-c[1]]); d=np.array([u[0],-u[1]])
            low=min(oy-(A['c']+A['u']*t)[1] for A_ in (A,B) for t in (A_['t0'],A_['t1']))
            strokes.append((p0,p0+d*100,low))
        if l=='k':
            st=runs_at(l,150)[0]; m|=C.stem(sc(st),0,ASC,top=term_line(l,'top'))
        if l=='y':
            b_,oy_=master(l); ylow=oy_-np.nonzero(b_.any(1))[0].max()
            longs=[t for t in strokes if t[2]<-60]; shorts=[t for t in strokes if t[2]>=-60]
            for p0,p1,low in longs: m|=C.diag(p0,p1,ylow,525)
            if longs:
                p0,p1,_=longs[0]; a=np.array([C.X(p0[0]),p0[1]]); bq=np.array([C.X(p1[0]),p1[1]]); d=unit(bq-a)
                xedge=lambda yy: a[0]+(yy-a[1])*d[0]/d[1]+ (T/2)/abs(d[1])
                for p0_,p1_,_ in shorts: m|=C.diag(p0_,p1_,0,525)&(C.xx<=xedge(C.yu))
            else:
                for p0,p1,low in strokes: m|=C.diag(p0,p1,0,525)
        elif l=='k':
            arm=[t for t in strokes if (t[1][1]-t[0][1])*(t[1][0]-t[0][0])>0]
            leg=[t for t in strokes if (t[1][1]-t[0][1])*(t[1][0]-t[0][0])<0]
            for p0,p1,_ in arm: m|=C.diag(p0,p1,0,525)
            if arm:
                p0,p1,_=arm[0]; a=np.array([C.X(p0[0]),p0[1]]); bq=np.array([C.X(p1[0]),p1[1]]); d=unit(bq-a)
                yarm=lambda xx: a[1]+(xx-a[0])*d[1]/d[0]
                for p0_,p1_,_ in leg: m|=C.diag(p0_,p1_,0,525)&(C.yu<=yarm(C.xx))
        elif l=='w':
            # slimmer where strokes meet (bottoms, centre top) so the joins don't spread into serifs
            tops=sorted(strokes,key=lambda t:C.X(t[0][0])+(525-t[0][1])*(t[1][0]-t[0][0])/(t[1][1]-t[0][1]))
            outer={id(tops[0]),id(tops[-1])}
            dl=0.16*T; yj=230
            for t_ in strokes:
                p0,p1,_=t_
                St=C.diag(p0,p1,0,525)
                a=np.array([C.X(p0[0]),p0[1]]); bq=np.array([C.X(p1[0]),p1[1]]); d=unit(bq-a)
                xc=lambda y: a[0]+(y-a[1])*d[0]/d[1]
                hw=(T/2)/abs(d[1])                             # half-width measured horizontally
                downright=(bq[0]-a[0])*(bq[1]-a[1])<0
                # bottom join: pull the outer edge in (left edge of '\\', right edge of '/')
                ramp=np.clip(1-C.yu/yj,0,1)*dl
                if not downright: cut=(C.xx<xc(C.yu)-hw+ramp)&(C.xx>xc(C.yu)-hw-15)
                else:             cut=(C.xx>xc(C.yu)+hw-ramp)&(C.xx<xc(C.yu)+hw+15)
                St&=~(cut&(C.yu<yj))
                if id(t_) not in outer:                          # centre top join of the inner strokes
                    ramp=np.clip(1-(525-C.yu)/yj,0,1)*dl
                    if not downright: cut=(C.xx>xc(C.yu)+hw-ramp)&(C.xx<xc(C.yu)+hw+15)
                    else:             cut=(C.xx<xc(C.yu)-hw+ramp)&(C.xx>xc(C.yu)-hw-15)
                    St&=~(cut&(C.yu>525-yj))
                m|=St
        else:
          for p0,p1,low in strokes:
            m|=C.diag(p0,p1,0,525)
        if l=='k':
            st_=runs_at(l,150)[0]
            m&=(C.xx>=C.X(sc(st_)))
            m|=C.stem(sc(st_),0,ASC,top=term_line(l,'top'))
    return m.astype(np.uint8),C.pad,C.noy
def fit_anchored(P,anchor,side,guess):
    """outer-edge points of a ring whose centre-line passes through the stem centre `anchor`
    (side=+1: anchor is the ring's left extreme, -1: right extreme). Returns centre-line ring."""
    x,y=P[:,0],P[:,1]; h,v=70.0,65.0
    def res(p):
        a,cy,b=p; cx=anchor+side*a
        return (np.sqrt(((x-cx)/(a+h))**2+((y-cy)/(b+v))**2)-1)*min(a+h,b+v)
    s=least_squares(res,guess,bounds=([20,-600,20],[600,900,600]))
    a,cy,b=s.x; return dict(cx=anchor+side*a,cy=cy,a=a,b=b,err=np.abs(s.fun).max())
def hook(l,stem_run,where):
    """quarter/half-ring hook from the stem: where='bottom-left'|'bottom-right'|'top-right'"""
    b,oy=master(l); xs_c=sc(stem_run)
    ys,xs=np.nonzero(b); yu=oy-ys
    if where.startswith('bottom'):
        if where=='bottom-left':
            sel=(xs>xs.min()+60)&(xs<xs_c)
            k=np.argmin(np.where(sel,yu,1e9))
        else: k=np.argmin(yu)
        x_ext,y_ext=xs[k],yu[k]
        # where the straight stem ends: lowest point of the stem's inner edge line
        edge=[L for L in ML[l] if abs(L['u'][1])>0.99 and abs(L['c'][0]-(stem_run[0] if where=='bottom-left' else stem_run[1]))<15]
        cy=min(oy-(L['c']+L['u']*t)[1] for L in edge for t in (L['t0'],L['t1'])) if edge else y_ext+200
        a=abs(xs_c-x_ext)
        if where=='bottom-left': cy=y_ext+65+a*0.95      # a near-circular turn into the hook
        bb=cy-(y_ext+65)
        cx=xs_c-a if where=='bottom-left' else xs_c+a
    else:
        k=np.argmax(yu); x_ext,y_ext=xs[k],yu[k]
        edge=[L for L in ML[l] if abs(L['u'][1])>0.99 and abs(L['c'][0]-stem_run[1])<15]
        cy=max(oy-(L['c']+L['u']*t)[1] for L in edge for t in (L['t0'],L['t1'])) if edge else y_ext-150
        a=abs(x_ext-xs_c); bb=(y_ext-65)-cy; cx=xs_c+a
    return dict(cx=cx,cy=cy,a=max(a,20),b=max(bb,20))
from skimage.morphology import skeletonize
def centre_pts(l,xr,yr):
    b,oy=master(l)
    dt=cv2.distanceTransform(b,cv2.DIST_L2,5); sk=skeletonize(b>0)
    ys,xs=np.nonzero(sk&(dt>40)); yu=oy-ys
    k=(xs>=xr[0])&(xs<=xr[1])&(yu>=yr[0])&(yu<=yr[1])
    return np.column_stack([xs[k],yu[k]]).astype(float)
def fit_centre(P,anchor,side,guess):
    x,y=P[:,0],P[:,1]
    def res(p):
        a,cy,b=p; cx=anchor+side*a
        return (np.sqrt(((x-cx)/a)**2+((y-cy)/b)**2)-1)*min(a,b)
    s=least_squares(res,guess,bounds=([15,-700,15],[700,1000,700]))
    a,cy,b=s.x; return dict(cx=anchor+side*a,cy=cy,a=a,b=b,err=np.abs(s.fun).max())
def build_hooked(l,f,T):
    C=Canvas(f,T); m=C.mask()>0; b,oy=master(l); x0,x1=ink_x(l)
    if l in 'gj':
        s=runs_at(l,300)[-1]; xs_=sc(s)
        if l=='j':
            ys_,xs2=np.nonzero(b); yb_=(oy-ys_).min()+65
            a_=(xs_-x0)*0.75; rf=dict(cx=xs_-a_,cy=yb_+a_*0.9,a=a_,b=a_*0.9)
        else:
            ys_,xs2=np.nonzero(b); yu2=oy-ys_; yb_=yu2.min()+65
            tl_=term_line(l,'bottom'); a_=(xs_-tl_['p'][0])*0.62
            rf=dict(cx=xs_-a_,cy=yb_+a_*0.9,a=a_,b=a_*0.9)
        r,g=C.ring(rf); keep=r&(C.yu<=g['cy'])&(C.xx<=C.X(xs_)+T/2)
        if l=='g':
            tl_=term_line(l,'bottom')
            P_=np.array([C.X(tl_['p'][0]),tl_['p'][1]]); u_=tl_['u']; n_=np.array([-u_[1],u_[0]])
            sd=(C.xx-P_[0])*n_[0]+(C.yu-P_[1])*n_[1]
            ref=(C.X(xs_)-P_[0])*n_[0]+(g['cy']-g['bc']-P_[1])*n_[1]
            keep&=(sd*np.sign(ref)>0)|(C.xx>g['cx'])        # keep the stem's side of the terminal cut
            keep&=~((C.xx<g['cx'])&(C.yu>g['cy']-0.45*g['bc']))   # the tail stays low; it doesn't curl back up
        else:
            keep&=(C.xx>=C.X(x0))
        m|=keep; m|=C.stem(xs_,g['cy'],525)
        if l=='j': m|=C.dot('j')
        if l=='g':
            import geomc2
            po,pi,_,_=geomc2.G['q']; c=geomc2.cl_params(po,pi)
            qb,qoy=master('q'); qs=runs_at('q',300)[-1]
            dx=xs_-sc(qs)                         # line the q bowl up with this stem
            L_=C.X(c['Lc']+dx); R_=C.X(c['Rc']+dx); cxq=(L_+R_)/2; acq=(R_-L_)/2
            bcq=c['bo']-T/2; cyq=c['cy'] if False else (qoy-c['cy'])
            img=np.zeros((C.H,C.W),np.uint8)
            cv2.ellipse(img,(int(round(cxq*16)),int(round(C.R(cyq)*16))),(int(round(acq*16)),int(round(bcq*16))),0,0,360,1,thickness=max(1,int(round(T))),shift=4)
            m|=(img>0)&(C.xx<=C.X(xs_))&(C.yu>=-5)
    elif l=='t':
        s=max(runs_at(l,250),key=lambda r:r[1]-r[0]); xs_=sc(s)
        ys_,xs2=np.nonzero(b); yu2=oy-ys_
        xt=xs2[yu2<150].max(); yb_=yu2.min()+65
        a_=(xt-xs_)*0.62; rf=dict(cx=xs_+a_,cy=yb_+a_*1.0,a=a_,b=a_*1.0)
        r,g=C.ring(rf)
        cutr=C.right_of(term_line(l,'bottom'))&(C.xx>=g['cx'])          # the terminal only trims the far end
        keep=r&(C.yu<=g['cy'])&(C.xx>=C.X(xs_)-T/2)&~cutr
        m|=keep; m|=C.stem(xs_,g['cy'],ASC_,top=term_line(l,'top'))
        br=max(runs_at(l,460),key=lambda r:r[1]-r[0]); th=131*T/TM
        m|=(C.yu>=527-th)&(C.yu<=527)&(C.xx>=C.X(br[0]))&(C.xx<=C.X(br[1]))
    elif l=='f':
        s=max(runs_at(l,250),key=lambda r:r[1]-r[0]); xs_=sc(s)
        P=centre_pts(l,(xs_-10,x1),(500,720))
        rf=fit_centre(P,xs_,+1,[(x1-xs_)/2,560,100])
        ys_,xs2=np.nonzero(b); yu2=oy-ys_; endcol=(xs2>=x1-35)&(yu2>540)
        y_end=yu2[endcol].min() if endcol.any() else 600
        r,g=C.ring(rf); keep=r&(C.yu>=g['cy'])&(C.xx>=C.X(xs_)-T/2)&(C.xx<=C.X(x1))&((C.yu>=y_end)|(C.xx<=C.X(xs_)+T/2))
        m|=keep; m|=C.stem(xs_,0,g['cy'])
        br=max(runs_at(l,455),key=lambda r:r[1]-r[0]); th=141*T/TM
        m|=(C.yu>=525-th)&(C.yu<=525)&(C.xx>=C.X(br[0]))&(C.xx<=C.X(br[1]))
    elif l=='a':
        s=runs_at(l,150)[-1]; xs_=sc(s)
        P=centre_pts(l,(x0,xs_+20),(330,600))
        top_c=P[:,1].max(); x_top=np.median(P[P[:,1]>=top_c-4,0])
        a_=xs_-x_top; b_=a_*0.85; cy_=top_c-b_
        ra=dict(cx=x_top,cy=cy_,a=a_,b=b_)
        r2,g2=C.ring(ra); tl_=term_line(l,'top')
        lowcut=(C.xx<g2['cx'])&(C.yu<tl_['p'][1]-25)          # the arc's left end stops at the terminal
        m|=r2&(C.yu>=g2['cy'])&(C.xx<=C.X(xs_)+T/2+1)&C.right_of(tl_)&~lowcut
        m|=(C.xx>=C.X(xs_)-T/2)&(C.xx<C.X(xs_)+T/2)&(C.yu<=g2['cy'])&(C.yu>=0)   # stem stops where the arc takes over
        Pb=centre_pts(l,(x0,s[0]),(-20,330))
        x,y=Pb[:,0],Pb[:,1]
        def res(p):
            cx,cy,a,bb=p; return (np.sqrt(((x-cx)/a)**2+((y-cy)/bb)**2)-1)*min(a,bb)
        cx,cy,a,bb=least_squares(res,[(x0+xs_)/2,150,(xs_-x0)/2,90]).x
        rb=dict(cx=cx,cy=cy,a=a,b=bb); r,_=C.ring(rb,0.9); m|=r&(C.xx<=C.X(xs_))
    m=cv2.morphologyEx(m.astype(np.uint8),cv2.MORPH_CLOSE,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(7,7)))
    return m,C.pad,C.noy
ASC_=700
