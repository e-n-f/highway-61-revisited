"""Elliptical constructions marked up by the user: c, n, m, r, u (inner ellipse + outward pen profile)."""
import numpy as np,cv2,pickle
import geomc3 as g3, geomc4 as g4, geomc5 as g5
from geomc3 import Canvas,master,runs_at,sc,ink_x,TM,edge_pts,term_line
from geomfit import fit_super
def fit2(P): p,e=fit_super(P,2.0); return p,e
def edge(l,xr,yr,kind):
    P=edge_pts(l,xr,yr,kind); b,oy=master(l)
    return np.column_stack([P[:,0],oy-P[:,1]])        # back to rows for fit_super/profile convention
# ---------- fits (master frame, rows) ----------
def fit_arch(l,c0,c1,s0,s1,xmax=None):
    """outer ellipse over an arch between stems (c0..c1 centres, s0 right edge of left stem, s1 left edge of right stem)"""
    b,oy=master(l)
    Po=edge(l,(s0+5,xmax if xmax else s1+150),(330,600),'top')
    Pi=edge(l,(s0+4,s1-4),(300,600),'bottom')
    po,eo=fit2(Po); pi,ei=fit2(Pi)
    return po,pi,oy,(eo,ei)
def fit_u():
    b,oy=master('u'); st=runs_at('u',450); s0,s1=st[0][1],st[1][0]
    Po=np.vstack([edge('u',(st[0][0]-5,s1),(-30,200),'bottom'),edge('u',(st[0][0]-5,st[0][0]+120),(-30,200),'left')])
    Pi=edge('u',(s0+4,s1-4),(-30,260),'top')
    po,eo=fit2(Po); pi,ei=fit2(Pi); return po,pi,oy,(eo,ei)
def fit_c():
    b,oy=master('c'); ys,xs=np.nonzero(b)
    cs,hh=cv2.findContours(b.copy(),cv2.RETR_CCOMP,cv2.CHAIN_APPROX_NONE)
    C_=max([c[:,0,:] for c in cs],key=len).astype(float)
    ctr=np.array([(xs.min()+xs.max())/2,(ys.min()+ys.max())/2])
    d=C_-ctr; ang=np.degrees(np.arctan2(-d[:,1],d[:,0]))%360; rad=np.hypot(d[:,0],d[:,1])
    O=[];I=[]
    for a0 in range(70,290,4):                 # away from the opening on the right
        k=(ang>=a0)&(ang<a0+4)
        if k.sum()<2: continue
        idx=np.nonzero(k)[0]; O.append(C_[idx[np.argmax(rad[idx])]]); I.append(C_[idx[np.argmin(rad[idx])]])
    po,_=fit2(np.array(O)); pi,_=fit2(np.array(I))
    return po,pi,oy
FIT={}
st=runs_at('n',150); FIT['n']=fit_arch('n',sc(st[0]),sc(st[1]),st[0][1],st[1][0])
st=runs_at('m',150); FIT['m1']=fit_arch('m',sc(st[0]),sc(st[1]),st[0][1],st[1][0],xmax=sc(st[1])); FIT['m2']=fit_arch('m',sc(st[1]),sc(st[2]),st[1][1],st[2][0])
st=runs_at('h',150); FIT['h']=fit_arch('h',sc(st[0]),sc(st[1]),st[0][1],st[1][0])
xr=runs_at('r',150)[0]; b_,oy_=master('r'); x1r=ink_x('r')[1]
_tl=term_line('r','right'); xcut=_tl['p'][0]-40
Po=edge('r',(xr[1]+5,xcut+30),(380,600),'top'); Pi=edge('r',(xr[1]+3,xcut),(300,560),'bottom')
FIT['r']=(fit2(Po)[0],fit2(Pi)[0],oy_,None)
FIT['u']=fit_u()
FIT['c']=fit_c()+((None,),)
PROF={k:g5.profile(v[0],v[1],v[2]) for k,v in FIT.items()}
def ring(C,key,inner_x=None,dx=0.0,side_theta=np.pi,junction=None):
    """inner ellipse + pen profile. inner_x=(xl,xr) in OUTPUT coords fixes the counter's width;
    otherwise the side at side_theta keeps its centre-line scaled by f."""
    pr=PROF[key]; th=pr['th']
    k=C.T/TM if junction is None else C.T/pr['t'][np.argmin(np.abs(th-junction))]   # leg joins flush
    t=pr['t']*k
    i90=np.argmin(np.abs(th-np.pi/2)); i270=np.argmin(np.abs(th-3*np.pi/2))
    if inner_x is not None:
        cx=(inner_x[0]+inner_x[1])/2; ai=max(3.0,(inner_x[1]-inner_x[0])/2)
    else:
        isd=np.argmin(np.abs(th-side_theta)); tr=pr['t'][isd]
        cx=C.X(pr['ci'][0]+dx); ai=max(3.0,C.f*(pr['ai']+tr/2)-tr*k/2)
    top=pr['otop']-pr['t'][i90]*k; bot=pr['obot']+pr['t'][i270]*k
    if pr['t'][i270]<0.15*pr['t'][i90]:                    # inner touches the outer at the bottom (n, m, r)
        bot=pr['ci'][1]-pr['bi']
    if pr['t'][i90]<0.15*pr['t'][i270]:                    # ... or at the top (u)
        top=pr['ci'][1]+pr['bi']
    cy=(top+bot)/2; bi=(top-bot)/2
    P=np.column_stack([cx+ai*np.cos(th),cy+bi*np.sin(th)])
    N=np.column_stack([np.cos(th)/ai,np.sin(th)/bi]); N/=np.linalg.norm(N,axis=1)[:,None]
    Q=P+N*t[:,None]
    out=C.mask(); inn=C.mask()
    toimg=lambda Z: np.round(np.column_stack([Z[:,0],C.noy-Z[:,1]])*16).astype(np.int32)
    cv2.fillPoly(out,[toimg(Q)],1,shift=4); cv2.fillPoly(inn,[toimg(P)],1,shift=4)
    oc_y=Q[:,1]; 
    return (out>0)&(inn==0), dict(cx=cx,cy=cy,ai=ai,bi=bi,Q=Q,P=P,ocy=(Q[:,1].max()+Q[:,1].min())/2,
                                    oxr=Q[:,0].max(),oxl=Q[:,0].min())
def build(l,f,T):
    C=Canvas(f,T); m=C.mask()>0; b,oy=master(l); ASC=700
    if l in 'nhm':
        st=runs_at(l,150)
        m|=C.stem(sc(st[0]),0,525 if l!='h' else ASC,top=term_line(l,'top') if l=='h' else None)
        keys={'n':['n'],'h':['h'],'m':['n','n']}[l]
        for i,key in enumerate(keys):
            a,bb=st[i],st[i+1]
            xl=C.X(sc(a))+T/2; xr=C.X(sc(bb))-T/2
            r,gg=ring(C,key,inner_x=(xl,xr),junction=0.0)
            m|=r&(C.xx>=C.X(sc(a)))&(C.yu>=gg['cy'])&(C.xx<xr+T)     # arch doesn't bulge past the leg
            # right leg: stem from the baseline up into the arch
            m|=(C.xx>=xr)&(C.xx<xr+T)&(C.yu>=0)&(C.yu<=gg['cy']+1)
    elif l=='u':
        st=runs_at('u',450)
        xl=C.X(sc(st[0]))+T/2; xr=C.X(sc(st[1]))-T/2
        r,gg=ring(C,'u',inner_x=(xl,xr),junction=np.pi)
        m|=r&(C.xx<=C.X(sc(st[1])))&(C.yu<=gg['cy'])&(C.xx>=xl-T)
        m|=(C.xx>=xl-T)&(C.xx<xl)&(C.yu>=gg['cy']-1)&(C.yu<=525)
        m|=C.stem(sc(st[1]),0,525)
    elif l=='r':
        s=runs_at('r',150)[0]; xs=sc(s)
        stn=runs_at('n',150); span=sc(stn[1])-sc(stn[0])        # same construction as the n's arch
        xl=C.X(xs)+T/2
        r,gg=ring(C,'n',inner_x=(xl,xl+f*span-T),junction=0.0)
        tl=term_line('r','right')
        m|=r&(C.xx>=C.X(xs))&~C.right_of(tl)&(C.yu<=530)&(C.yu>=gg['cy'])
        m|=C.stem(xs,0,525)
    elif l=='c':
        r,gg=ring(C,'c',side_theta=np.pi,junction=np.pi); m|=r
        pr=PROF['c']; po,pi,oyc,_=FIT['c']
        mcx=pr['ci'][0]; mcy=pr['ci'][1]
        lines=[]
        for L in g3.ML['c']:
            u=L['u']; ang=np.degrees(np.arctan2(abs(u[1]),abs(u[0])))
            if not (3<ang<87): continue
            px,py=L['c'][0],oy-L['c'][1]
            th=np.arctan2((py-mcy)/pr['bi'],(px-mcx)/pr['ai'])
            p=np.array([gg['cx']+gg['ai']*np.cos(th), gg['cy']+gg['bi']*np.sin(th)])
            lines.append((p,(-u[1])/u[0]))
        if len(lines)>=2:
            lines.sort(key=lambda t_:t_[0][1]); (pb,sb),(pt,st_)=lines[0],lines[-1]
            m&=~((C.xx>gg['cx'])&(C.yu<pt[1]+st_*(C.xx-pt[0]))&(C.yu>pb[1]+sb*(C.xx-pb[0])))
    return m.astype(np.uint8),C.pad,C.noy
