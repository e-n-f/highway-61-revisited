"""Two-ellipse rings (1966 counters) shared between letters, per the curvature findings."""
import numpy as np,cv2,pickle
import geomc3 as g3
from geomc3 import Canvas,master,runs_at,sc,ink_x,term_line,TM
G=pickle.load(open('geomfit.pkl','rb'))
import os
def rp_of(key):
    po,pi,oy,_=G[key]
    return dict(Lc=((po[0]-po[2])+(pi[0]-pi[2]))/2, Rc=((po[0]+po[2])+(pi[0]+pi[2]))/2,
                cy=oy-po[1], bo=po[3], top=(pi[1]-pi[3])-(po[1]-po[3]), bot=(po[1]+po[3])-(pi[1]+pi[3]))
OFF=pickle.load(open('ring_off.pkl','rb')) if os.path.exists('ring_off.pkl') else {}
def ring2(C,rp,dx=0.0,dy=0.0,tag=None):
    """outer and inner ellipse (the 1966 construction) - even sides, wider counter; dx,dy move it (master units)"""
    T=C.T; k=T/TM
    if tag is not None and tag in OFF: dx+=OFF[tag][0]; dy+=OFF[tag][1]
    L=C.X(rp['Lc']+dx); R=C.X(rp['Rc']+dx); cx=(L+R)/2
    ao=(R-L)/2+T/2; ai=max(2,(R-L)/2-T/2)
    cyu=rp['cy']+dy; bo=rp['bo']; top=rp['top']*k; bot=rp['bot']*k
    bi=bo-(top+bot)/2; cyi=cyu-(top-bot)/2
    o=C.mask(); i=C.mask()
    cv2.ellipse(o,(int(round(cx*16)),int(round(C.R(cyu)*16))),(int(round(ao*16)),int(round(bo*16))),0,0,360,1,-1,shift=4)
    cv2.ellipse(i,(int(round(cx*16)),int(round(C.R(cyi)*16))),(int(round(ai*16)),int(round(bi*16))),0,0,360,1,-1,shift=4)
    return (o>0)&(i==0), dict(cx=cx,cy=cyu,ao=ao,ai=ai,bo=bo,bi=bi,cyi=cyi,L=L,R=R)
def stemrun(l,y): return max(runs_at(l,y),key=lambda r:r[1]-r[0])
def build(l,f,T):
    C=Canvas(f,T); m=C.mask()>0; b,oy=master(l); x0,x1=ink_x(l); ASC=700
    P=rp_of('p'); pst=sc(runs_at('p',300)[0])          # p's stem centre
    if l in 'bdpq':
        rp=rp_of(l); r,gg=ring2(C,rp); m|=r
        st=runs_at(l,300)[0] if l in 'bp' else runs_at(l,300)[-1]; xs=sc(st)
        up=l in 'bd'; TL=[L for L in g3.term_line.__globals__['ML'][l] if 3<np.degrees(np.arctan2(abs(L['u'][1]),abs(L['u'][0])))<87 and abs(L['u'][0])>0.5]
        Lt=(min if up else max)(TL,key=lambda L_:L_['c'][1]); tl=dict(p=np.array([Lt['c'][0],oy-Lt['c'][1]]),u=np.array([Lt['u'][0],-Lt['u'][1]]))
        if up: m|=C.stem(xs,0,ASC,top=tl)
        else:  m|=C.stem(xs,-400,525,bot=tl)
        near=(C.xx>=C.X(xs)-T*1.2)&(C.xx<C.X(xs)+T*1.2)
        if up: m&=~(near&(C.yu<0))
        else:  m&=~(near&(C.yu>525))
    elif l=='r':
        xs=sc(runs_at('r',150)[0]); r,gg=ring2(C,P,dx=xs-pst,tag='r')
        tl=term_line('r','right')
        m|=r&(C.yu>=gg['cy'])&(C.xx>=C.X(xs))&~C.right_of(tl)
        m|=C.stem(xs,0,525); m&=~(C.yu>525)|((C.xx>=C.X(xs)-T/2)&(C.xx<C.X(xs)+T/2)&False)
        m|=r&(C.yu>=gg['cy'])&(C.yu<=530)&(C.xx>=C.X(xs))&~C.right_of(tl)
    elif l in 'nhm':
        st=runs_at(l,150)
        m|=C.stem(sc(st[0]),0,525 if l!='h' else ASC,top=term_line(l,'top') if l=='h' else None)
        for rs in st[1:]:
            xr=sc(rs); r,gg=ring2(C,P,dx=xr-P['Rc'],tag=l+str(st.index(rs)))
            m|=r&(C.yu>=gg['cy'])&(C.xx>=C.X(sc(st[0])))
            m|=(C.xx>=C.X(xr)-T/2)&(C.xx<C.X(xr)+T/2)&(C.yu>=0)&(C.yu<=gg['cy'])
    elif l=='t':
        s=stemrun('t',250); xs=sc(s); ys_,xs2=np.nonzero(b); ybot=(oy-ys_).min()
        dy=(ybot+P['bo'])-P['cy']; r,gg=ring2(C,P,dx=xs-P['Lc'],dy=dy,tag='t')
        tl=term_line('t','bottom'); cut=C.right_of(tl)&(C.xx>=gg['cx'])
        m|=r&(C.yu<=gg['cy'])&(C.xx>=C.X(xs)-T/2)&~cut
        m|=C.stem(xs,gg['cy'],ASC,top=term_line('t','top'))
        br=max(runs_at('t',460),key=lambda r_:r_[1]-r_[0]); th=131*T/TM
        m|=(C.yu>=527-th)&(C.yu<=527)&(C.xx>=C.X(br[0]))&(C.xx<=C.X(br[1]))
    elif l=='j':
        s=stemrun('j',300); xs=sc(s); ys_,xs2=np.nonzero(b); ybot=(oy-ys_).min()
        dy=(ybot+P['bo'])-P['cy']; r,gg=ring2(C,P,dx=xs-P['Rc'],dy=dy,tag='j')
        m|=r&(C.yu<=gg['cy'])&(C.xx<=C.X(xs)+T/2)&(C.xx>=C.X(x0))
        m|=C.stem(xs,gg['cy'],525); m|=C.dot('j')
    elif l in 'af':
        O=rp_of('o'); ys_,xs2=np.nonzero(b); yu_=oy-ys_
        if l=='a':
            s=runs_at('a',150)[-1]; xs=sc(s); top=yu_.max()
            dy=(top-O['bo'])-O['cy']; r,gg=ring2(C,O,dx=xs-O['Rc'],dy=dy,tag='a')
            tl=term_line('a','top'); low=(C.xx<gg['cx'])&(C.yu<tl['p'][1]-25)
            m|=r&(C.yu>=gg['cy'])&(C.xx<=C.X(xs)+T/2)&C.right_of(tl)&~low
            m|=C.stem(xs,0,gg['cy'])
            # bowl as before
            img,pad,noy=g3.build_hooked('a',f,T); bowl_only=img>0
            m|=bowl_only&(C.yu<gg['cy']-T*0.4)&(C.xx<C.X(xs)-T/2+1)
        else:
            s=stemrun('f',250); xs=sc(s); top=yu_.max()
            dy=(top-O['bo'])-O['cy']; r,gg=ring2(C,O,dx=xs-O['Lc'],dy=dy,tag='f')
            endcol=(xs2>=x1-35)&(yu_>540); y_end=yu_[endcol].min() if endcol.any() else 600
            m|=r&(C.yu>=gg['cy'])&(C.xx>=C.X(xs)-T/2)&(C.xx<=C.X(x1))&((C.yu>=y_end)|(C.xx<=C.X(xs)+T/2))
            m|=C.stem(xs,0,gg['cy'])
            br=max(runs_at('f',455),key=lambda r_:r_[1]-r_[0]); th=141*T/TM
            m|=(C.yu>=525-th)&(C.yu<=525)&(C.xx>=C.X(br[0]))&(C.xx<=C.X(br[1]))
    elif l=='g':
        img3,pad3,noy3=g3.build_hooked('g',f,T); tail=(img3>0)&(C.yu<-5)
        s=runs_at('g',300)[-1]; xs=sc(s); qs=sc(runs_at('q',300)[-1])
        r,gg=ring2(C,rp_of('q'),dx=xs-qs)
        m|=tail|(r&(C.yu>=-5)&(C.xx<=C.X(xs)+T/2))|C.stem(xs,-60,525)
    elif l=='c':
        O=rp_of('o'); xs_=np.nonzero(b.any(0))[0]; po=G['o'][0]
        dx=xs_[0]-(po[0]-po[2]); r,gg=ring2(C,O,dx=dx,tag='c')
        m|=r
        ML=g3.ML
        off=OFF.get('c',(0,0))
        mcx=(O['Lc']+O['Rc'])/2+dx+off[0]; ma=(O['Rc']-O['Lc'])/2; mcy=O['cy']+off[1]; mb=O['bo']-O['top']/2
        nb=gg['bo']-(O['top']*T/TM)/2; na=(gg['R']-gg['L'])/2
        lines=[]
        for L in ML['c']:
            u=L['u']; ang=np.degrees(np.arctan2(abs(u[1]),abs(u[0])))
            if not (3<ang<87): continue
            px,py=L['c'][0],oy-L['c'][1]
            th=np.arctan2((py-mcy)/mb,(px-mcx)/ma)                 # where on the ring the terminal sits
            p=np.array([gg['cx']+na*np.cos(th), gg['cy']+nb*np.sin(th)])
            slope=(-u[1])/u[0]
            lines.append((p,slope))
        if len(lines)>=2:
            lines.sort(key=lambda t:t[0][1])
            (pb,sb),(pt,st_)=lines[0],lines[-1]
            below_top=C.yu< pt[1]+st_*(C.xx-pt[0])
            above_bot=C.yu> pb[1]+sb*(C.xx-pb[0])
            m&=~((C.xx>gg['cx'])&below_top&above_bot)          # the opening between the two terminal cuts
    return m.astype(np.uint8),C.pad,C.noy
