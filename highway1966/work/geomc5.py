"""Bowls as an inner ellipse plus an outward pen whose width varies around it (thickness
profile measured from the 1966 drawing and kept the same in every width, scaled by weight)."""
import numpy as np,cv2,pickle
import geomc3 as g3, geomc4 as g4
from geomc3 import Canvas,master,runs_at,sc,ink_x,TM,edge_pts,fit_ellipse
from geomfit import fit_super
G=pickle.load(open('geomfit.pkl','rb'))
def fit_a_bowl():
    b,oy=master('a'); x0,x1=ink_x('a'); s=runs_at('a',150)[-1]
    cs,hh=cv2.findContours(b.copy(),cv2.RETR_CCOMP,cv2.CHAIN_APPROX_NONE)
    hole=max([c[:,0,:] for i,c in enumerate(cs) if hh[0][i][3]!=-1],key=len)
    Hp=hole[hole[:,0]<s[0]-8].astype(float)
    outer=max([c[:,0,:] for i,c in enumerate(cs) if hh[0][i][3]==-1],key=len).astype(float)
    yu=oy-outer[:,1]
    Po=outer[(outer[:,0]<s[0]-8)&(yu<330)&((yu<250)|(outer[:,0]<x0+120))]
    po,_=fit_super(Po,2.0); pi,_=fit_super(Hp,2.0)
    return po,pi,oy
BOWL={}
for k in 'bdpq': BOWL[k]=(G[k][0],G[k][1],G[k][2])
BOWL['a']=fit_a_bowl()
def profile(po,pi,oy,n=720):
    """inner ellipse (x, y-up) and the thickness along its outward normal to the outer ellipse"""
    ci=np.array([pi[0],oy-pi[1]]); ai,bi=pi[2],pi[3]; co=np.array([po[0],oy-po[1]]); ao,bo=po[2],po[3]
    th=np.linspace(0,2*np.pi,n,endpoint=False)
    P=np.column_stack([ci[0]+ai*np.cos(th),ci[1]+bi*np.sin(th)])
    N=np.column_stack([np.cos(th)/ai,np.sin(th)/bi]); N/=np.linalg.norm(N,axis=1)[:,None]
    # ray P+t*N meets outer ellipse
    d=P-co; A=(N[:,0]/ao)**2+(N[:,1]/bo)**2; B=2*(d[:,0]*N[:,0]/ao**2+d[:,1]*N[:,1]/bo**2); Cc=(d[:,0]/ao)**2+(d[:,1]/bo)**2-1
    t=(-B+np.sqrt(np.maximum(B*B-4*A*Cc,0)))/(2*A)
    return dict(ci=ci,ai=ai,bi=bi,th=th,t=t,otop=co[1]+bo,obot=co[1]-bo)
PROF={k:profile(*v) for k,v in BOWL.items()}
def bowl(C,key,dx=0.0,right_to=None):
    """draw the bowl of `key` in canvas C (x mapped, weight scaled), shifted by dx master units"""
    pr=PROF[key]; k=C.T/TM; th=pr['th']; t=pr['t']*k
    i0=np.argmin(np.abs(th-0)); i90=np.argmin(np.abs(th-np.pi/2)); i270=np.argmin(np.abs(th-3*np.pi/2))
    tr=pr['t'][i0]                                       # side thickness in the drawing
    cx=C.X(pr['ci'][0]+dx)
    ai=max(3.0,C.f*(pr['ai']+tr/2)-tr*k/2)               # centre-line of the side narrows with the width
    cx0,ai0=cx,ai
    if right_to is not None:                                # widen only the counter so the stem cuts it flat
        left=cx-ai; ai=max(ai,(right_to-left)/2); cx=left+ai
    top=pr['otop']-pr['t'][i90]*k; bot=pr['obot']+pr['t'][i270]*k   # outer top/bottom stay where drawn
    cy=(top+bot)/2; bi=(top-bot)/2
    P0=np.column_stack([cx0+ai0*np.cos(th),cy+bi*np.sin(th)])       # outline drawn from the original counter
    N=np.column_stack([np.cos(th)/ai0,np.sin(th)/bi]); N/=np.linalg.norm(N,axis=1)[:,None]
    Q=P0+N*t[:,None]
    P=np.column_stack([cx+ai*np.cos(th),cy+bi*np.sin(th)])
    out=C.mask(); inn=C.mask()
    toimg=lambda Z: np.round(np.column_stack([Z[:,0],C.noy-Z[:,1]])*16).astype(np.int32)
    cv2.fillPoly(out,[toimg(Q)],1,shift=4); cv2.fillPoly(inn,[toimg(P)],1,shift=4)
    return (out>0)&(inn==0), dict(cx=cx,cy=cy,ai=ai,bi=bi)
def build(l,f,T):
    C=Canvas(f,T); b,oy=master(l)
    if l in 'bdpq':
        img,pad,noy=g4.build(l,f,T); stemetc=img>0
        ring,_=g4.ring2(C,g4.rp_of(l))                    # remove the old bowl, keep stem/cuts
        st=runs_at(l,300)[0] if l in 'bp' else runs_at(l,300)[-1]; xs=sc(st)
        stem=(C.xx>=C.X(xs)-T/2)&(C.xx<C.X(xs)+T/2)
        keep=stemetc&stem
        bw,_=bowl(C,l); up=l in 'bd'
        near=(C.xx>=C.X(xs)-T*1.2)&(C.xx<C.X(xs)+T*1.2)
        bw&=~(near&((C.yu<0) if up else (C.yu>525)))
        return (keep|bw).astype(np.uint8),C.pad,C.noy
    if l=='g':
        img,pad,noy=g4.build('g',f,T); s=runs_at('g',300)[-1]; xs=sc(s); qs=sc(runs_at('q',300)[-1])
        stem=(C.xx>=C.X(xs)-T/2)&(C.xx<C.X(xs)+T/2)
        keep=(img>0)&(stem|(C.yu<-5))
        bw,_=bowl(C,'q',dx=xs-qs); bw&=(C.yu>=-5)
        return (keep|bw).astype(np.uint8),C.pad,C.noy
    if l=='r':
        img,pad,noy=g4.build('r',f,T); xs=sc(runs_at('r',150)[0]); ps=sc(runs_at('p',300)[0])
        stem=(C.xx>=C.X(xs)-T/2)&(C.xx<C.X(xs)+T/2)
        bw,gg=bowl(C,'p',dx=xs-ps); tl=g3.term_line('r','right')
        arc=bw&(C.yu>=gg['cy'])&(C.xx>=C.X(xs))&~C.right_of(tl)&(C.yu<=530)
        return (((img>0)&stem)|arc).astype(np.uint8),C.pad,C.noy
    if l=='a':
        img,pad,noy=g3.build_hooked('a',f,T); s=runs_at('a',150)[-1]; xs=sc(s)
        stem=(C.xx>=C.X(xs)-T/2)&(C.xx<C.X(xs)+T/2)
        bw,gg=bowl(C,'a')
        # keep the arc and stem from before; replace the bowl
        old_bowl_zone=(C.yu<gg['cy']+ (PROF['a']['otop']-PROF['a']['ci'][1]))&(C.xx<C.X(xs)-T/2)
        keep=(img>0)&~(old_bowl_zone&(C.yu<300))
        return (keep|bw&(C.xx<=C.X(xs)+T/2)).astype(np.uint8),C.pad,C.noy
