import numpy as np,cv2,pickle
from family import graph,prune,smooth_path
from narrowmk import M,ML,constrain_narrow
import narrowmk
def unit(v): return v/np.hypot(*v)
def skeleton_paths(l):
    b,ox,oy,_=M[l]
    # drop the dots of i/j before skeletonising
    n,lab,st,_=cv2.connectedComponentsWithStats(b)
    main=np.argmax(st[1:,4])+1
    dots=[(st[i],np.nonzero(lab==i)) for i in range(1,n) if i!=main and st[i][4]>500]
    bm=(lab==main).astype(np.uint8)
    paths,nodes,dt=graph(bm); paths=prune(paths,nodes,dt)
    out=[]
    for p in paths:
        P=smooth_path(p['pts'],p['closed'])[:, ::-1].astype(float)   # (x,row)
        ends=[]
        for k,(key,idx) in enumerate([('a',0),('z',-1)]):
            ends.append(bool(p[key]) and nodes[p[key]]['deg']==1)
        out.append((P,p['closed'],ends))
    return out,dots,oy
def build(l,f,T,pad=40):
    paths,dots,oym=skeleton_paths(l)
    H=1300; noy=900.0; W=int(600*max(f,1)+2*pad+300)
    img=np.zeros((H,W),np.uint8)
    mp=lambda P: np.column_stack([pad+P[:,0]*f, noy-(oym-P[:,1])])
    for P,closed,ends in paths:
        Q=mp(P)
        if len(Q)<2: continue
        # extend free ends straight on, so the cuts (not the round brush) define them
        if not closed:
            if ends[0] and len(Q)>4: d=unit(Q[0]-Q[min(4,len(Q)-1)]); Q=np.vstack([Q[0]+d*T,Q])
            if ends[1] and len(Q)>4: d=unit(Q[-1]-Q[max(-5,-len(Q))]); Q=np.vstack([Q,Q[-1]+d*T])
        cv2.polylines(img,[np.round(Q*16).astype(np.int32)],closed,1,thickness=max(1,int(round(T))),lineType=cv2.LINE_8,shift=4)
    b,ox,oy,_=M[l]
    for st,(ys,xs) in dots:
        cx=pad+xs.mean()*f; cy=noy-(oym-ys.mean()); r=(st[2]+st[3])/4*(T/140.0)
        cv2.circle(img,(int(round(cx*16)),int(round(cy*16))),int(round(r*16)),1,-1,shift=4)
    # hard cuts at every marked stroke end: terminals, baseline, x-height, ascender/descender ends
    H_,W_=img.shape; yy,xx=np.mgrid[0:H_,0:W_]; P=np.stack([xx,yy],-1).astype(float)
    for L in ML.get(l,[]):
        u=L['u']; ang=np.degrees(np.arctan2(abs(u[1]),abs(u[0])))
        horiz=ang<3; term=(3<=ang<=60) and not (l in 'kvwxyz' and ang>=35)
        if not (horiz or term): continue
        p0=L['c']+u*L['t0']; p1=L['c']+u*L['t1']
        q0=np.array([pad+p0[0]*f, noy-(oym-p0[1])]); q1=np.array([pad+p1[0]*f, noy-(oym-p1[1])])
        uu=unit(q1-q0) if horiz else (u if u@(q1-q0)>=0 else -u)
        nn=np.array([-uu[1],uu[0]]); c=(q0+q1)/2
        half=abs((q1-q0)@uu)/2+(T*0.6 if term else 2)
        side=L['side'] if nn@L['n']>0 else -L['side']
        d=(P-c)@nn; t=(P-c)@uu
        zone=(np.abs(t)<=half)&(d*side<0)&(np.abs(d)<=1.6*T)
        img[zone]=0
    return img,pad,noy
def build_all(s,f,T,letters):
    narrowmk.T_CUR=T
    out={}
    for l in letters:
        img,pad,noy=build(l,f,T)
        nb,_=constrain_narrow(l,img,pad,noy,f)
        out[l]=(nb,pad,noy,'brush-'+s)
    return out
