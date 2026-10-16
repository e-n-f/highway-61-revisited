import cv2,numpy as np,pickle
from skimage.morphology import skeletonize
import hires
from hires import page,boxes,LAY
def clusters(d):
    sk=skeletonize(d>0).astype(np.uint8)*255
    segs=cv2.HoughLinesP(sk,1,np.pi/720,50,minLineLength=50,maxLineGap=30)
    segs=[] if segs is None else [s[0] for s in segs]
    cl=[]
    for x1,y1,x2,y2 in segs:
        p=np.array([x1,y1],float); q=np.array([x2,y2],float)
        th=np.arctan2(y2-y1,x2-x1)%np.pi
        for c in cl:
            dth=min(abs(th-c['th']),np.pi-abs(th-c['th']))
            n=np.array([-np.sin(c['th']),np.cos(c['th'])])
            if dth<np.radians(2.5) and abs(n@p-c['rho'])<14 and abs(n@q-c['rho'])<14:
                c['segs'].append((p,q)); break
        else:
            n=np.array([-np.sin(th),np.cos(th)]); cl.append({'th':th,'rho':n@p,'segs':[(p,q)]})
    out=[]
    for c in cl:
        P=np.array([v for s in c['segs'] for v in s]); m=P.mean(0); u,s,vt=np.linalg.svd(P-m); dv=vt[0]
        iv=[]
        for p,q in c['segs']:
            a,b=sorted([(p-m)@dv,(q-m)@dv]); iv.append([a,b])
        iv.sort(); merged=[]
        for a,b in iv:
            if merged and a<=merged[-1][1]+40: merged[-1][1]=max(merged[-1][1],b)
            else: merged.append([a,b])
        for a,b in merged:
            if b-a>50: out.append((m+dv*a,m+dv*b))
    return out
def letter_lines():
    res={}
    for p in ['0031','0032','0033']:
        d=cv2.imread('hr/mask%s.png'%p,0)
        segs=clusters(d)
        o,ang,k=page(p); bx=boxes(o,LAY[p])
        H,W=o.shape; M=cv2.getRotationMatrix2D((W/2,H/2),ang,1)
        def tr(pt): return M@np.array([pt[0],pt[1],1.0])
        segs=[(tr(a),tr(b)) for a,b in segs]
        S=525/457.75; m=30
        for l,(x,y,w,h) in bx.items():
            ex=60
            mine=[]
            for a,b in segs:
                # clip segment to expanded box
                t=np.linspace(0,1,200); P=a[None]+(b-a)[None]*t[:,None]
                inside=(P[:,0]>x-ex)&(P[:,0]<x+w+ex)&(P[:,1]>y-ex)&(P[:,1]<y+h+ex)
                if inside.sum()<3: continue
                Q=P[inside]
                A=((Q[0][0]-(x-m))*S,(Q[0][1]-(y-m))*S); B=((Q[-1][0]-(x-m))*S,(Q[-1][1]-(y-m))*S)
                mine.append((np.array(A),np.array(B)))   # (col,row) in letter bitmap frame
            res[l]=mine
    return res
if __name__=='__main__':
    R=letter_lines(); pickle.dump(R,open('letterlines.pkl','wb'))
    print({l:len(v) for l,v in sorted(R.items())})
    # visual check on master bitmaps
    M=pickle.load(open('hr_lc_sq.pkl','rb')); tiles=[]
    for l in 'abcdefghijklmnopqrstuvwxyz':
        b,ox,oy,_=M[l]; im=cv2.cvtColor(255-b*255,cv2.COLOR_GRAY2BGR)
        for A,B in R.get(l,[]): cv2.line(im,tuple(int(v) for v in A),tuple(int(v) for v in B),(0,0,255),4)
        tiles.append(cv2.resize(im,(260,int(im.shape[0]*260/im.shape[1]))))
    Hm=max(t.shape[0] for t in tiles); rows=[]
    for i in range(0,26,9):
        rows.append(np.hstack([np.pad(t,((0,Hm-t.shape[0]),(0,10),(0,0)),constant_values=255) for t in tiles[i:i+9]]))
    Wm=max(r.shape[1] for r in rows)
    cv2.imwrite('ll_check.png',cv2.resize(np.vstack([np.pad(r,((0,10),(0,Wm-r.shape[1]),(0,0)),constant_values=255) for r in rows]),None,fx=0.6,fy=0.6))
