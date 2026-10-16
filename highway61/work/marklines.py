import cv2,numpy as np,pickle
from skimage.morphology import skeletonize
def lines_from_mask(d):
    sk=skeletonize(d>0).astype(np.uint8)*255
    segs=cv2.HoughLinesP(sk,1,np.pi/720,50,minLineLength=60,maxLineGap=30)
    segs=[] if segs is None else [s[0] for s in segs]
    L=[]
    for x1,y1,x2,y2 in segs:
        th=np.arctan2(y2-y1,x2-x1)%np.pi; n=np.array([-np.sin(th),np.cos(th)]); rho=n@np.array([x1,y1])
        L.append([th,rho,(x1,y1),(x2,y2)])
    # cluster
    cl=[]
    for th,rho,p,q in L:
        for c in cl:
            dth=min(abs(th-c['th']),np.pi-abs(th-c['th']))
            if dth<np.radians(2.5):
                # compare rho in c's frame
                n=np.array([-np.sin(c['th']),np.cos(c['th'])])
                if abs(n@np.array(p,float)-c['rho'])<14 and abs(n@np.array(q,float)-c['rho'])<14:
                    c['pts']+= [p,q]; break
        else: cl.append({'th':th,'rho':rho,'pts':[p,q]})
    out=[]
    for c in cl:
        P=np.array(c['pts'],float); m=P.mean(0); u,s,vt=np.linalg.svd(P-m); dvec=vt[0]
        t=(P-m)@dvec; a=m+dvec*t.min(); b=m+dvec*t.max()
        if np.hypot(*(b-a))>60: out.append((a,b))
    return out
if __name__=='__main__':
    allL={}
    for p in ['0031','0032','0033']:
        d=cv2.imread('hr/mask%s.png'%p,0)
        Ls=lines_from_mask(d); allL[p]=Ls; print(p,len(Ls))
        im=cv2.cvtColor(cv2.imread('hr/%s.png'%p,0),cv2.COLOR_GRAY2BGR)
        for a,b in Ls: cv2.line(im,tuple(int(v) for v in a),tuple(int(v) for v in b),(0,0,255),5)
        cv2.imwrite('hr/lines%s_s.png'%p,cv2.resize(im,None,fx=0.3,fy=0.3,interpolation=cv2.INTER_AREA))
    pickle.dump(allL,open('marklines.pkl','wb'))
