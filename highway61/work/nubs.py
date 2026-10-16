import pickle,numpy as np,cv2
def denub(b,r=5,amin=25,amax=400,asp=2.2,report=False):
    se=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*r+1,2*r+1))
    b=b.copy(); found=[]
    for mode in ('bump','dent'):
        if mode=='bump': t=cv2.morphologyEx(b,cv2.MORPH_OPEN,se); diff=(b&(1-t)).astype(np.uint8)
        else: t=cv2.morphologyEx(b,cv2.MORPH_CLOSE,se); diff=(t&(1-b)).astype(np.uint8)
        n,lab,st,cen=cv2.connectedComponentsWithStats(diff)
        for i in range(1,n):
            a=st[i][4]
            if not (amin<a<amax): continue
            pts=np.column_stack(np.nonzero(lab==i)).astype(np.float32)
            (cx,cy),(w,h),ang=cv2.minAreaRect(pts)
            if max(w,h)/max(1,min(w,h))>asp: continue
            m=lab==i
            b[m]= 0 if mode=='bump' else 1
            found.append((mode,int(cen[i][0]),int(cen[i][1]),int(a)))
    return b,found
if __name__=='__main__':
    for f in ['hr_lc_clean.pkl','caps_E.pkl']:
        d=pickle.load(open(f,'rb'))
        for k in sorted(d):
            b,ox,oy,t=d[k]; nb,fd=denub(b)
            if fd: print(k,[(m,x,int(oy-y),a) for m,x,y,a in fd])
