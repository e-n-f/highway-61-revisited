import pickle,numpy as np,cv2
def find(b,r=6):
    se=cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2*r+1,2*r+1))
    oc=cv2.morphologyEx(cv2.morphologyEx(b,cv2.MORPH_OPEN,se),cv2.MORPH_CLOSE,se)
    diff=(b^oc).astype(np.uint8)
    n,lab,st,cen=cv2.connectedComponentsWithStats(diff)
    return [(int(cen[i][0]),int(cen[i][1]),int(st[i][4])) for i in range(1,n) if 40<st[i][4]<2000],oc,lab,st
if __name__=='__main__':
    for f in ['hr_lc_clean.pkl','caps_E.pkl']:
        d=pickle.load(open(f,'rb'))
        for k in sorted(d):
            b,ox,oy,_=d[k]; fs,_,_,_=find(b)
            if fs: print(k,[(x,int(oy-y),a) for x,y,a in fs])
