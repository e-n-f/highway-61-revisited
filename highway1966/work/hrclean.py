import pickle,numpy as np,cv2
from build import smooth
def clean(b):
    b=b.copy()
    # remove specks
    n,lab,st,_=cv2.connectedComponentsWithStats(b)
    b=np.isin(lab,[i for i in range(1,n) if st[i][4]>300]).astype(np.uint8)
    # fill small holes
    inv=(1-b).astype(np.uint8); n,lab,st,_=cv2.connectedComponentsWithStats(inv)
    for i in range(1,n):
        if st[i][4]<900: b[lab==i]=1
    return smooth(b,1.5)
if __name__=='__main__':
    for src,dst in [('hr_caps_E.pkl','caps_E.pkl'),('hr_lc.pkl','hr_lc_clean.pkl')]:
        d=pickle.load(open(src,'rb'))
        pickle.dump({k:(clean(v[0]),)+tuple(v[1:]) for k,v in d.items()},open(dst,'wb'))
    from lighten import lighten
    lc=pickle.load(open('hr_lc_clean.pkl','rb'))
    L=lighten(lc,14); pickle.dump(L,open('lcE_E.pkl','wb'))
    from src import stem_width,vruns
    for l in 'nlo': b,ox,oy,_=L[l]; print(l,stem_width(b,oy,0.4),vruns(b))
