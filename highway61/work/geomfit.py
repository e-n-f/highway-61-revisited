import pickle,numpy as np,cv2
from scipy.optimize import least_squares
M=pickle.load(open('master_mk.pkl','rb')); ML=pickle.load(open('master_lines.pkl','rb'))
def fit_super(P,n_fix=None):
    x,y=P[:,0].astype(float),P[:,1].astype(float)
    p0=[x.mean(),y.mean(),(x.max()-x.min())/2+1,(y.max()-y.min())/2+1,2.0]
    def res(p):
        cx,cy,a,b,n=p
        if n_fix: n=n_fix
        return ((np.abs((x-cx)/a)**n+np.abs((y-cy)/b)**n)**(1/n)-1)*min(a,b)
    s=least_squares(res,p0,bounds=([-1e4,-1e4,1,1,1.5],[1e4,1e4,1e4,1e4,5]))
    return s.x,np.abs(s.fun).max()
def stem_edges(l):
    oy=M[l][2]; out=[]
    for L in ML[l]:
        if abs(L['u'][1])>0.99: out.append(L['c'][0])
    return sorted(set(round(v) for v in out))
def parts(l):
    b,ox,oy,_=M[l]
    cs,h=cv2.findContours(b.copy(),cv2.RETR_CCOMP,cv2.CHAIN_APPROX_NONE)
    outer=[c[:,0,:] for i,c in enumerate(cs) if h[0][i][3]==-1]
    holes=[c[:,0,:] for i,c in enumerate(cs) if h[0][i][3]!=-1]
    return b,oy,max(outer,key=len),(max(holes,key=len) if holes else None)
if __name__=='__main__':
    res={}
    for l in 'obdpqce':
        b,oy,O,Hh=parts(l); se=stem_edges(l)
        if l in 'bp': cut=lambda P: P[P[:,0]>max(se)+12] if se else P
        elif l in 'dq': cut=lambda P: P[P[:,0]<min(se)-12] if se else P
        else: cut=lambda P: P
        Po=cut(O)
        # bowls: keep only the x-height zone
        Po=Po[(oy-Po[:,1]>-20)&(oy-Po[:,1]<560)]
        Ph=cut(Hh) if Hh is not None else None
        po,eo=fit_super(Po,2.0); 
        txt='%s stems %s | outer cx %.0f cy %.0f a %.0f b %.0f err %.1f'%(l,se,po[0],oy-po[1],po[2],po[3],eo)
        if Ph is not None and len(Ph)>20:
            ph,eh=fit_super(Ph,2.0); txt+=' | inner cx %.0f cy %.0f a %.0f b %.0f err %.1f'%(ph[0],oy-ph[1],ph[2],ph[3],eh)
            res[l]=(po,ph,oy,se)
        print(txt)
    pickle.dump(res,open('geomfit.pkl','wb'))
