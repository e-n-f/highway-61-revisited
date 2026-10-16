import pickle,numpy as np
from edgefix import edges,fix
def find_dents(b,side,jump=4,maxlen=30):
    e=edges(b,side); H=len(e); out=[]
    r=1
    while r<H-1:
        if not np.isnan(e[r]) and not np.isnan(e[r-1]) and abs(e[r]-e[r-1])>=jump:
            a=r-1
            for z in range(r+1,min(H,r+maxlen)):
                if np.isnan(e[z]): break
                # back on the line from e[a-6]..e[z+6]?
                if z+6<H and a-6>=0 and not np.isnan(e[a-6]) and not np.isnan(e[z+6]):
                    pred=e[a]+(e[z]-e[a])*0  # placeholder
                if abs(e[z]-e[z-1])>=jump and np.sign(e[z]-e[z-1])!=np.sign(e[r]-e[r-1]):
                    # check endpoints consistent (smooth continuation)
                    if abs((e[z]-e[a]))<=2+0.5*(z-a):
                        out.append((a,z)); r=z; break
        r+=1
    return out
if __name__=='__main__':
    import sys
    for f in ['hr_lc_clean.pkl','caps_E.pkl']:
        d=pickle.load(open(f,'rb')); changed=False
        for k in sorted(d):
            b,ox,oy,t=d[k]
            for side in 'LR':
                ds=find_dents(b,side)
                if ds:
                    print(f[:6],k,side,[(int(oy-z),int(oy-a)) for a,z in ds])
                    for a,z in ds: b=fix(b,side,a,z,pad=4)
            d[k]=(b,ox,oy,t)
        if len(sys.argv)>1: pickle.dump(d,open(f,'wb'))
