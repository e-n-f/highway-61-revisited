import pickle,numpy as np
from build import smooth
SER={'B':85,'C':95,'D':102,'E':113,'F':123}
for s,T in SER.items():
    d=pickle.load(open('lcE_%s.pkl'%s,'rb'))
    for l in 'bdhkltf':
        b,ox,oy,tag=d[l]; b=b.copy()
        for y in range(530,700):
            r=int(round(oy-y))
            row=b[r]; xs=np.nonzero(row)[0]
            if len(xs)==0: continue
            # runs
            splits=np.where(np.diff(xs)!=1)[0]
            starts=np.r_[xs[0],xs[splits+1]]; ends=np.r_[xs[splits],xs[-1]]
            for a,z in zip(starts,ends):
                w=z-a+1
                if T<w<=T+40:
                    c=(a+z)/2; row[a:z+1]=0; row[int(round(c-T/2)):int(round(c+T/2))]=1
        d[l]=(smooth(b,2),ox,oy,tag)
    pickle.dump(d,open('lcE_%s.pkl'%s,'wb'))
    b,ox,oy,_=d['h']; xs=np.nonzero(b[int(oy-580)])[0]; print(s,'h@580',xs[-1]-xs[0]+1,'target',T)
