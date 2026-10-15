import sys,pickle,numpy as np
import hires
from hires import run
from hrclean import clean
from dents import find_dents
from edgefix import fix
P={'B':['0004','0005','0006','0007'],'C':['0008','0009','0010','0011'],'D':['0012','0013','0014','0015'],'F':['0020','0021','0022','0023']}
rows=[['ABC','DEF','GHI'],['JKL','MNO','PQR'],['STU','VWX','YZ1'],['234','567','890']]
s=sys.argv[1]
for p,r in zip(P[s],rows): hires.LAY[p]=r
out,S=run(P[s],'Q','CGOSJ0235689','BDEFHIKLMNPRTUVWXYZ14',690)
pickle.dump(out,open('hr_caps_%s.pkl'%s,'wb'))
d={}
for k,(b,ox,oy,t) in out.items():
    b=clean(b)
    for side in 'LR':
        for a,z in find_dents(b,side):
            try: b=fix(b,side,a,z,pad=4); print('dent',k,side)
            except Exception: pass
    d[k]=(b,ox,oy,t)
pickle.dump(d,open('caps_%s.pkl'%s,'wb'))
for k in sorted(d):
    b,ox,oy,_=d[k]; ys,xs=np.nonzero(b); print(k,xs.max()-xs.min()+1,round(oy-ys.min()),round(oy-ys.max()),end=' | ')
