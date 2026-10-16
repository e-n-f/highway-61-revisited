import pickle,numpy as np
import hires
from hires import run
from hrclean import clean
from dents import find_dents
from edgefix import fix
pages=['em035','em036','em037','em038']
rows=[['ABC','DEF','GHI'],['JKL','MNO','PQR'],['STU','VWX','YZ1'],['234','567','890']]
for p,r in zip(pages,rows): hires.LAY[p]=r
out,S=run(pages,'Q','CGOSJ0235689','BDEFHIKLMNPRTUVWXYZ14',690)
d={}
for k,(b,ox,oy,t) in out.items():
    b=clean(b)
    for side in 'LR':
        for a,z in find_dents(b,side):
            try: b=fix(b,side,a,z,pad=4); print('dent',k,side)
            except Exception: pass
    d[k]=(b,ox,oy,t)
pickle.dump(d,open('caps_EM.pkl','wb'))
from src import stem_width
for l in 'HIT':
    b,ox,oy,_=d[l]; print(l,stem_width(b,oy,300/525))
for k in sorted(d):
    b,ox,oy,_=d[k]; ys,xs=np.nonzero(b); print(k,xs.max()-xs.min()+1,round(oy-ys.min()),round(oy-ys.max()),end=' | ')
