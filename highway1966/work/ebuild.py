import pickle,numpy as np
from capsbuild import *
PE={16:['ABC','DEF','GHI'],17:['JKL','MNO','PQR'],18:['STU','VWX','YZ1'],19:['234','567','890']}
E=extract(PE)
CP=float(np.median([d['cap'] for d in E.values() if d['cap']])); print('cap px',CP)
A7={}
for l,d in E.items():
    b,ox,oy=rasterize(d,CP); A7[l]=(b,ox,oy,'seriesE')
pickle.dump(A7,open('A7.pkl','wb'))
