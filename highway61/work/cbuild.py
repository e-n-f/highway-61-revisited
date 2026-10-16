import pickle,numpy as np
from capsbuild import *
PC={8:['ABC','DEF','GHI'],9:['JKL','MNO','PQR'],10:['STU','VWX','YZ1'],11:['234','567','890']}
E=extract(PC)
CP=float(np.median([d['cap'] for d in E.values() if d['cap']])); print('cap px',CP)
A3={}
for l,d in E.items():
    b,ox,oy=rasterize(d,CP); A3[l]=(b,ox,oy,'seriesC')
pickle.dump(A3,open('A3.pkl','wb'))
