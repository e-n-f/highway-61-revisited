import pickle,numpy as np,sys
from capsbuild import *
sets={'B':{4:['ABC','DEF','GHI'],5:['JKL','MNO','PQR'],6:['STU','VWX','YZ1'],7:['234','567','890']},
      'F':{20:['ABC','DEF','GHI'],21:['JKL','MNO','PQR'],22:['STU','VWX','YZ1'],23:['234','567','890']}}
s=sys.argv[1]
E=extract(sets[s])
CP=float(np.median([d['cap'] for d in E.values() if d['cap']])); print(s,'cap px',CP)
out={}
for l,d in E.items():
    b,ox,oy=rasterize(d,CP); out[l]=(b,ox,oy,'series'+s)
pickle.dump(out,open('caps_%s.pkl'%s,'wb'))
