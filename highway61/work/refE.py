import pickle,numpy as np
from src import raster
from build import smooth
L=pickle.load(open('L.pkl','rb'))
out={}
for l in 'abcdefghijklmnopqrstuvwxyz':
    d=L[l]
    b,ox,oy,S=raster(d['img'],d['box'],d['base'],115.5,thr=128,margin=8,sigma=0.8)
    b=smooth(b,4)
    out[l]=(b,ox,oy,'ref1966')
pickle.dump(out,open('lcRef.pkl','wb'))
from src import stem_width
for l in 'nlh':
    b,ox,oy,_=out[l]; print(l,stem_width(b,oy,0.4))
ys,xs=np.nonzero(out['h'][0]); print('h top',out['h'][2]-ys.min())
