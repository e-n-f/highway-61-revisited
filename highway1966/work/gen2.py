import sys,pickle,time
from family import glyph
SER={'B':(0.654,85),'C':(0.78,95),'D':(0.89,102),'E':(1.0,113),'F':(1.10,123)}
s=sys.argv[1]; f,stem=SER[s]
t=time.time(); out={}
for l in 'abcdefghijklmnopqrstuvwxyz':
    b,ox,oy=glyph(l,f,stem); out[l]=(b,ox,oy,'even-'+s)
pickle.dump(out,open('lcE_%s.pkl'%s,'wb')); print(s,'%.0fs'%(time.time()-t))
