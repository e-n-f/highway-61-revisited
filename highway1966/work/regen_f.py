import pickle
from family import glyph
SER={'B':(0.654,85),'C':(0.78,95),'D':(0.89,102),'E':(1.0,113),'F':(1.10,123)}
for s,(f,stem) in SER.items():
    d=pickle.load(open('lcE_%s.pkl'%s,'rb'))
    b,ox,oy=glyph('f',f,stem); d['f']=(b,ox,oy,'even-'+s)
    pickle.dump(d,open('lcE_%s.pkl'%s,'wb'))
