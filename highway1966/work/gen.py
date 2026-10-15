import sys,pickle,time
from family import glyph,stem_for,WIDTHS
name=sys.argv[1]; f=WIDTHS[name]
t=time.time(); out={}
for l in 'abcdefghijklmnopqrstuvwxyz':
    b,ox,oy=glyph(l,f,stem_for(f)); out[l]=(b,ox,oy,'fam-'+name)
pickle.dump(out,open('lc_%s.pkl'%name,'wb')); print(name,'%.0fs'%(time.time()-t))
