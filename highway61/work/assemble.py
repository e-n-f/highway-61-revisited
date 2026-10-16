import pickle,numpy as np,sys
import geomc2,geomc3,geomc4,geomc5,geomc6,geomc7,geomc8,geomc9,narrowmk
from lighten import lighten
SRC={'o':'g2','e':'g7','j':'g4','f':'g9',**{l:'g5' for l in 'bdpq'},'a':'g9','g':'g8',**{l:'g6' for l in 'cnmhru'},**{l:'g3' for l in 'ilkvwxy'},'t':'g9'}
def make(l,f,T):
    k=SRC[l]
    if k=='g2': return geomc2.build(l,f,T)
    if k=='g4': return geomc4.build(l,f,T)
    if k=='g5': return geomc5.build(l,f,T)
    if k=='g6': return geomc6.build(l,f,T)
    if k=='g7': return geomc7.build(l,f,T)
    if k=='g8': return geomc8.build(l,f,T)
    if k=='g9': return geomc9.build(l,f,T)
    if k=='g3': return geomc3.build(l,f,T)
    return geomc3.build_hooked(l,f,T)
def series(s,f,T,base):
    d=pickle.load(open(base,'rb'))
    for l in SRC:
        ff=0.66 if (s=='B' and l=='x') else f
        img,pad,noy=make(l,ff,T); d[l]=(img,pad,noy,'geom-'+s)
    d['z']=(narrowmk.rebuild_z(np.zeros((1300,int(650*max(f,1)+280)),np.uint8),[],900.0,40,f,T),40,900.0,'geom-'+s)
    if s!='EMC':
        d.update(lighten({'w':d['w']},4))
        if s=='B': d.update(lighten({'x':d['x']},4))
    return d
if __name__=='__main__':
    s=sys.argv[1]
    if s=='EMC':
        M=pickle.load(open('master_mk.pkl','rb'))
        d=series('EMC',1.0,140,'lcE_EMC.pkl'); d['s']=M['s']
    else:
        f,T=(1.0,123) if s=='E' else narrowmk.SER[s]; d=series(s,f,T,'lcE_%s.pkl'%s)
    pickle.dump(d,open('lcE_%s.pkl'%s,'wb')); print(s,'ok')
