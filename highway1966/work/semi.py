import pickle,numpy as np
import build
from build import synth2, smooth, L
A=pickle.load(open('A.pkl','rb'))   # narrow sign lowercase (unlightened) for width targets
T=100; R=17; STEM_SPEC=137
build.ERODE_R=R; build.STEM_T=T; build.DOT=103
build.DK=(1+230/186.0)/2
S=525/115.5
def width(b): xs=np.nonzero(b.any(0))[0]; return xs[-1]-xs[0]+1
out={}
for l in 'abcdefghijklmnopqrstuvwxyz':
    Wstd=L[l]['box'][2]*S
    Wn=width(A[l][0])
    Wt=(Wstd+Wn)/2
    fx=0.92 if l in 'ijl' else (Wt-T)/(Wstd-STEM_SPEC)
    fx=min(max(fx,0.6),1.15)
    build.FX[l]=fx
    b,ox,oy=synth2(l); b=smooth(b,4)
    out[l]=(b,ox,oy,'semi')
    print(l,'std %d narrow %d target %d got %d fx %.2f'%(Wstd,Wn,Wt,width(b),fx))
pickle.dump(out,open('A5.pkl','wb'))
