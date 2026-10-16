"""Choose the lower-case width factor of a series so that its alphabet (sum of the a-z advances) is as wide as the
FHWA 2000 font's, both measured in cap heights:  python3 lcwidth.py <series> [f0 f1]

The advances are the build's: the condensed outline's width plus the 1966 base font's side bearings.  The total is
close to linear in f (only the core, width less stroke, scales), so two evaluations and a secant step, then one
more to confirm, settle it."""
import sys, os, string, numpy as np
from fontTools.ttLib import TTFont
from fontTools.pens.boundsPen import BoundsPen
import lccondense, strokeop
W=os.path.dirname(os.path.abspath(__file__))
def target(S):
    f=TTFont(os.path.join(W,'..','fhwa2000','Series%s2000-Normal.otf'%S)); g=f.getGlyphSet(); cm=f.getBestCmap()
    p=BoundsPen(g); g['H'].draw(p)
    return sum(g[cm[ord(c)]].width for c in string.ascii_lowercase)/p.bounds[3]
def total(S,f):
    base=TTFont(os.path.join(W,'..','outputs','HighwayEven1966-%s.ttf'%S)); gs=base.getGlyphSet(); hm=base['hmtx']
    out=lccondense.series(S,f); tot=0
    for c in string.ascii_lowercase:
        P=np.vstack(strokeop.to_polygons(out[c])); bp=BoundsPen(gs); gs[c].draw(bp)
        tot+=round(bp.bounds[0]+(P[:,0].max()-P[:,0].min())+hm[c][0]-bp.bounds[2])
    return tot/690
if __name__=='__main__':
    S=sys.argv[1]; T=target(S)
    f0,f1=(float(sys.argv[2]),float(sys.argv[3])) if len(sys.argv)>3 else (lccondense.F[S],lccondense.F[S]*1.05)
    t0=total(S,f0); print(S,'target',round(T,3),'f',f0,round(t0,3),flush=True)
    t1=total(S,f1); print(S,'f',f1,round(t1,3),flush=True)
    f2=f1+(T-t1)*(f1-f0)/(t1-t0); f2=round(f2*1000)/1000
    t2=total(S,f2); print(S,'f',f2,round(t2,3),'RESULT',flush=True)
