"""Check the stroke-changed outlines in the built fonts: each glyph's flattened contours must form a valid shape, and
no corner may touch a segment shorter than 6 units that the source outline did not have.
python3 strokecheck.py E|EM"""
import sys, numpy as np
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen
from shapely.geometry import Polygon
from shapely.validation import explain_validity
S=sys.argv[1]
F=TTFont(f'../../otf/Highway61Revisited-{S}.otf'); gs=F.getGlyphSet()
names=[c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ']+['zero','one','two','three','four','five','six','seven','eight','nine']+list('abcdefghijklmnopqrstuvwxyz')
for n in names:
    rp=RecordingPen(); gs[n].draw(rp); rings=[]; cur=[]; segs=[]
    for op,a in rp.value:
        if op=='moveTo': cur=[np.array(a[0])]; start=cur[0]
        elif op=='lineTo': segs.append(('L',cur[-1],np.array(a[0]))); cur.append(np.array(a[0]))
        elif op=='curveTo':
            p0=cur[-1]; p1,p2,p3=map(np.array,a)
            for t in np.linspace(0,1,65)[1:]: cur.append((1-t)**3*p0+3*(1-t)**2*t*p1+3*(1-t)*t*t*p2+t**3*p3)
            segs.append(('C',p0,p3))
        elif op=='closePath':
            R=np.array(cur)
            if len(R)>1 and np.hypot(*(R[-1]-R[0]))<0.01: R=R[:-1]     # CFF rounding of the closing point
            rings.append(R)
    bad=[]
    for r in rings:
        P=Polygon(r)
        if not P.is_valid: bad.append(explain_validity(P))
    short=[(s[1].round(1),s[2].round(1)) for s in segs if np.hypot(*(s[2]-s[1]))<6]
    if bad or short: print(n,'INVALID '+'; '.join(bad) if bad else '', 'short:',short[:4])
