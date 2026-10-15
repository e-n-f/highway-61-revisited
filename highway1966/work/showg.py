import sys
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import DecomposingRecordingPen
from PIL import Image,ImageDraw
from fontTools.pens.basePen import BasePen
f=TTFont(sys.argv[1]); gs=f.getGlyphSet()
names=sys.argv[3:]
sc=0.6; W=int(sum(gs[n].width for n in names)*sc)+40; H=int(1100*sc)
im=Image.new('RGB',(W,H),'white'); d=ImageDraw.Draw(im)
from fontTools.pens.pointPen import SegmentToPointPen
class P(BasePen):
    def __init__(s,gs,dx): super().__init__(gs); s.dx=dx; s.pts=[]; s.cur=None; s.poly=[]; s.polys=[]
    def T(s,p): return (20+(p[0]+s.dx)*sc, H-(p[1]+330)*sc)
    def _moveTo(s,p): s.poly=[s.T(p)]; s.pts.append((s.T(p),1))
    def _lineTo(s,p): s.poly.append(s.T(p)); s.pts.append((s.T(p),1))
    def _curveToOne(s,a,b,c):
        import numpy as np
        p0=np.array(s._getCurrentPoint() if False else s.last)
        for t in np.linspace(0,1,12)[1:]:
            q=(1-t)**3*p0+3*(1-t)**2*t*np.array(a)+3*(1-t)*t*t*np.array(b)+t**3*np.array(c)
            s.poly.append(s.T(q))
        s.pts.append((s.T(c),1))
    def _qCurveToOne(s,a,b):
        import numpy as np
        p0=np.array(s.last)
        for t in np.linspace(0,1,10)[1:]:
            q=(1-t)**2*p0+2*(1-t)*t*np.array(a)+t*t*np.array(b); s.poly.append(s.T(q))
        s.pts.append((s.T(a),0)); s.pts.append((s.T(b),1))
    def _closePath(s): s.polys.append(s.poly)
    def moveTo(s,p): s.last=p; super().moveTo(p)
    def lineTo(s,p): super().lineTo(p); s.last=p
    def qCurveTo(s,*pts):
        super().qCurveTo(*pts); s.last=pts[-1]
x=0
for n in names:
    p=P(gs,x); gs[n].draw(p)
    for poly in p.polys: d.line(poly+[poly[0]],fill='black',width=1)
    for (pt,on) in p.pts:
        r=2; d.ellipse([pt[0]-r,pt[1]-r,pt[0]+r,pt[1]+r],fill='red' if on else 'blue')
    for yy in (0,525,680,-230): d.line([(0,H-(yy+330)*sc),(W,H-(yy+330)*sc)],fill=(200,200,255))
    x+=gs[n].width
im.save(sys.argv[2])
