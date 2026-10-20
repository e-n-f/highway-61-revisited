"""Specimen of the constructed lower case: python3 lcspecimen.py -> outputs/specimen-lc-E.png
Top: words in Highway 61 Revisited E (constructed letters in black, the rest still from Highway Even 1966 in grey).
Bottom: each constructed letter's outline (red) over the scan of the 1966 sheet."""
import os, numpy as np, cv2
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import DecomposingRecordingPen
import lcE
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','outputs')
DIST=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..')
F=TTFont(os.path.join(DIST,'otf','Highway61Revisited-E.otf')); gs=F.getGlyphSet(); cm=F.getBestCmap(); hm=F['hmtx']
DONE=set(lcE.dims())|set(lcE.polys())
def glyph_poly(name,s):
    pen=DecomposingRecordingPen(gs); gs[name].draw(pen); polys=[]; cur=None
    for op,a in pen.value:
        if op=='moveTo': cur=[np.array(a[0])]
        elif op=='lineTo': cur.append(np.array(a[0]))
        elif op=='curveTo':
            p0=cur[-1]; p1,p2,p3=map(np.array,a)
            for t in np.linspace(0,1,13)[1:]: cur.append((1-t)**3*p0+3*(1-t)**2*t*p1+3*(1-t)*t*t*p2+t**3*p3)
        elif op=='closePath': polys.append(np.array(cur)*s)
    return polys
def line(img,text,x,y,s):
    """PIL text, constructed letters (and capitals) in black, the rest in grey"""
    from PIL import Image, ImageDraw, ImageFont
    font=ImageFont.truetype(os.path.join(DIST,'otf','Highway61Revisited-E.otf'),int(1000*s))
    im=Image.fromarray(img); d=ImageDraw.Draw(im)
    for ch in text:
        col=(0,0,0) if (ch in DONE or ch.isupper() or ch==' ') else (150,150,150)
        d.text((x,y),ch,font=font,fill=col,anchor='ls'); x+=d.textlength(ch,font=font)
    img[:]=np.array(im)
def main():
    W=2400; img=np.full((700,W,3),255,np.uint8)
    for i,t in enumerate(['Highway 61 Revisited', 'jumpy quiz of lowercase', 'Exit 42 Mendocino', 'abcdefghijklmnopqrstuvwxyz']):
        line(img,t,40,170+i*170,0.15 if i<3 else 0.11)
    from lcmodel import overlay
    import lcgeom
    tiles=[]
    ser=lcE.series()
    for ch in 'abcdefghijklmnopqrstuvwxyz':
        if ch not in ser: continue
        g=ser[ch][0]
        if ch in lcE.dims(): g=__import__('shapely').affinity.translate(g,lcE.dims()[ch].get('ox',0)*lcE.SQ)
        fn='/tmp/claude-0/-home-user-standard-alphabets/d2c3cfcf-66d0-5f41-8524-bbf06168e396/scratchpad/_lc_%s.png'%ch; overlay(ch,g,fn,scale=0.3); tiles.append(cv2.imread(fn))
    h=max(t.shape[0] for t in tiles); rows=[]; cur=[]; w=0
    for t in tiles+[None]:
        if t is None or w+t.shape[1]>W:
            r=np.hstack([np.vstack([x,np.full((h-x.shape[0],x.shape[1],3),255,np.uint8)]) for x in cur])
            rows.append(np.hstack([r,np.full((h,W-r.shape[1],3),255,np.uint8)])); cur=[]; w=0
            if t is None: break
        cur.append(t); w+=t.shape[1]
    cv2.imwrite(os.path.join(OUT,'specimen-lc-E.png'),np.vstack([img]+rows))
main()
