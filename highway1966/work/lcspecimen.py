"""Specimen of the constructed lower case: python3 lcspecimen.py -> outputs/specimen-lc-E.png
Top: words in Highway 61 Revisited E (constructed letters in black, the rest still from Highway Even 1966 in grey).
Bottom: each constructed letter's outline (red) over the scan of the 1966 sheet."""
import os, numpy as np, cv2
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import DecomposingRecordingPen
import lcE
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','outputs')
F=TTFont(os.path.join(OUT,'Highway61Revisited-E.otf')); gs=F.getGlyphSet(); cm=F.getBestCmap(); hm=F['hmtx']
DONE=set(lcE.dims())
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
    font=ImageFont.truetype(os.path.join(OUT,'Highway61Revisited-E.otf'),int(1000*s))
    im=Image.fromarray(img); d=ImageDraw.Draw(im)
    for ch in text:
        col=(0,0,0) if (ch in DONE or ch.isupper() or ch==' ') else (150,150,150)
        d.text((x,y),ch,font=font,fill=col,anchor='ls'); x+=d.textlength(ch,font=font)
    img[:]=np.array(im)
def main():
    W=2000; img=np.full((700,W,3),255,np.uint8)
    for i,t in enumerate(['one gang', 'an age ago', 'Orange Neon', 'abcdefghijklmnopqrstuvwxyz']):
        line(img,t,40,150+i*160,0.22 if i<3 else 0.11)
    from lcmodel import overlay
    import lcgeom
    tiles=[]
    for ch,d in lcE.dims().items():
        g=lcE.LETTERS[ch]({k:(v if k.endswith('th') else v*lcE.SQ) for k,v in d.items()})
        fn='/tmp/claude-0/-home-user-standard-alphabets/d2c3cfcf-66d0-5f41-8524-bbf06168e396/scratchpad/_lc_%s.png'%ch; overlay(ch,g,fn,scale=0.6); tiles.append(cv2.imread(fn))
    h=max(t.shape[0] for t in tiles)
    row=np.hstack([np.vstack([t,np.full((h-t.shape[0],t.shape[1],3),255,np.uint8)]) for t in tiles])
    if row.shape[1]<W: row=np.hstack([row,np.full((h,W-row.shape[1],3),255,np.uint8)])
    cv2.imwrite(os.path.join(OUT,'specimen-lc-E.png'),np.vstack([img,row[:, :W]]))
main()
