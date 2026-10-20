"""Comparison specimen, Highway 61 Revisited against the FHWA 2000 Standard Alphabets fonts:
python3 fhwacompare.py -> ../outputs/compare-fhwa2000.png

For each series: a line in Highway 61 Revisited (black), the same line in the FHWA 2000 font (blue), both at the same
cap height (the FHWA fonts' caps are 285 units of their 1000-unit em, ours 690), and an overlay of every letter and
figure, the pair centred on each other on a common baseline: Highway 61 Revisited only in red, FHWA 2000 only in blue,
both in black.  Words are set one by one with the same space (0.3 cap height) in both, since the FHWA subsets carry
no space glyph (nor apostrophe or curly quotes, hence the straight double quote)."""
import os, numpy as np
from PIL import Image, ImageDraw, ImageFont
from fontTools.ttLib import TTFont
from fontTools.pens.boundsPen import BoundsPen
W=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(W,'..','outputs'); DIST=os.path.join(W,'..','..'); FH=os.path.join(W,'..','fhwa2000')
PAIRS=[('B','B'),('C','C'),('D','D'),('E','E'),('EM','ModE'),('F','F')]
CAP=120                                         # cap height, pixels
TEXT=['Exit 42 Mendocino Highway','Is "Hi-Fi!" said Jo. 1-2-3']
CHARS='ABCDEFGHIJKLMNOPQRSTUVWXYZ abcdefghijklmnopqrstuvwxyz 0123456789'
def cap_units(fn):
    f=TTFont(fn); g=f.getGlyphSet(); p=BoundsPen(g); g['H'].draw(p); return p.bounds[3]/f['head'].unitsPerEm
def font(fn): return ImageFont.truetype(fn,int(round(CAP/cap_units(fn))))
def words(d,x,y,text,f,col):
    sp=0.3*CAP
    for w in text.split(' '):
        d.text((x,y),w,font=f,fill=col,anchor='ls'); x+=d.textlength(w,font=f)+sp
def glyph_mask(ch,f,size):
    im=Image.new('L',size,0); ImageDraw.Draw(im).text((size[0]//2,int(size[1]*0.72)),ch,font=f,fill=255,anchor='ms')
    a=np.array(im)>127
    xs=np.where(a.any(0))[0]
    if len(xs): a=np.roll(a,size[0]//2-(xs[0]+xs[-1])//2,axis=1)
    return a
def main():
    rows=[]
    lab=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',30)
    for ours,theirs in PAIRS:
        fo=font(os.path.join(DIST,'otf',f'Highway61Revisited-{ours}.otf'))
        ft=font(os.path.join(FH,f'Series{theirs}2000-Normal.otf'))
        cell=(int(CAP*1.25),int(CAP*1.9)); per=26
        groups=[g for g in CHARS.split(' ')]
        nrow=sum((len(g)+per-1)//per for g in groups)
        H=40+int(CAP*1.45)+len(TEXT)*2*int(CAP*1.55)-int(CAP*0.9)+nrow*cell[1]+20
        im=Image.new('RGB',(per*cell[0]+40,H),'white'); d=ImageDraw.Draw(im)
        name='E(M)' if ours=='EM' else ours
        d.text((20,40),f'Highway 61 Revisited {name}  (black, red)   vs   FHWA 2000 Series {"E Modified" if theirs=="ModE" else theirs}  (blue)',
               font=lab,fill=(110,110,110),anchor='ls')
        y=40+int(CAP*1.45)
        for t in TEXT:
            words(d,20,y,t,fo,(0,0,0)); y+=int(CAP*1.55)
            words(d,20,y,t,ft,(30,80,200)); y+=int(CAP*1.55)
        y-=int(CAP*0.9)
        a=np.array(im)
        for g in groups:
            for k in range(0,len(g),per):
                for i,ch in enumerate(g[k:k+per]):
                    mo=glyph_mask(ch,fo,cell); mt=glyph_mask(ch,ft,cell)
                    x0=20+i*cell[0]; sub=a[y:y+cell[1],x0:x0+cell[0]]
                    sub[mo&~mt]=(225,45,45); sub[mt&~mo]=(40,100,225); sub[mo&mt]=(0,0,0)
                y+=cell[1]
        rows.append(Image.fromarray(a))
    Wd=max(r.width for r in rows); out=Image.new('RGB',(Wd,sum(r.height for r in rows)),'white'); y=0
    for r in rows: out.paste(r,(0,y)); y+=r.height
    out.save(os.path.join(OUT,'compare-fhwa2000.png'))
if __name__=='__main__': main()
