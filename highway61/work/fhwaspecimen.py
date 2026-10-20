"""Two-column specimen, Highway 61 Revisited against the FHWA 2000 Standard Highway Signs fonts:
python3 fhwaspecimen.py -> ../../specimen/comparison-1966-vs-2000.png

Series B C D on the left, E E(M) F on the right (E(M) against Series E Modified 2000).  Each series: sign copy, the
complete lower case, capitals, and figures with punctuation, five lines in each font, set with the fonts' own
advances and kerning (no tracking) at equal cap height.  The FHWA subsets have no space glyph, so they get the
Highway 61 Revisited space of the same series (in cap heights); nor do they have curly quotes, so their line ends
with a straight double quote."""
import os
from PIL import Image, ImageDraw, ImageFont
from fontTools.ttLib import TTFont
from fontTools.pens.boundsPen import BoundsPen
W=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(W,'..','outputs'); DIST=os.path.join(W,'..','..'); FH=os.path.join(W,'..','fhwa2000')
CAP=64; LH=int(CAP*1.7); LAB=26
COLS=[[('B','B'),('C','C'),('D','D')],[('E','E'),('EM','ModE'),('F','F')]]
LINES=['Bayshore Fwy  Business District  280','Page Mill Rd  Oregon Expwy  Exit 2',
       'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ','0123456789  .,-!/{q}']
def info(fn):
    f=TTFont(fn); g=f.getGlyphSet(); p=BoundsPen(g); g['H'].draw(p); upm=f['head'].unitsPerEm
    sp=f['hmtx']['space'][0] if 'space' in f.getGlyphOrder() else None
    return p.bounds[3]/upm, (sp/p.bounds[3] if sp else None)
class Setter:
    """sets a line at cap height CAP, the words with the font's advances and the given space (cap heights)"""
    def __init__(s,fn,space=None):
        cap,sp=info(fn); s.f=ImageFont.truetype(fn,int(round(CAP/cap))); s.space=(space if space is not None else sp)*CAP
    def width(s,d,t): return sum(d.textlength(w,font=s.f) for w in t.split(' '))+s.space*t.count(' ')
    def draw(s,d,x,y,t):
        for w in t.split(' '):
            d.text((x,y),w,font=s.f,fill=(0,0,0),anchor='ls'); x+=d.textlength(w,font=s.f)+s.space
def main():
    lab=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',LAB)
    labB=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',int(LAB*1.15))
    tit=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',int(LAB*1.15))
    sets={}
    for col in COLS:
        for o,t in col:
            ours=os.path.join(DIST,'otf',f'Highway61Revisited-{o}.otf')
            so=Setter(ours); st=Setter(os.path.join(FH,f'Series{t}2000-Normal.otf'),info(ours)[1])
            sets[o]=(so,st)
    tmp=ImageDraw.Draw(Image.new('RGB',(10,10)))
    gut=int(max(tmp.textlength(w,font=lab) for w in ('Highway 61 Rev.','2000 FHWA')))+30   # room for the row labels
    colw=[gut+max(max(so.width(tmp,l.replace('{q}','”')),st.width(tmp,l.replace('{q}','"'))) for o,_ in col for so,st in [sets[o]] for l in LINES)+40
          for col in COLS]
    secH=int(LAB*2.2)+2*len(LINES)*LH+int(LAB*0.8)
    H=int(LAB*2.4)+3*secH
    colw=[int(c) for c in colw]
    im=Image.new('RGB',(sum(colw)+60,H),'white'); d=ImageDraw.Draw(im)
    d.text((20,int(LAB*1.4)),'Highway 61 Revisited (1961/1966-based) vs. 2000 FHWA Standard Highway Signs fonts '
           '(capitals scaled to equal height, natural spacing)',font=tit,fill=(40,40,40),anchor='ls')
    x0=20
    for ci,col in enumerate(COLS):
        y=int(LAB*2.4)
        for o,t in col:
            so,st=sets[o]
            d.text((x0,y+int(LAB*1.6)),'Series '+('E(M)' if o=='EM' else o),font=labB,fill=(0,0,0),anchor='ls')
            yy=y+int(LAB*2.2)+int(CAP*1.15)
            for who,s,q in (('Highway 61 Rev.',so,'”'),('2000 FHWA',st,'"')):
                for l in LINES:
                    d.text((x0,yy),who,font=lab,fill=(140,140,140),anchor='ls')
                    s.draw(d,x0+gut,yy,l.replace('{q}',q)); yy+=LH
            y+=secH
            d.line((x0,y-int(LAB*0.4),x0+colw[ci]-40,y-int(LAB*0.4)),fill=(190,190,190),width=2)
        if ci==0: d.line((x0+colw[0]-10,int(LAB*2.4),x0+colw[0]-10,H-10),fill=(190,190,190),width=2)
        x0+=colw[ci]+20
    im.save(os.path.join(DIST,'specimen','comparison-1966-vs-2000.png'))
    print(im.size)
if __name__=='__main__': main()
