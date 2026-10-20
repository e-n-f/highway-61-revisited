"""One-line specimen of each font: python3 linespecimen.py -> ../../specimen/specimen-lines-h61.png (the Highway 61 Revisited
fonts) and ../../specimen/specimen-lines-fhwa.png (the FHWA 2000 fonts), a line per font.

The text is the capitals, the lower case, the figures and , . " - ! /, set at equal cap height with each font's own
advances and kerning.  (No apostrophe: the FHWA fonts have none.)  The FHWA fonts' space (added by fhwa2otf.py) is
replaced by the Highway 61 Revisited space of the same series, in cap heights, as in fhwaspecimen.py."""
import os
from PIL import Image, ImageDraw, ImageFont
from fhwaspecimen import info
W=os.path.dirname(os.path.abspath(__file__)); DIST=os.path.join(W,'..','..')
FH=os.path.join(W,'..','fhwa2000')
TEXT='ABCDEFGHIJKLMNOPQRSTUVWXYZ abcdefghijklmnopqrstuvwxyz 0123456789 ,."-!/'
SER=[('A',None),('B','B'),('C','C'),('D','D'),('E','E'),('EM','ModE'),('F','F')]
CAP=100; PAD=30
def line(fn,text,space,label):
    """the text at cap height CAP, words set with the font, spaces of `space` cap heights; label above"""
    cap,_=info(fn); f=ImageFont.truetype(fn,int(round(CAP/cap)))
    lab=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',22)
    d=ImageDraw.Draw(Image.new('L',(1,1))); words=text.split(' ')
    w=sum(d.textlength(t,font=f) for t in words)+space*CAP*(len(words)-1)
    im=Image.new('L',(int(w)+2*PAD,int(CAP*2.1)+40),255); d=ImageDraw.Draw(im)
    d.text((PAD,28),label,font=lab,fill=110,anchor='ls'); x=PAD; y=40+int(CAP*1.35)
    for t in words:
        d.text((x,y),t,font=f,fill=0,anchor='ls'); x+=d.textlength(t,font=f)+space*CAP
    return im
def sheet(ims,fn):
    out=Image.new('L',(max(i.width for i in ims),sum(i.height for i in ims)),255); y=0
    for i in ims: out.paste(i,(0,y)); y+=i.height
    out.save(os.path.join(DIST,'specimen',fn))
def main():
    ours,fhwa=[],[]
    for S,F in SER:
        fn=os.path.join(DIST,'otf','Highway61Revisited-%s.otf'%S); sp=info(fn)[1]
        ours.append(line(fn,TEXT,sp,'Highway 61 Revisited '+('E(M)' if S=='EM' else S)))
        fh=os.path.join(FH,'Series%s2000-Normal.otf'%F) if F else None
        if fh and os.path.exists(fh):
            fhwa.append(line(fh,TEXT,sp,'FHWA 2000 Series '+('E Modified' if F=='ModE' else F)))
    sheet(ours,'specimen-lines-h61.png'); sheet(fhwa,'specimen-lines-fhwa.png')
if __name__=='__main__': main()
