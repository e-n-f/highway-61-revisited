import sys
from PIL import Image,ImageDraw,ImageFont
f=sys.argv[1]; out=sys.argv[2]; size=int(sys.argv[3]); lines=sys.argv[4:]
F=ImageFont.truetype(f,size)
W=max(int(F.getlength(t)) for t in lines)+60; lh=int(size*1.3)
im=Image.new('L',(W,lh*len(lines)+20),255); d=ImageDraw.Draw(im)
for i,t in enumerate(lines): d.text((30,10+lh*i),t,font=F,fill=0)
im.save(out)
