import sys
from PIL import Image,ImageDraw,ImageFont
f=sys.argv[1]; out=sys.argv[2]
lines=sys.argv[3:] or ['abcdefghijklm','nopqrstuvwxyz','bayshore fwy','business district']
F=ImageFont.truetype(f,180)
W=2600;H=230*len(lines)+40
im=Image.new('L',(W,H),255); d=ImageDraw.Draw(im)
for i,t in enumerate(lines): d.text((30,20+230*i),t,font=F,fill=0)
im.save(out)
