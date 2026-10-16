"""Specimen of all six series: python3 lcseries.py -> ../../specimen-lc-series.png
(the lower case, a line of sign copy and a line with punctuation, in each Highway 61 Revisited font)"""
import os
from PIL import Image, ImageDraw, ImageFont
W=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(W,'..','outputs'); DIST=os.path.join(W,'..','..')
def main():
    H=520
    im=Image.new('RGB',(2700,6*H+40),'white'); d=ImageDraw.Draw(im)
    lab=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',36)
    for i,S in enumerate(['B','C','D','E','EM','F']):
        f=ImageFont.truetype(os.path.join(DIST,'Highway61Revisited-%s.otf'%S),140); y=i*H
        d.text((20,y+50),'Highway 61 Revisited '+('E(M)' if S=='EM' else S),font=lab,fill=(120,120,120),anchor='ls')
        d.text((20,y+200),'abcdefghijklmnopqrstuvwxyz',font=f,fill='black',anchor='ls')
        d.text((20,y+340),'Exit 42 Mendocino Highway',font=f,fill='black',anchor='ls')
        d.text((20,y+480),'It’s “Hi-Fi!” said Jo. 1-2-3',font=f,fill='black',anchor='ls')
    im.save(os.path.join(DIST,'specimen-lc-series.png'))
if __name__=='__main__': main()
