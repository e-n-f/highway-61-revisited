"""Corners next to short pieces in a built font's lower case (notches, ears, steps):
python3 fontnotch.py <font.otf> [letters]  -> prints (point, turn in degrees, lengths of the two pieces)"""
import sys, numpy as np
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen
F=TTFont(sys.argv[1]); gs=F.getGlyphSet()
letters=sys.argv[2] if len(sys.argv)>2 else 'abcdefghijklmnopqrstuvwxyz'
for ch in letters:
    rp=RecordingPen(); gs[ch].draw(rp); cons=[]; cur=None
    for op,a in rp.value:
        if op=='moveTo': cur=[]; start=np.array(a[0]); last=start
        elif op=='lineTo': p=np.array(a[0]); cur.append((last,last,p,p,np.hypot(*(p-last)))); last=p
        elif op=='curveTo':
            c1,c2,p=map(np.array,a); L=np.hypot(*(c1-last))+np.hypot(*(c2-c1))+np.hypot(*(p-c2))
            cur.append((last,c1,c2,p,L)); last=p
        elif op=='closePath':
            if np.hypot(*(last-start))>1e-6: cur.append((last,last,start,start,np.hypot(*(start-last))))
            cons.append(cur)
    out=[]
    for c in cons:
        for i in range(len(c)):
            a,b=c[i-1],c[i]
            ta=a[3]-(a[2] if np.hypot(*(a[3]-a[2]))>1e-6 else a[0]); tb=(b[1] if np.hypot(*(b[1]-b[0]))>1e-6 else b[3])-b[0]
            if np.hypot(*ta)<1e-9 or np.hypot(*tb)<1e-9: continue
            turn=np.degrees(np.arccos(np.clip(ta@tb/np.hypot(*ta)/np.hypot(*tb),-1,1)))
            if turn>3 and min(a[4],b[4])<6: out.append((tuple(np.round(b[0],1)),round(turn,1),round(a[4],1),round(b[4],1)))
    if out: print(ch,out)
