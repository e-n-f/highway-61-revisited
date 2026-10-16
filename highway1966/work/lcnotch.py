"""Find notches and ledges in the lower-case outlines: corners (tangent turns over 3 degrees) where one of the two
segments meeting there is short.  python3 lcnotch.py [letters]"""
import sys, numpy as np, lcE, lcoutline
def tangents(s):
    if s[0]=='L':
        d=s[2]-s[1]; d=d/np.linalg.norm(d); return d,d,np.linalg.norm(s[2]-s[1])
    b=lcoutline.beziers(s)
    t0=b[0][1]-b[0][0]; t1=b[-1][3]-b[-1][2]
    L=sum(np.linalg.norm(x[3]-x[0]) for x in b)
    return t0/np.linalg.norm(t0),t1/np.linalg.norm(t1),L
def check(ch,g,prims,short=20.0,ang=3.0):
    out=[]
    for segs in lcoutline.contours(g,prims):
        T=[tangents(s) for s in segs]
        for i in range(len(segs)):
            a,b=T[i-1],T[i]
            turn=np.degrees(np.arccos(np.clip(a[1]@b[0],-1,1)))
            if turn>ang and min(a[2],b[2])<short:
                out.append((tuple(np.round(segs[i][1],1)),round(turn,1),round(a[2],1),round(b[2],1)))
    return out
if __name__=='__main__':
    ser=lcE.series()
    for ch in (sys.argv[1] if len(sys.argv)>1 else sorted(ser)):
        r=check(ch,*ser[ch])
        if r: print(ch,r)
