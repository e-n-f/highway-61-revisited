"""Straight-line lower-case letters (k v w x y z): each edge of the scan outline is fitted as a line (the scan's
vertices found by polygon simplification, the edge's middle 80% fitted, the 1.5-unit printing spread taken off along
the normal) and the vertices are the neighbouring lines' intersections, snapped to 1/32 square.
python3 lcpoly.py <letters>  ->  lc1966/poly.json"""
import sys, json, os, numpy as np, cv2
from shapely.geometry import Polygon
from lcview import units, SQ
from lcmodel import SPREAD
FN='lc1966/poly.json'
MINLEN={'z':20.0}
def _line(seg):
    m=len(seg); seg=seg[int(m*0.1):max(int(m*0.9),int(m*0.1)+2)]
    ctr=seg.mean(0); _,_,vt=np.linalg.svd(seg-ctr); return ctr,vt[0]
def _inter(l1,l2):
    (p1,u1),(p2,u2)=l1,l2
    t=np.linalg.solve(np.array([u1,-u2]).T,p2-p1); return p1+u1*t[0]
def fit_letter(ch,eps=6.0,minlen=14.0,minang=6.0,which=0):
    cs=sorted(units(ch),key=len,reverse=True); c=cs[which]
    if not Polygon(c).exterior.is_ccw: c=c[::-1]
    ap=cv2.approxPolyDP(c.astype(np.float32).reshape(-1,1,2),eps,True).reshape(-1,2)
    idx=sorted(int(np.argmin(np.hypot(*(c-v).T))) for v in ap)
    def seg(i,j): return c[i:j+1] if j>i else np.vstack([c[i:],c[:j+1]])
    while True:
        n=len(idx); L=[_line(seg(idx[k],idx[(k+1)%n])) for k in range(n)]
        for k in range(n):
            if L[k][1]@(c[idx[(k+1)%n]]-c[idx[k]])<0: L[k]=(L[k][0],-L[k][1])
        ang=lambda u,v: np.degrees(np.arccos(np.clip(abs(u@v),0,1)))
        V=[_inter(L[k-1],L[k]) for k in range(n)]
        # merge nearly collinear neighbours, drop short edges
        k=next((k for k in range(n) if ang(L[k-1][1],L[k][1])<minang),None)
        if k is not None: idx.pop(k); continue
        lens=[np.hypot(*(V[(k+1)%n]-V[k])) for k in range(n)]
        k=int(np.argmin(lens))
        if lens[k]<minlen and n>3:
            idx.pop((k+1)%n); continue
        break
    out=[]
    for k in range(n):
        p,u=L[k]; nrm=np.array([u[1],-u[0]]); out.append((p-nrm*SPREAD,u))
    return np.array([_inter(out[k-1],out[k]) for k in range(n)])
if __name__=='__main__':
    out=json.load(open(FN)) if os.path.exists(FN) else {}
    for ch in sys.argv[1]:
        V=fit_letter(ch,minlen=MINLEN.get(ch,14.0))/SQ
        V=np.round(V*32)/32
        for lv in (0,6,8,-2.1875): V[np.abs(V[:,1]-lv)<0.05,1]=lv   # on the baseline, x-height, ascender line
        out[ch]=V.tolist()
        print(ch,len(V),[tuple(np.round(v/SQ,3)) for v in V])
    json.dump(out,open(FN,'w'),indent=1)
