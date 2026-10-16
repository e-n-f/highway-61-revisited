"""Fit the lcE constructions to the scans: python3 lcfitall.py <letters> [--snap]"""
import sys, json, os, numpy as np
import lcE
from lcmodel import fit, overlay, SQ
SCR='/tmp/claude-0/-home-user-standard-alphabets/d2c3cfcf-66d0-5f41-8524-bbf06168e396/scratchpad/'
FN='lc1966/fit.json'
def load():
    return json.load(open(FN)) if os.path.exists(FN) else {}
INIT={
 'o':dict(OL=0.387,Oyl=3.0,OR=5.586,Oyr=3.0,Oxt=2.99,OT=6.08,Oxb=2.99,OB=-0.06,
          IL=1.99,Iyl=3.0,IR=4.02,Iyr=3.0,Ixt=3.0,IT=4.55,Ixb=3.0,IB=1.46),
 'n':dict(w=1.58,x2=3.5,oL=1.0,oyl=4.0,oyr=3.92,oxt=2.79,T=6.09,ccx=2.56,ccy=3.4,ca=1.1,cb=1.05),

 'a':dict(xs=3.39,xr=5.0,aL=-0.0,ayl=4.71,ayr=3.84,axt=2.689,aT=6.087,ux=2.414,uy=2.059,ur=2.595,
          BL=-0.1,Byl=1.99,BR=4.36,Byr=1.99,BR2=4.36,Byr2=1.99,Bxt=2.13,BT=4.05,Bxb=2.13,BB=-0.07,cx=2.44,cy=2.05,ca=1.0,cb=0.707,tx=0.678,ty=4.776,tth=-63.1+180),
 'e':dict(OL=-0.06,Oyl=3.0,OR=5.22,Oyr=2.58,Oxt=2.75,OT=6.12,Oxb=2.53,OB=-0.07,xc=5.05,yb=2.53,yt=3.82,
          ux=2.532,uy=3.68,ua=0.836,ub=0.99,lx=2.66,ly=2.42,la=1.0,lb=1.0,tx=4.14,ty=1.335,tth=142.8),
 'g':dict(xs=3.35,xr=5.0,toa=3.44,toy=-0.075,tob=3.0,tix=1.54,tiy=1.325,tia=2.07,tib=2.97,
          BL=-0.17,Byl=2.98,BR=4.46,Byr=2.98,BR2=4.46,Byr2=2.98,Bxt=2.145,BT=6.08,Bxb=2.145,BB=-0.11,cx=2.494,cy=2.992,ca=1.068,cb=1.514,tx=2.097,ty=-2.36,tth=-66.3+180),
}
TIES={'o':[('Oxb','Oxt'),('Ixb','Ixt'),('Ixt','Oxt')]}
def sc(n): return 1.0 if n.endswith('th') else SQ
def resolve(p,ties):
    for _ in range(5):
        for a,b in ties.items():
            if b in p: p[a]=p[b]
def run(ch,snap=False,fixed=None):
    init=load().get(ch,INIT[ch]); fixed=fixed or {}
    names=list(INIT[ch])
    ties=dict(TIES.get(ch,[]))
    free=[n for n in names if n not in ties and n not in fixed]
    def build(pu):
        p=dict(pu)
        resolve(p,ties)
        return lcE.LETTERS[ch](p)
    p0={n:init[n]*sc(n) for n in names}
    p,g,st=fit(ch,build,p0,free,fixed={k:v*sc(k) for k,v in fixed.items()})
    res={n:p[n]/sc(n) for n in free}; res.update(fixed)
    resolve(res,ties)
    overlay(ch,g,SCR+f'fit_{ch}.png',scale=1.3)
    return res,st
if __name__=='__main__':
    allr=load()
    for ch in sys.argv[1]:
        r,st=run(ch); allr[ch]=r
        print(ch,{k:round(v,3) for k,v in r.items()})
    json.dump(allr,open(FN,'w'),indent=1)
def snap(r,q=32):
    return {k:(round(v*2)/2 if k.endswith('th') else round(v*q)/q) for k,v in r.items()}
def evaluate(ch,r,tag='snap'):
    """fit statistics of a fixed parameter set"""
    from lcmodel import scan_pts, boundary_pts, SPREAD
    from scipy.spatial import cKDTree
    g=lcE.LETTERS[ch]({k:v*sc(k) for k,v in r.items()})
    S=scan_pts(ch); M=boundary_pts(g.buffer(SPREAD,join_style=2))
    dS,_=cKDTree(M).query(S); dM,_=cKDTree(S).query(M); d=np.r_[dS,dM]
    overlay(ch,g,SCR+f'{tag}_{ch}.png',scale=1.3)
    return dict(med=np.median(d),p95=np.percentile(d,95),mx=max(dS.max(),dM.max()))
