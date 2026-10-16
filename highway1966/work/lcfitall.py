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

 'l':dict(w=1.57,ytop=8.0,s=0.4),
 'i':dict(x0=0.356,w=1.57,dx=1.14,dy=7.17,dr=0.78),
 'h':dict(w=1.57,x2=3.37,oL=0.97,oyl=3.96,oyr=3.96,oxt=2.75,T=6.11,ccx=2.55,ccy=3.44,ca=1.02,cb=1.05,ytop=8.0,s=0.4),
 'm':dict(w=1.57,x2=3.39,x3=6.83,oL=0.97,oyl=3.96,oyr=3.96,oxt=2.75,T=6.11,ccx=2.55,ccy=3.44,ca=1.02,cb=1.05,
          roL=0.97,royl=3.96,royr=3.96,roxt=2.75,rT=6.11,rccx=2.55,rccy=3.44,rca=1.02,rcb=1.05),
 'u':dict(w=1.57,x2=3.45,oL=0.97,oyl=3.96,oyr=3.96,oxt=2.75,T=6.11,ccx=2.55,ccy=3.44,ca=1.02,cb=1.05),
 'r':dict(w=1.57,x2=3.4,oL=0.97,oyl=3.96,oyr=3.96,oxt=2.75,T=6.11,ccx=2.55,ccy=3.44,ca=1.02,cb=1.05,tx=3.75,ty=5.0,tth=75.0),

 'd':dict(xs=3.357,xr=5.0,ytop=8.0,s=0.4,BL=-0.165,Byl=2.844,BR=4.615,Byr=3.323,BR2=4.941,Byr2=1.693,Bxt=1.994,BT=6.057,Bxb=2.039,BB=-0.073,cx=2.496,cy=2.994,ca=1.092,cb=1.529),
 'q':dict(xs=3.357,xr=5.0,ybot=-2.35,sb=0.364,BL=-0.165,Byl=2.844,BR=4.615,Byr=3.323,BR2=4.941,Byr2=1.693,Bxt=1.994,BT=6.057,Bxb=2.039,BB=-0.073,cx=2.496,cy=2.994,ca=1.092,cb=1.529),
 'b':dict(w=1.56,ytop=8.0,s=0.4,BR=5.05,Byr=2.85,BL=0.39,Byl=3.32,BL2=0.06,Byl2=1.69,Bxt=3.0,BT=6.06,Bxb=2.96,BB=-0.07,cx=2.5,cy=3.0,ca=1.09,cb=1.53),
 'p':dict(w=1.56,ybot=-2.35,sb=0.354,BR=5.05,Byr=2.85,BL=0.39,Byl=3.32,BL2=0.06,Byl2=1.69,Bxt=3.0,BT=6.06,Bxb=2.96,BB=-0.07,cx=2.5,cy=3.0,ca=1.09,cb=1.53),

 'c':dict(OL=-0.05,Oyl=3.0,OR=5.2,Oyr=2.9,Oxt=2.75,OT=6.1,Oxb=2.8,OB=-0.07,ix=2.69,iy=3.04,ia=1.14,ib=1.6,
          ux=4.42,uy=4.42,uth=45.0,lx=4.6,ly=1.66,lth=131.5),

 't':dict(x0=1.0,x1=2.53,ytop=8.0,s=0.4,bx0=0.0,bx1=3.92,by=4.53,hx=2.62,hB=-0.02,hyl=1.5,hR=3.95,hyr=1.6,
          icy=2.2,ia=0.55,ib=0.72,tx=3.75,ty=1.67,tth=8.0),
 'f':dict(x0=1.45,x1=3.01,bx0=0.76,bx1=4.0,by=4.39,hx=2.94,hB=-2.0,hyl=-0.26,hR=4.03,hyr=-1.42,
          icy=-0.2,ia=0.5,ib=0.5,tx=4.0,ty=-1.0,tth=90.0),
 'j':dict(x0=2.0,x1=3.54,hx=4.24,hB=-2.93,hyl=0.0,hR=4.9,hyr=-1.6,icy=0.0,ia=1.25,ib=1.21,
          tx=4.84,ty=-1.8,tth=90.0,dx=2.733,dy=7.185,dr=0.85),

 's':dict(UL=-0.012,Uyl=4.174,UR=4.870,Uyr=3.478,Uxt=1.971,UT=6.133,Uxb=3.478,UB=1.913,LL=0.000,Lyl=1.913,LR=4.904,Lyr=1.855,Lxt=1.449,LT=4.174,Lxb=2.667,LB=-0.081,uix=3.014,uiy=3.826,uia=1.357,uib=0.754,lix=1.913,liy=2.203,lia=1.345,lib=0.835,utx=4.325,uty=4.858,ltx=0.197,lty=1.229,utth=65.000,ltth=67.000),
}
TIES={'o':[('Oxb','Oxt'),('Ixb','Ixt'),('Ixt','Oxt')]}
def sc(n): return 1.0 if (n.endswith('th') or n in ('s','sb')) else SQ
def resolve(p,ties):
    for _ in range(5):
        for a,b in ties.items():
            if b in p: p[a]=p[b]
def run(ch,snap=False,fixed=None):
    init=load().get(ch,INIT[ch]); fixed=fixed or {}
    names=list(INIT[ch])+['ox']
    init=dict(init); init.setdefault('ox',0.0)
    ties=dict(TIES.get(ch,[]))
    free=[n for n in names if n not in ties and n not in fixed]
    def build(pu):
        p=dict(pu)
        resolve(p,ties)
        from shapely import affinity
        return affinity.translate(lcE.LETTERS[ch](p),p.get('ox',0.0))
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
    return {k:(round(v*2)/2 if k.endswith('th') else round(v*1000)/1000 if k in ('s','sb') else round(v*q)/q) for k,v in r.items()}
def evaluate(ch,r,tag='snap'):
    """fit statistics of a fixed parameter set"""
    from lcmodel import scan_pts, boundary_pts, SPREAD
    from scipy.spatial import cKDTree
    from shapely import affinity
    g=affinity.translate(lcE.LETTERS[ch]({k:v*sc(k) for k,v in r.items()}),r.get('ox',0)*SQ)
    S=scan_pts(ch); M=boundary_pts(g.buffer(SPREAD,join_style=2))
    dS,_=cKDTree(M).query(S); dM,_=cKDTree(S).query(M); d=np.r_[dS,dM]
    overlay(ch,g,SCR+f'{tag}_{ch}.png',scale=1.3)
    return dict(med=np.median(d),p95=np.percentile(d,95),mx=max(dS.max(),dM.max()))
