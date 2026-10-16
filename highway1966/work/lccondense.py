"""Condensed (B, C, D) and expanded (F) lower case from the E(M) constructions, by the rule the 1961 capitals follow:
each series keeps the Series E core (width less stroke) scaled by f = 1 - 0.175 k (B 0.475, C 0.65, D 0.825, F 1.175)
with its own stroke, which steps by 690/64 per series (capitals B 86.25 ... F 129.38); the lower-case strokes step
the same way from the Series E lower case.  Both vertical and horizontal strokes keep the full weight; terminal cuts
and slanted stem ends keep their angles.  series(S) -> {letter: segment contours} (strokeop format)."""
import numpy as np, lcE, strokeop
SQ=690/8; X=6*SQ
import importlib
A={S:importlib.import_module('dim%s1'%S).A for S in 'BCDEF'}     # the 1961 table strokes
F={'B':0.475,'C':0.65,'D':0.825,'E':1.0,'F':1.175}
V_EM=(1+9/16)*SQ                       # E(M) lower-case stem, the source of the constructions
D_EM=(690/5-A['E'])/2                  # E(M) <-> E stroke offset
T_CORE=45.0                            # offset to the thin core (leaves about 45 units of every stroke)
def stroke(S): return V_EM-2*D_EM+(A[S]-A['E'])
def pair(c):
    """(cut letter, its ellipses, uncut letter, its ellipses).  Straight-line letters: the uncut letter has its
    slanted stem top (k) squared off"""
    from shapely.geometry import Polygon, box
    if c in lcE.dims():
        gc,pc=lcE.build(c); gu,pu=lcE.build(c,True); return gc,pc,gu,pu
    V=np.array(lcE.polys()[c])*SQ; gc=Polygon(V); gu=gc
    for a,b in zip(V,np.roll(V,-1,axis=0)):
        if min(a[1],b[1])>X+10 and 0.05<abs((b[1]-a[1])/(b[0]-a[0]+1e-9))<2:
            gu=gu|box(min(a[0],b[0]),min(a[1],b[1])-1,max(a[0],b[0]),max(a[1],b[1]))
    return gc,[],lcE.clean(gu),[]
def cuts(c):
    """the straight cuts of a letter: lines on the cut letter's outline that lie inside the uncut letter.
    Returns [(p0, p1, n)], n the unit normal pointing into the removed part"""
    import lcoutline
    from shapely.geometry import Point
    gc,pc,gu,_=pair(c)
    out=[]
    for segs in lcoutline.contours(gc,pc):
        for sg in segs:
            if sg[0]!='L': continue
            p0,p1=np.array(sg[1]),np.array(sg[2]); m=(p0+p1)/2
            if gu.boundary.distance(Point(m))<0.05: continue
            u=(p1-p0)/np.linalg.norm(p1-p0); n=np.array([-u[1],u[0]])
            if not (gu.contains(Point(m+n)) and not gc.contains(Point(m+n))): n=-n
            out.append((p0,p1,n))
    return out
def _polygon(C):
    from shapely.geometry import Polygon
    rings=strokeop.to_polygons(C,step=np.radians(0.25))
    g=Polygon()
    for r in rings:
        if len(r)>2: g=g.symmetric_difference(Polygon(r).buffer(0))
    return g
def series(S):
    f=F[S]; T=stroke(S)
    gx=(T-V_EM*f)/2; ey=(V_EM-T)/2              # horizontal offset after scaling; vertical erosion
    t=ey                                        # stub length scale
    k=(X+2*ey)/X
    def ymap(y): return y-ey if y<0 else (k*y-ey if y<X else y+ey)
    P=lambda q: np.array([q[0]*f,ymap(q[1])])
    def edge(q,n):                             # a point on the outline with outward normal n
        return P(q)+strokeop._support(n,(gx,0.0))+strokeop._support(n,(0.0,-ey))
    import lcoutline
    from shapely.geometry import Polygon
    out={}
    for c,(g,prims) in lcE.series().items():
        if c not in lcE.dims() and c!='k':              # straight-line letters without cuts: condensed directly
            out[c]=strokeop.condense2(strokeop.from_lc(g,prims),f,gx,ey,X,dot_scale=T/V_EM); continue
        gc,pc,gu,pu=pair(c)
        cs=cuts(c); gl=gc; stubs=[]; DX={}
        for p0,p1,n in cs:
            u=(p1-p0)/np.linalg.norm(p1-p0); L=np.linalg.norm(p1-p0)
            Dx=0.3*L+abs(gx)+ey+10                   # covers the turn of the cut back to its angle, and the offsets
            for _ in range(8):                       # shorter if it would reach another part of the letter
                st=(gu&Polygon([p0-u*2,p1+u*2,p1+u*2+n*Dx,p0-u*2+n*Dx])).difference(gc)
                gap=2*max(gx,0)/f+6                  # clearance that survives the scaling and the growth
                if st.intersection(gc.buffer(gap)).area<1.15*L*gap+30: break
                Dx*=0.75
            DX[(tuple(p0),tuple(p1))]=Dx
            stubs.append(st); gl=gl|st
        gl=lcE.clean(gl)
        C=strokeop.condense2(strokeop.from_lc(gl,pc+pu),f,gx,ey,X,mitre=4.0,dot_scale=T/V_EM)
        geom=_polygon(C)
        prims2=[s[3] for cc in C for s in cc if s[0]=='E']
        for (p0,p1,n),st in zip(cs,stubs):
            # the cut keeps its angle and passes through the mapped edge point at the middle of the original cut; what
            # is removed beyond it is the condensed image of the stub, condensed on its own together with a band of
            # the letter on the near side of the cut (so that its edge at the cut does not move)
            u=(p1-p0)/np.linalg.norm(p1-p0); m=edge((p0+p1)/2,n)
            Dx=DX[(tuple(p0),tuple(p1))]
            band=gl&Polygon([p0-u*2-n*Dx,p1+u*2-n*Dx,p1+u*2+n*Dx,p0-u*2+n*Dx])
            band=lcE.clean(band)
            if band.is_empty: continue
            Cs=strokeop.condense2(strokeop.from_lc(band,pc+pu),f,gx,ey,X,mitre=4.0)
            Sg=_polygon(Cs).buffer(3.0,join_style=2)
            H=Polygon([m-u*3000,m+u*3000,m+u*3000+n*3000,m-u*3000+n*3000])
            geom=geom.difference(H&Sg)
            # ears: whatever is left beyond the cut close to its two ends (found exactly on the output's outline)
            V=[np.array(q) for pg in getattr(geom,'geoms',[geom]) for r in [pg.exterior,*pg.interiors] for q in r.coords]
            half=np.linalg.norm(p1-p0)*max(f,1)/2+40
            V=[q for q in V if abs((q-m)@n)<0.05 and abs((q-m)@u)<half]
            if len(V)>=2:
                tt=[(q-m)@u for q in V]
                for e in (V[int(np.argmin(tt))],V[int(np.argmax(tt))]):
                    E=Polygon([e-u*15,e+u*15,e+u*15+n*15,e-u*15+n*15])
                    geom=geom.difference(E)
        geom=lcE.clean(geom)
        if geom.geom_type=='MultiPolygon':            # stray bits of a stub (a dot is far larger)
            from shapely.geometry import MultiPolygon
            ps=[q for q in geom.geoms if q.area>1000]; geom=ps[0] if len(ps)==1 else MultiPolygon(ps)
        out[c]=[[tuple(sg) for sg in segs] for segs in lcoutline.contours(geom,prims2)]
    return out
if __name__=='__main__':
    for S in 'BCDF': print(S,'f',F[S],'stroke',round(stroke(S),2))
