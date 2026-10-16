"""Condensed (B, C, D) and expanded (F) lower case from the E(M) constructions, by the rule the 1961 capitals follow:
each series keeps the Series E core (width less stroke) scaled by f = 1 - 0.175 k (B 0.475, C 0.65, D 0.825, F 1.175)
with its own stroke, which steps by 690/64 per series (capitals B 86.25 ... F 129.38); the lower-case strokes step
the same way from the Series E lower case.  Both vertical and horizontal strokes keep the full weight; terminal cuts
and slanted stem ends keep their angles.  series(S) -> {letter: segment contours} (strokeop format)."""
import numpy as np, lcE, strokeop
SQ=690/8; X=6*SQ
import importlib
A={S:importlib.import_module('dim%s1'%S).A for S in 'ABCDEF'}     # the 1961 table strokes
F={'A':0.450,'B':0.475,'C':0.65,'D':0.825,'E':1.0,'F':1.175}
F_LC={'A':0.521,'B':0.55,'C':0.82,'D':0.919,'E':1.0,'F':1.124}
                                       # lower-case width factors: the whole a-z as wide as the FHWA 2000 font's,
                                       # in cap heights (lcwidth.py); the capitals keep the 1961 factors F
                                       # Series A (no FHWA font): B's ratio F_LC/F, 0.55/0.475, times A's 0.450
                                       # Series E: 1, E(M)'s forms exactly (as the capitals); its spacing
                                       # (spacing.LC_TIGHT) makes up the alphabet width instead
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
TIDY=4.0
def _tidy(geom):
    """opening with mitred offsets of TIDY: removes hooks and slivers narrower than 2 TIDY, keeps corners (mitre limit
    high enough for the acute ones) and arcs (concentric offsets) where they are.  (A closing as well would mend
    nicks, but it fills in beside acute terminal tips and would soften the notches meant to be there.)"""
    from shapely.geometry import Polygon, MultiPolygon
    k=dict(join_style=2,mitre_limit=100.0)
    def keep(g):
        ps=[p for p in getattr(g,'geoms',[g]) if p.area>50]
        return ps[0] if len(ps)==1 else MultiPolygon(ps)
    g=geom.buffer(-TIDY,**k).buffer(TIDY,**k)
    return keep(g)
S_SPINE={'A':-9.0,'B':-12.0}   # the s spine's tangent points moved towards the lobes' extremes (units): the most
                                # condensed spines come out steeper and thinner (narrowest 0.77, 0.76 of the stroke
                                # against 0.81-0.83 elsewhere); this makes them 0.89, 0.87
def _s_direct(f,gx,ey,k,S=None):
    """the s condensed by mapping its construction rather than its outline: ellipses bounding ink widen by gx and
    lose ey in height, those bounding counters the opposite; x scales by f and heights by k (less ey), as in
    condense2; the spine's tangent points scale; each terminal line keeps its angle and passes through the mapped
    middle of the original cut, moved out by the offsets.  The construction then re-solves its tangencies exactly
    (offsetting the outline would part the non-axis joins: an arc moves by changing its semi-axes)."""
    import lcgeom, lcoutline
    d=lcE.dims()['s']; p={kk:(v if kk.endswith('th') else v*SQ) for kk,v in d.items()}
    ink=lambda e: (e[0]*f,e[1]*k-ey,e[2]*f+gx,e[3]*k-ey)
    hole=lambda e: (e[0]*f,e[1]*k-ey,e[2]*f-gx,e[3]*k+ey)
    u1,u2,um,l1,l2,lm=lcE.s_ellipses(p)
    q=dict(p); q['ell']=(ink(u1),ink(u2),hole(um),ink(l1),ink(l2),hole(lm))
    q['C1x'],q['C1a']=p['C1x']*f,p['C1a']*f-gx; q['C2x'],q['C2a']=p['C2x']*f,p['C2a']*f-gx
    q['s1'],q['s2']=p['s1']*f-S_SPINE.get(S,0.0),p['s2']*f+S_SPINE.get(S,0.0)
    for p0,p1,n in cuts('s'):
        m=(p0+p1)/2; pre='ut' if m[1]>X/2 else 'lt'
        mm=np.array([m[0]*f,m[1]*k-ey])+strokeop._support(n,(gx,0.0))+strokeop._support(n,(0.0,-ey))
        q[pre+'x'],q[pre+'y']=mm; q[pre+'ref']=tuple(mm-n*30.0)      # a point on the ink side of the cut
    lcgeom.REG.clear(); g=lcE.clean(lcE.letter_s_ell(q)); prims=list(lcgeom.REG)
    return [[tuple(sg) for sg in segs] for segs in lcoutline.contours(g,prims)]
def series(S,f=None,only=None):
    f=F_LC[S] if f is None else f; T=stroke(S)
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
        if only and c not in only: continue
        if c=='s':                                      # built from its mapped construction
            out[c]=_s_direct(f,gx,ey,k,S); continue
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
                near=Polygon([p0-u*15-n*15,p1+u*15-n*15,p1+u*15+n*15,p0-u*15+n*15])
                clear=st.distance(gc.difference(near)) if not st.is_empty else 1e9    # e.g. the e's bar
                if st.intersection(gc.buffer(gap)).area<1.15*L*gap+30 and clear>=3: break
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
            # ears: whatever is left beyond the cut close to its two ends (found exactly on the output's outline): the
            # ends of the straight run of outline along the cut line nearest the cut's middle (another part of the
            # letter that merely crosses the line further on, like the inside of a wider t's foot, is not an ear)
            half=np.linalg.norm(p1-p0)*max(f,1)/2+40
            runs=[]
            for pg in getattr(geom,'geoms',[geom]):
                for r in [pg.exterior,*pg.interiors]:
                    Q=np.array(r.coords)
                    for a,b in zip(Q[:-1],Q[1:]):
                        if abs((a-m)@n)<0.05 and abs((b-m)@n)<0.05:
                            ta,tb=sorted(((a-m)@u,(b-m)@u))
                            if tb-ta>1e-6 and ta<half and tb>-half: runs.append([ta,tb])
            runs.sort(); merged=[]
            for ta,tb in runs:
                if merged and ta<=merged[-1][1]+0.05: merged[-1][1]=max(merged[-1][1],tb)
                else: merged.append([ta,tb])
            if merged:
                ta,tb=min(merged,key=lambda iv: 0 if iv[0]<=0<=iv[1] else min(abs(iv[0]),abs(iv[1])))
                for e in (m+u*ta,m+u*tb):
                    E=Polygon([e-u*15,e+u*15,e+u*15+n*15,e-u*15+n*15])
                    geom=geom.difference(E)
        geom=lcE.clean(geom)
        if geom.geom_type=='MultiPolygon':            # stray bits of a stub: keep the body and the dot (above X)
            from shapely.geometry import MultiPolygon
            body=max(geom.geoms,key=lambda q:q.area)
            ps=[body]+[q for q in geom.geoms if q is not body and c in 'ij' and q.bounds[1]>X]
            geom=ps[0] if len(ps)==1 else MultiPolygon(ps)
        if gx<0:                                     # the outline was eroded sideways (Series A): its re-cuts and
            geom=_tidy(geom)                         # junctions leave nicks and hooks a few units wide
        out[c]=[[tuple(sg) for sg in segs] for segs in lcoutline.contours(geom,prims2)]
    return out
if __name__=='__main__':
    for S in 'BCDF': print(S,'f',F[S],'stroke',round(stroke(S),2))
