"""Spacing by the standard tables (spacingcheck.py has them): side bearings and class kerning such that every pair
the tables cover is spaced exactly as they say.

The tables give the space between extreme edges from the preceding character and the class of the following one
only, so each character gets its own right side bearing and every character of a following class the same left
side bearing; least squares chooses them (pulled weakly towards the current ones, which fixes the split between
left and right), and the remainder for each (preceding character, following class) becomes a kerning value.

Capitals and figures: Table I/II codes, space = factor x series stroke (E(M): Series E's spaces).  Lower case: the
1966 E(M) table (capitals and lower case preceding), scaled for the other series by their lower-case stroke over
E(M)'s and tightened to the FHWA 2000 alphabet widths (LC_TIGHT).  Pairs the tables do not cover (lower case before a capital, letters and figures, punctuation) are left to
the side bearings."""
import numpy as np
import spacingcheck as sc
NAME=sc.NAME
def classes():
    """(preceding characters, following classes as strings) for the three tables"""
    cap_f=['BDEFHIKLMNPRU','CGOQSXZ','AJTVWY']; num_f=['15','236890','47']
    return cap_f,num_f,sc.LCCOLS
# Lower-case spaces relative to the E(M) table scaled by lower-case stroke, solved by lcspace.py so that each
# series' a-z advances add up to the FHWA 2000 font's (as the lower-case widths are, lcwidth.py).  E(M) is the
# table itself; Series A (no FHWA font) takes Series B's factor.
LC_TIGHT={'B':0.574,'C':0.587,'D':0.691,'E':0.473,'F':0.765}
def lc_scale(S):
    """lower-case table scale for series S: E(M) 1, others their lower-case stroke over E(M)'s times LC_TIGHT"""
    import lccondense
    if S=='EM': return 1.0
    return lccondense.stroke(S)/lccondense.V_EM*LC_TIGHT.get('B' if S=='A' else S,1.0)
_JOV_EM=[]
def jov_em():
    """how far E(M)'s j's tail reaches left of its stem (units; 102.4)"""
    if not _JOV_EM:
        import lcE
        from shapely.geometry import box
        g,_=lcE.series()['j']; _JOV_EM.append((g&box(-1e4,50,1e4,1e4)).bounds[0]-g.bounds[0])
    return _JOV_EM[0]
def geometry(poly,chars):
    """what the j column's spaces need from a font's glyphs (poly(char) -> shapely geometry): how far its j's tail
    reaches left of the stem, and for each character with a descender, how far that stops short of its right edge"""
    from shapely.geometry import box
    g=poly('j'); jov=(g&box(-1e4,50,1e4,1e4)).bounds[0]-g.bounds[0]; desc={}
    for ch in chars:
        g=poly(ch); lo=g&box(-1e4,-1e4,1e4,-20)
        if not lo.is_empty: desc[ch]=g.bounds[2]-lo.bounds[2]
    return {'jov':jov,'desc':desc}
def targets(S,lc=None,geo=None):
    """{(preceding char, following class name): space in font units}; lc overrides the lower-case scale.  geo
    (geometry()): the j column is spaced by E(M)'s stem-to-stem gap, scaled.  The table gives the space to the j's
    extreme edge, the end of its tail, and the tail does not shrink with the spaces (Series E's reaches as far as
    E(M)'s, whose spaces it has at 0.4), so scaling the space alone left a gap before every j (i-j 1.55 of i-i in
    E, against E(M)'s 0.96).  Characters with a descender keep at least the scaled table space at its level."""
    cap_f,num_f,lc_f=classes(); st=sc.STROKE[S]; out={}
    for a,codes in sc.CODES.items():
        for g,cls in enumerate(cap_f): out[(a,cls)]=sc.FACTOR[codes[g]]*st
    for a,codes in sc.NCODES.items():
        for g,cls in enumerate(num_f): out[(a,cls)]=sc.FACTOR[codes[g]]*st
    r=lc_scale(S) if lc is None else lc
    for row,vals in zip(sc.LCROWS,sc.LCRATIO):
        for cls,v in zip(lc_f,vals):
            for a in row:
                t=v*sc.XH*r
                if cls=='j' and geo is not None:
                    t2=(v*sc.XH+jov_em())*r-geo['jov']
                    if a in geo['desc']: t2=max(t2,t-geo['desc'][a])
                    t=t2
                out[(a,cls)]=t
    return out
def solve(S,cur,lc=None,geo=None):
    """cur: {char: (lsb, rsb)} now; geo: geometry().  Returns ({char: (lsb, rsb)}, {(char, class): kern})"""
    T=targets(S,lc,geo); pre=sorted({a for a,_ in T}); cls=sorted({c for _,c in T})
    ip={a:i for i,a in enumerate(pre)}; ic={c:len(pre)+i for i,c in enumerate(cls)}; n=len(pre)+len(cls)
    rows=[]; y=[]
    for (a,c),v in T.items():
        r=np.zeros(n); r[ip[a]]=1; r[ic[c]]=1; rows.append(r); y.append(v)
    w=0.05                                          # weak pull towards the current bearings
    for a in pre:
        r=np.zeros(n); r[ip[a]]=w; rows.append(r); y.append(w*cur[a][1])
    for c in cls:
        r=np.zeros(n); r[ic[c]]=w; rows.append(r); y.append(w*np.mean([cur[ch][0] for ch in c]))
    x,*_=np.linalg.lstsq(np.array(rows),np.array(y),rcond=None)
    rsb={a:x[ip[a]] for a in pre}; lsb={ch:x[ic[c]] for c in cls for ch in c}
    sb={ch:(round(lsb.get(ch,cur[ch][0])),round(rsb.get(ch,cur[ch][1]))) for ch in cur}
    kern={}
    for (a,c),v in T.items():
        k=round(v-sb[a][1]-sb[c[0]][0])
        if k: kern[(a,c)]=k
    return sb,kern
def fea(kern,gname):
    """feature file text: one class per following class, a pair rule per (preceding glyph, class)"""
    cls=sorted({c for _,c in kern}); cn={c:'@F%d'%i for i,c in enumerate(cls)}
    lines=['%s = [%s];'%(cn[c],' '.join(gname(ch) for ch in c)) for c in cls]
    lines.append('feature kern {')
    for (a,c),k in sorted(kern.items()): lines.append('    pos %s %s %d;'%(gname(a),cn[c],k))
    lines.append('} kern;')
    return '\n'.join(lines)
