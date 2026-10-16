"""Spacing by the standard tables (spacingcheck.py has them): side bearings and class kerning such that every pair
the tables cover is spaced exactly as they say.

The tables give the space between extreme edges from the preceding character and the class of the following one
only, so each character gets its own right side bearing and every character of a following class the same left
side bearing; least squares chooses them (pulled weakly towards the current ones, which fixes the split between
left and right), and the remainder for each (preceding character, following class) becomes a kerning value.

Capitals and figures: Table I/II codes, space = factor x series stroke (E(M): Series E's spaces).  Lower case: the
1983 E(M) table (capitals and lower case preceding), scaled for the other series by their lower-case stroke over
E(M)'s.  Pairs the tables do not cover (lower case before a capital, letters and figures, punctuation) are left to
the side bearings."""
import numpy as np
import spacingcheck as sc
NAME=sc.NAME
def classes():
    """(preceding characters, following classes as strings) for the three tables"""
    cap_f=['BDEFHIKLMNPRU','CGOQSXZ','AJTVWY']; num_f=['15','236890','47']
    return cap_f,num_f,sc.LCCOLS
def targets(S):
    """{(preceding char, following class name): space in font units}"""
    import lccondense
    cap_f,num_f,lc_f=classes(); st=sc.STROKE[S]; out={}
    for a,codes in sc.CODES.items():
        for g,cls in enumerate(cap_f): out[(a,cls)]=sc.FACTOR[codes[g]]*st
    for a,codes in sc.NCODES.items():
        for g,cls in enumerate(num_f): out[(a,cls)]=sc.FACTOR[codes[g]]*st
    r=1.0 if S=='EM' else lccondense.stroke(S)/lccondense.V_EM
    for row,vals in zip(sc.LCROWS,sc.LC375):
        for cls,v in zip(lc_f,vals):
            for a in row: out[(a,cls)]=v*sc.XH/375*r
    return out
def solve(S,cur):
    """cur: {char: (lsb, rsb)} now.  Returns ({char: (lsb, rsb)}, {(char, class): kern})"""
    T=targets(S); pre=sorted({a for a,_ in T}); cls=sorted({c for _,c in T})
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
