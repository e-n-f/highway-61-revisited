"""dimension tables from the 1961 pages: parse rows (inches) and cross-check rows of different letter heights"""
from fractions import Fraction as Fr
def inch(s): return float(sum(Fr(p) for p in s.split()))
def units(row, h): return {k: inch(v) * 690.0 / h for k, v in row.items()}
def check(rows, tol=0.75):
    """rows: {height: {letter: 'inches'}}; report letters whose rows disagree by more than tol font units"""
    out = []
    hs = sorted(rows)
    for k in rows[hs[-1]]:
        vals = {h: inch(rows[h][k]) * 690.0 / h for h in hs if k in rows[h]}
        spread = max(vals.values()) - min(vals.values())
        if spread > tol: out.append((k, {h: round(v, 1) for h, v in vals.items()}))
    return out
