"""The 'matched' version's lower-case width factor and spacing factor (variant.py):
python3 lcmatch.py [series...]  -> matched.json            (B C D E F EM; Series A follows Series B's)

For each FHWA 2000 Series <S> font, two totals in cap heights: the sum of the a-z ink widths (bounding boxes) and
the sum of the a-z advances.  The s and the w are left out of both sums: the FHWA narrow series draw them in other
proportions from the rest (in Series D the s is 0.71 of ModE's, the w 1.11, against about 0.9 for the others), and
this version keeps ours in proportion.

The width factor f (lccondense.series) is chosen so that our a-z ink widths, s and w left out, add up to FHWA's
(in cap heights); then the spacing factor (as spacing.LC_TIGHT: the 1966 table scaled by the stroke and by this) so
that our advances add up to FHWA's.  Both totals are close to linear, so secant steps settle them.  (The multipliers
over Series E Modified 2000 are recorded too.  Matching those ratios against our E(M) instead left every series about
3% narrow: our E(M) lower case is 3.3% narrower than ModE's.)  Needs the FHWA fonts (fhwa2otf.py)."""
import sys, os, json, string, importlib, fcntl, numpy as np
from fontTools.ttLib import TTFont
from fontTools.pens.boundsPen import BoundsPen
import lccondense, strokeop, spacing, variant
W = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(W, '..', '..')
LC = [c for c in string.ascii_lowercase if c not in 'sw']        # the letters measured


def metrics(fn):
    """{c: (ink width, advance)}"""
    f = TTFont(fn); g = f.getGlyphSet(); cm = f.getBestCmap(); out = {}
    for c in string.ascii_lowercase:
        n = cm[ord(c)]; p = BoundsPen(g); g[n].draw(p); out[c] = (p.bounds[2] - p.bounds[0], g[n].width)
    return out


FHWA = {'EM': 'ModE'}                     # E(M) is matched to Series E Modified 2000
def fhwa_totals(S):
    """(sum of ink widths, sum of advances) of FHWA Series S 2000's measured letters, in cap heights"""
    fn = os.path.join(W, '..', 'fhwa2000', 'Series%s2000-Normal.otf' % FHWA.get(S, S))
    f = TTFont(fn); g = f.getGlyphSet(); p = BoundsPen(g); g['H'].draw(p); cap = p.bounds[3]; X = metrics(fn)
    return sum(X[c][0] for c in LC) / cap, sum(X[c][1] for c in LC) / cap


def multipliers(S):
    """(ink-width multiplier, advance multiplier) of FHWA Series S 2000 over Series E Modified 2000"""
    X, M = (metrics(os.path.join(W, '..', 'fhwa2000', 'Series%s2000-Normal.otf' % s)) for s in (FHWA.get(S, S), 'ModE'))
    return (sum(X[c][0] for c in LC) / sum(M[c][0] for c in LC), sum(X[c][1] for c in LC) / sum(M[c][1] for c in LC))


def ink(S, f):
    out = lccondense.series(S, f, only=LC); tot = 0
    for c in LC:
        P = np.vstack(strokeop.to_polygons(out[c])); tot += P[:, 0].max() - P[:, 0].min()
    return tot


def save(S, d):
    with open(variant.JSON + '.lock', 'w') as lk:          # (several series can be solved at once)
        fcntl.flock(lk, fcntl.LOCK_EX)
        M = json.load(open(variant.JSON)) if os.path.exists(variant.JSON) else {}
        M.setdefault(S, {}).update(d)
        json.dump(M, open(variant.JSON, 'w'), indent=1, sort_keys=True)
    variant._M.clear()


def solve(S):
    mw, ma = multipliers(S); tw, ta = fhwa_totals(S)
    Tw = tw * 690; Ta = ta * 690                         # our cap height is 690
    print(S, 'FHWA totals (cap heights): ink width %.4f, advance %.4f' % (tw, ta), flush=True)
    f0 = lccondense.F_LC.get(S, 1.0); f1 = f0 * 1.05; w0 = ink(S, f0); w1 = ink(S, f1)
    for _ in range(3):
        f0, w0, f1 = f1, w1, f1 + (Tw - w1) * (f1 - f0) / (w1 - w0)
        w1 = ink(S, f1)
        if abs(w1 - Tw) < 0.5: break
    f = round(f1, 4); print(S, 'f %.4f: ink %.1f, target %.1f' % (f, ink(S, f), Tw), flush=True)
    save(S, {'f': f, 'ink_total': round(tw, 4), 'advance_total': round(ta, 4),
             'ink_multiplier': round(mw, 5), 'advance_multiplier': round(ma, 5)})
    # the spacing factor, with the build's own glyph boxes and starting bearings (build1961 imported, nothing written)
    sys.argv = ['build1961.py', S, 'matched']; sys.modules.pop('build1961', None); b = importlib.import_module('build1961')
    base = lccondense.stroke(S) / lccondense.V_EM
    def total(k):
        sb, _ = spacing.solve(S, b._cur, base * k, geo=b._geo)
        return sum(round(sb[c][0] + b._box[c][1] - b._box[c][0] + sb[c][1]) for c in LC)
    k0 = spacing.LC_TIGHT.get(S, 1.0); k1 = k0 * 0.8; t0, t1 = total(k0), total(k1)
    for _ in range(5):
        if abs(t1 - Ta) < 1 or t1 == t0: break
        k0, t0, k1 = k1, t1, k1 + (Ta - t1) * (k1 - k0) / (t1 - t0)
        t1 = total(k1)
    k = round(k1, 4); print(S, 'tight %.4f: advances %d, target %.1f' % (k, total(k), Ta), flush=True)
    save(S, {'tight': k})


if __name__ == '__main__':
    for S in sys.argv[1:] or 'BCDEF': solve(S)
