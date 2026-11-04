"""The 'natural' version's lower-case spacing factor (variant.py): python3 lcnatural.py [series...] -> natural.json

The Natural lower case is narrowed by the capitals' 1961 width factor, and its spaces go with it: the a-z's total
space (advances less ink widths) is in the same proportion to its total ink width as E(M)'s
(otf/Highway61Revisited-EM.otf, 0.44).  (Scaling the 1966 table by the stroke alone, as the capitals' spaces are,
had left the narrow series' lower case much looser than E(M)'s: Series B 0.56.)  The spacing factor multiplies the
stroke-scaled 1966 table, as spacing.LC_TIGHT does; the total is linear in it, so secant steps settle it."""
import sys, os, json, string, importlib, fcntl
import lccondense, spacing, variant, lcmatch
W = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(W, '..', '..'); LC = string.ascii_lowercase


def save(S, d):
    with open(variant.NATURAL_JSON + '.lock', 'w') as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        M = json.load(open(variant.NATURAL_JSON)) if os.path.exists(variant.NATURAL_JSON) else {}
        M.setdefault(S, {}).update(d)
        json.dump(M, open(variant.NATURAL_JSON, 'w'), indent=1, sort_keys=True)
    variant._N.clear()
    try: os.remove(variant.NATURAL_JSON + '.lock')
    except OSError: pass


def solve(S):
    E = lcmatch.metrics(os.path.join(ROOT, 'otf', 'Highway61Revisited-EM.otf'))
    ratio = sum(E[c][1] - E[c][0] for c in LC) / sum(E[c][0] for c in LC)
    sys.argv = ['build1961.py', S, 'natural']; sys.modules.pop('build1961', None); b = importlib.import_module('build1961')
    ink = sum(b._box[c][1] - b._box[c][0] for c in LC); Ta = ink * (1 + ratio)
    base = lccondense.stroke(S) / lccondense.V_EM
    def total(k):
        sb, _ = spacing.solve(S, b._cur, base * k, geo=b._geo)
        return sum(round(sb[c][0] + b._box[c][1] - b._box[c][0] + sb[c][1]) for c in LC)
    k0, k1 = 1.0, 0.8; t0, t1 = total(k0), total(k1)
    for _ in range(5):
        if abs(t1 - Ta) < 1 or t1 == t0: break
        k0, t0, k1 = k1, t1, k1 + (Ta - t1) * (k1 - k0) / (t1 - t0)
        t1 = total(k1)
    k = round(k1, 4); t = total(k)
    print(S, 'E(M) space/ink %.4f; tight %.4f: advances %d, target %.1f, space/ink %.4f' % (ratio, k, t, Ta, (t - ink) / ink),
          flush=True)
    save(S, {'tight': k})


if __name__ == '__main__':
    for S in sys.argv[1:] or 'ABCDEF': solve(S)
