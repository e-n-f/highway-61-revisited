"""Series F capitals and figures from the 1961 dimensioned drawings (pages F-1..F-6), as exact vector regions."""
import dimF1 as d1, dimF2 as d2, dimF3 as d3, dimF4 as d4, dimF5 as d5, dimF6 as d6
def _g(fn):
    def f():
        r = fn(); r = r[0] if isinstance(r, tuple) else r
        return r.geom if hasattr(r, 'geom') else r
    return f
REG = {}
for m in (d1, d2, d3, d4):
    for ch in m.__dict__:
        if ch.startswith('letter_'): REG[ch[7:]] = _g(getattr(m, ch))
for ch in 'YZ': REG[ch] = _g(getattr(d5, 'letter_' + ch))
for ch in '1234': REG[ch] = _g(getattr(d5, 'fig' + ch))
for ch in '567890': REG[ch] = _g(getattr(d6, 'fig' + ch))
INFO = {'F-%d' % (i + 1): m.CHECK for i, m in enumerate((d1, d2, d3, d4, d5, d6))}
