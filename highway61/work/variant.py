"""Which version of the lower-case widths and spacing is being built (build1961.py <series> [variant]).

  ''        the main fonts: one lower-case width factor per series (lccondense.F_LC) and one spacing factor
            (spacing.LC_TIGHT), solved so that the whole a-z is as wide as the FHWA 2000 font's
  matched   one lower-case width factor and one spacing factor per series, solved so that the a-z ink widths and
            the a-z advances add up, in cap heights, to the FHWA 2000 Series <S> font's (s and w left out of both
            sums: FHWA draws them in other proportions); lcmatch.py solves them into matched.json.  E(M) too,
            against Series E Modified 2000, widened with its strokes as drawn (lccondense.PROP)
  natural   the lower case narrowed like the capitals, by the capitals' 1961 width factor (lccondense.F), and its
            spaces with it: one spacing factor per series (on the stroke-scaled 1966 table, as LC_TIGHT) keeps the
            a-z's spaces in E(M)'s proportion to its ink widths; lcnatural.py solves them into natural.json"""
import os, json
NAME = os.environ.get('H61_VARIANT', '')
NAMES = ('', 'matched', 'natural')
FAMILY = {'': '', 'matched': ' Matched', 'natural': ' Natural'}     # 'Highway 61 Revisited Matched B'
PS = {'': '', 'matched': 'Matched', 'natural': 'Natural'}          # 'Highway61RevisitedMatched-B'
JSON = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'matched.json')
NATURAL_JSON = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'natural.json')
_N = {}
def natural():
    """natural.json: {series: {'tight': spacing factor (as spacing.LC_TIGHT), ...}}"""
    if not _N and os.path.exists(NATURAL_JSON): _N.update(json.load(open(NATURAL_JSON)))
    return _N
_M = {}
def matched():
    """matched.json: {series: {'f': lower-case width factor, 'tight': spacing factor (as spacing.LC_TIGHT), ...}}"""
    if not _M: _M.update(json.load(open(JSON)))
    return _M
