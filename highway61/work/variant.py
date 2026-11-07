"""Which version of the lower-case widths and spacing is being built (build1961.py <series> [variant]).

  ''        the main fonts: one lower-case width factor per series (lccondense.F_LC) and one spacing factor
            (spacing.LC_TIGHT), solved so that the whole a-z is as wide as the FHWA 2000 font's
  matched   one lower-case width factor and one spacing factor per series, solved so that the a-z ink widths and
            the a-z advances add up, in cap heights, to the FHWA 2000 Series <S> font's (s and w left out of both
            sums: FHWA draws them in other proportions); lcmatch.py solves them into matched.json.  E(M) too,
            against Series E Modified 2000, widened with its strokes as drawn (lccondense.PROP)
(A 'natural' version, the lower case narrowed by the capitals' width factor, was tried and dropped: HANDOFF.md)"""
import os, json
NAME = os.environ.get('H61_VARIANT', '')
NAMES = ('', 'matched')
FAMILY = {'': '', 'matched': ' Matched'}     # 'Highway 61 Revisited Matched B'
PS = {'': '', 'matched': 'Matched'}          # 'Highway61RevisitedMatched-B'
JSON = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'matched.json')
_M = {}
def matched():
    """matched.json: {series: {'f': lower-case width factor, 'tight': spacing factor (as spacing.LC_TIGHT), ...}}"""
    if not _M: _M.update(json.load(open(JSON)))
    return _M
