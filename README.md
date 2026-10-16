# standard-alphabets

Standard Alphabets for Highway Signs: **Highway 61 Revisited**, a family of digital fonts reconstructed from the
U.S. Bureau of Public Roads' *Standard Alphabets for Highway Signs*. The capitals and figures come from the 1961
dimensioned drawings, and the lower case from the 1966 grid sheets. Each is rebuilt as exact geometry rather than
traced. The goal is to improve on the official digital fonts, not to copy them.

| Font | Files |
|---|---|
| Highway 61 Revisited B, C, D, E, F | `highway1966/outputs/Highway61Revisited-{B,C,D,E,F}.otf` / `.ttf` |
| Highway 61 Revisited E(M) | `highway1966/outputs/Highway61Revisited-EM.otf` / `.ttf` |

The `.otf` files are the exact versions: cubic Béziers with fractional coordinates. The `.ttf` files are
quadratic with integer coordinates. Each font has A–Z, a–z, 0–9, `. , ! ' " ‘ ’ “ ” -` and the space.

Font units: 1000 per em. Cap height 690 (an 8-inch letter, so 1 inch = 86.25 units), x-height 517.5, ascender 690,
descenders about −200 to −260.

Specimens: `highway1966/outputs/specimen-lc-series.png` (all six series) and
`highway1966/outputs/comparison-1966-vs-2000.png` (against the FHWA 2000 fonts). The detailed working history, with
every construction decision, is in `highway1966/HANDOFF.md`.

## How the letterforms are represented

### Capitals and figures: the 1961 dimensioned drawings

The 1961 book draws every capital and figure of Series A–F as straight lines and circular arcs joined by tangency.
Each drawing has a table of dimensions for 8-, 12- and 18-inch letters. Series B to F are reconstructed here; A is
not yet done.

- **One construction per glyph** (`work/dim<S><page>.py`, e.g. `dimD3.py` for page 3 of Series D, registered in
  `work/series<S>1961.py`). Each is a short program that reads the page's table (inches × 690/8 or × 690/18). It
  builds the letter from half-planes, disks, sectors and polygons, combined with set operations, following the
  book's rules: tangents and major arcs first, shorter connecting arcs fitted afterwards. Unlabelled centres follow
  from tangency to their neighbours, and round tops and bottoms overshoot by about 10.8 units.
- **Vector back end** (`work/vgeom.py`): the coordinates are linear expressions, so a comparison such as `g.X < 120`
  is an exact half-plane. Regions are shapely geometry.
- **Exact outlines** (`work/outline.py`): the lines and circular arcs are recovered from the region's boundary, and
  arcs are drawn as cubic Béziers of at most 45°. The maximum deviation from the construction is 0.012 units.
- **Checks**:
  - `work/scanfit.py` fits each construction onto its page scan (`work/p1961/`). Every glyph is within about 1–5 units
    at the 95th percentile, most within 1–2.
  - `work/tancheck.py` lists joins that are not tangent.
- **Table corrections**: where a printed dimension cannot be tangent to its neighbours, the shorter connecting arc
  is refitted and the extreme point the table fixes is kept. Every refit is listed in `HANDOFF.md`.

### Lower case: the 1966 grid sheets

The 1966 edition gives the Series E(M) lower case only as letters drawn on a grid of 1-inch squares (pp. 33–35, in
`work/lc1966/`), with no dimensions. Each letter is therefore a parametric construction fitted to the scans.

- **Grid** (`work/lcgrid.py`, `work/lcscan.py`): every printed grid line is fitted, so scan points map to
  grid squares. One square is 86.25 units, and the x-height is 6 squares.
- **Curve primitive: the quad oval** (`work/lcgeom.py`). Four axis-aligned quarter-ellipses are joined at their axis
  points, so the curve is tangent-continuous. It is fixed by its left, right, top and bottom extreme points; an
  ordinary ellipse is the special case. Stems, slanted ascender tops, descender ends, terminal cuts and diagonals
  are straight lines.
- **One construction per letter** (`work/lcE.py`), with parameters in grid squares:
  - **n h m r u**: stems and a shoulder of two quarter-ellipses, with a counter tangent to both stems.
  - **a b d g p q**: a bowl of two quad ovals, and a counter ellipse cut by the stem.
  - **t f j**: one hook construction, mirrored.
  - **s**: two outer quad ovals, with counters that rest tangentially on the spine.
  - **k v w x y z**: straight-line polygons (`work/lc1966/poly.json`).
- **Fitting** (`work/lcfitall.py`): a two-sided distance between the model and the scan outline, allowing 1.5 units
  of printing gain per side. The adopted values (`lcE.dims()`) are the fits in `work/lc1966/fit.json`, snapped to
  1/32 square and angles to 0.5°. Most letters fit within 2–4 units at the 95th percentile.
- **Adjustments made after review**, applied on top of the fits:
  - the bowls of b d p q are lightened to the o's thickness (`lcE.BOWL_T`);
  - the w's strokes are thinned (`build1961.W_LIGHT`).
- **Exact outlines** (`work/lcoutline.py`): every boundary edge is classified as straight or as lying on one of the
  recorded quarter-ellipses. Elliptical arcs are drawn as cubic Béziers of at most 45° of parameter.

### Punctuation

The punctuation is derived from each font's own letters, so it follows the series' stroke (`work/punct.py`):

| Glyph | Made from |
|---|---|
| period | the i's dot, on the baseline |
| comma | the y's tail |
| apostrophe, closing single quote | the comma, raised to the cap line |
| opening single quote | the comma turned 180° |
| exclamation point | the i turned upside down |
| double quotes | two single quotes, half a stroke apart |
| hyphen | the H's crossbar |

## Extrapolation to other widths and weights

The 1961 drawings give capitals for Series A–F but no E(M) weight, and the 1966 sheets give the lower case only in
E(M). The missing combinations are derived with one exact geometric operation, the stroke offset (`work/strokeop.py`):

- Every segment moves sideways by d. Lines move along their normal. Circular and elliptical arcs become the
  concentric arc whose radius or semi-axes change by d.
- Corners are re-intersected, with a mitre limit of 4 and bevels beyond it. Segments that vanish are dropped.
- A piecewise-linear vertical map then puts the baseline, x-height, cap line and descenders back on their lines.

Arcs therefore stay exact axis-aligned ellipses, and the outlines stay exact.

**E(M) capitals.** The 1966 sheets ask for "Series E uppercase, but with the stroke width widened to approximately
one-fifth of the letter height". The Series E constructions are offset outward by d = (690/5 − 118.59)/2 = 9.70.
They are then squeezed back to the baseline and cap line by 690/(690 + 2d) (`build1961.py EM`).

**Series B, C, D, E and F lower case.** These come from the E(M) constructions by the rule the 1961 capitals follow
(`work/lccondense.py`, `strokeop.condense2`). The width less the stroke scales by a width factor f. The stroke steps
with the series, as the capitals' table strokes do:

1. Scale x by f.
2. Offset horizontally by (T − V·f)/2, which brings the vertical strokes to the series stroke T (V is the E(M)
   stem, 134.8).
3. Offset vertically by −(V − T)/2, with a stretch that lands the baseline and x-height back on their lines, which
   brings the horizontal strokes to T.

| Series | Capital stroke (1961 table) | Lower-case stroke T | Capitals' width factor | Lower-case width factor |
|---|---|---|---|---|
| B | 86.25 | 83.0 | 0.475 | 0.55 |
| C | 97.03 | 93.8 | 0.65 | 0.82 |
| D | 107.81 | 104.6 | 0.825 | 0.919 |
| E | 118.59 | 115.4 | 1 | 1.035 |
| E(M) | 138 (one-fifth of 690) | 134.8 | — | 1 (as constructed) |
| F | 129.38 | 126.1 | 1.175 | 1.124 |

- **Width factors.** The capitals' factors are measured from the 1961 capitals. The lower-case factors
  (`lccondense.F_LC`) are solved by `work/lcwidth.py` so that each lower-case alphabet (the sum of the a–z
  advances, in cap heights) is as wide as the FHWA 2000 font's. This opens up Series B noticeably.
- **Terminal cuts and slanted stem ends keep their angles.** Each letter is also built without its cuts
  (`lcE.build(c, nocut=True)`), with a short stub beyond each cut. The whole letter is condensed, and each cut is then
  remade at its original angle through the mapped edge point. The diagonals of v w x y z change angle, as the
  capitals' diagonals do from series to series.
- **Dots** of i and j scale with the stroke, with their tops kept on the ascender line.
- **Spacing.** Side bearings come from the matching Highway Even 1966 font (`outputs/HighwayEven1966-<S>.ttf`), and
  for E(M) from `outputs/Highway1966-EM-Constructed.ttf`. These are the earlier traced fonts of this project, whose
  spacing was matched to the FHWA 2000 alphabets.

## Regenerating the fonts

Requirements: Python 3.10 or later (tested with 3.13) and

```
pip install fonttools shapely numpy scipy pillow opencv-python-headless
```

These versions were used: fontTools 4.67, shapely 2.2, numpy 2.5, scipy 1.18, Pillow 12.3, OpenCV 5.0. The
specimens use the DejaVu Sans font for their labels (`/usr/share/fonts/truetype/dejavu/`).

All scripts run from `highway1966/work/` and write to `highway1966/outputs/`:

```
cd highway1966/work
for S in B C D E EM F; do python3 build1961.py $S; done     # about 15–60 s per series; they can run in parallel
for S in B C D E EM F; do python3 strokecheck.py $S; done   # every glyph must be a valid shape (no INVALID lines)
```

`build1961.py <S>` writes `Highway61Revisited-<S>.otf` and `.ttf`. Rebuilding from a fresh clone reproduces the
committed fonts exactly, apart from the timestamp in the `head` table.

The build reads these inputs, all in the repository:

- the construction modules: `dim*.py`, `series*1961.py`, `vgeom.py`, `lcE.py`, `lccondense.py`, `punct.py`, `strokeop.py`
- the lower-case parameters: `work/lc1966/fit.json` and `poly.json`
- the base fonts for side bearings and glyph order: `outputs/HighwayEven1966-{B,C,D,E,F}.ttf` and
  `outputs/Highway1966-EM-Constructed.ttf`

The scans (`work/p1961/`, `work/lc1966/p3?.jpg`) are needed only for refitting and checking.

Specimens:

```
python3 lcseries.py          # outputs/specimen-lc-series.png: all six series
python3 specimen1961.py D    # outputs/specimen-1961-D.png: one series' capitals and figures
```

Refitting and retuning, only when changing the designs:

```
python3 lcfitall.py nhm      # refit the named lower-case letters to the 1966 scans (rewrites lc1966/fit.json)
python3 lcwidth.py B         # solve a series' lower-case width factor; copy the result into lccondense.F_LC
```

### Comparison with the FHWA 2000 fonts

The 2000 FHWA *Standard Highway Signs* alphabets publication embeds its fonts (Series B–F 2000 and Series E Modified
2000) as Type 1 font resources. Given the publication converted to PostScript (for example with `pdftops`):

```
python3 fhwa2otf.py Alphabets.ps   # extracts highway1966/fhwa2000/*.pfa verbatim and converts them to .otf
python3 fhwaspecimen.py            # outputs/comparison-1966-vs-2000.png: two-column specimen at equal cap height
python3 fhwacompare.py             # outputs/compare-fhwa2000.png: per-glyph overlays
```

The extracted fonts carry a URW copyright notice, so they are not committed (`highway1966/fhwa2000/.gitignore`).
`fhwaspecimen.py`, `fhwacompare.py` and `lcwidth.py` need them to be regenerated first.
