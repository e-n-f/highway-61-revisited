# Highway 61 Revisited

The **Highway 61 Revisited** family of fonts is a faithful revival of the
[1961 Standard Alphabets for Highway Signs](https://www.trafficsign.us/oldmutcd/oldstds/stdalphabets1961.pdf), including the long-forgotten Series A,
and of the [1966 Standard Lower Case Alphabets for Highway Signs](https://archive.org/details/standard-alphabets-for-highway-signs-1966/page/n29/mode/2up),
with the addition of a pleasing, consistent generalization of the
1966 lower case to the widths and weights of the rest of the series. The "/" character comes from the [1958 Manual for Signing and Pavement Marking of the National System of Interstate and Defense Highways](https://mutcd.kittelson.com/wp-content/uploads/2022/11/1958InterstateSignsMarkings-1.pdf).

<img width="6987" height="1750" alt="Image" src="https://github.com/user-attachments/assets/b08d4dba-1283-4cf2-a15d-7a7e753b94ff" />

It is a protest against the ugly and inconsistent lower case and punctuation
of the narrow widths of the [2000 Standard Alphabets for Traffic Control Devices](https://mutcd.fhwa.dot.gov/shse/alphabets.pdf),
which since 2009 have been required on American street name signs.

<img width="6789" height="1500" alt="Image" src="https://github.com/user-attachments/assets/d14a2244-b723-4200-8226-2b9df44d0b3e" />

The 1961 and 1966 letterforms are spaced and kerned as specified in
the original documents; the other lower cases have tighter spacing
like those from 2000.

----------

| Font | Files |
|---|---|
| Highway 61 Revisited A, B, C, D, E, F | `Highway61Revisited-{A,B,C,D,E,F}.otf` / `.ttf` |
| Highway 61 Revisited E(M) | `Highway61Revisited-EM.otf` / `.ttf` |

The `.otf` files are the exact versions: cubic Béziers with fractional coordinates. The `.ttf` files are
quadratic with integer coordinates. Each font has A–Z, a–z, 0–9, `. , ! ' " ‘ ’ “ ” -` and the space.

Font units: 1000 per em. Cap height 690 (an 8-inch letter, so 1 inch = 86.25 units), x-height 517.5, ascender 690,
descenders about −200 to −260.

Specimens: `specimen-lc-series.png` (all six series) and `comparison-1966-vs-2000.png` (against the FHWA 2000
fonts).  The fonts and these two specimens are at the top level of the repository; everything used to make them
is under `highway61/`. The detailed working history, with
every construction decision, is in `highway61/HANDOFF.md`.

## How the letterforms are represented

Paths in this section are relative to `highway61/`.

### Capitals and figures: the 1961 dimensioned drawings

The 1961 book draws every capital and figure of Series A–F as straight lines and circular arcs joined by tangency.
Each drawing has a table of dimensions for 8-, 12- and 18-inch letters. All six series, A to F, are reconstructed here.

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
| A | 64.69 | 61.5 | 0.450 | 0.521 |
| B | 86.25 | 83.0 | 0.475 | 0.55 |
| C | 97.03 | 93.8 | 0.65 | 0.82 |
| D | 107.81 | 104.6 | 0.825 | 0.919 |
| E | 118.59 | 115.4 | 1 | 1 |
| E(M) | 138 (one-fifth of 690) | 134.8 | — | 1 (as constructed) |
| F | 129.38 | 126.1 | 1.175 | 1.124 |

- **Series A** was dropped in 1966 and never had a lower case or an FHWA digital font. Its lower case follows the
  same rule: stroke 61.5 by the series step, and a width factor keeping Series B's ratio of lower-case to capital factor.
  Its punctuation and space are spaced as Series B's, with every side bearing scaled by 0.766 and the space by 0.814:
  across B to F the capitals' mean side bearing is linear in the stroke, and that is its value at Series A's stroke.
  Letters and figures follow the spacing tables like the other series (below).
- **Width factors.** The capitals' factors are measured from the 1961 capitals. The lower-case factors
  (`lccondense.F_LC`) are solved by `work/lcwidth.py` so that each lower-case alphabet (the sum of the a–z
  advances, in cap heights) is as wide as the FHWA 2000 font's. This opens up Series B noticeably. Series E is
  the exception: its factor is 1, so that its lower case has E(M)'s forms exactly (as its capitals have), and its
  letterspacing makes up the alphabet width instead.
- **Terminal cuts and slanted stem ends keep their angles.** Each letter is also built without its cuts
  (`lcE.build(c, nocut=True)`), with a short stub beyond each cut. The whole letter is condensed, and each cut is then
  remade at its original angle through the mapped edge point. The diagonals of v w x y z change angle, as the
  capitals' diagonals do from series to series.
- **Dots** of i and j scale with the stroke, with their tops kept on the ascender line.
- **The e's bar and the a's bowl top** get heavier with condensation (`lcE.condensed_ae`). In E(M) they are
  lighter than the top and bottom strokes, as drawn in 1966. Two early-1960s Series C signs on the Oregon Expressway
  in Palo Alto (Bayshore Fwy, Business District) have them as heavy as the other strokes. The extra weight goes into
  the e's lower opening and into the a's gap under the arch, with a slight further move down for the bar and up for
  the bowl top; the e's eye and the a's bowl counter keep their size. The change is in proportion to the condensation
  (1 − the lower-case width factor), full at Series C's, and never heavier than the stroke. Series E, E(M) and F are
  unchanged.
- **Spacing** follows the spacing tables of the 1966 and 1983 *Standard Alphabets* (`work/spacingcheck.py` has
  them). For capitals and figures, Tables I and II give a code (1–4) from the preceding character and the group of
  the following one (B D E F H I K L M N P R U | C G O Q S X Z | A J T V W Y, and 1 5 | 2 3 6 8 9 0 | 4 7); Table VI
  gives the space for each code, between the two characters' extreme edges, as the series stroke times 1.5, 1.2,
  0.8 or 0.4. E(M) uses Series E's spaces, as both editions do, and Series A, which has no table, uses the same
  factors of its own stroke. For the lower case, the 1983 table gives the space from the preceding character
  (capital or lower case, in 22 classes) and the following lower-case class (a c d e g o q | b h i k l m n p r u |
  f w | j | s t | v y | x | z). It is for E(M), whose lower case follows it as it stands. The other series get it
  scaled by their lower-case stroke over E(M)'s and then tightened uniformly (`spacing.LC_TIGHT`, 0.51–0.81) so
  that each lower-case alphabet (the sum of the a–z advances) is as wide as the FHWA 2000 font's, as the
  lower-case widths are; set in text, a line of lower case runs the same length as in the official font, so these
  can stand in for it. Series A, with no FHWA font, takes Series B's factor.

  The tables are not additive, so they are met with side bearings and class kerning (`work/spacing.py`): each
  character gets its own right side bearing, every character of a following class the same left side bearing,
  chosen by least squares (pulled weakly towards the earlier bearings), and what remains for each (character,
  following class) pair is a GPOS `kern` value. All 2128 table pairs of each font come out within a unit of the
  table, as scaled (`python3 spacingcheck.py`). Pairs the tables do not cover (a lower-case letter before a capital, letters
  next to figures, punctuation, the space) are spaced by the side bearings alone. Punctuation keeps its side
  bearings from the matching Highway Even 1966 font (`highway61/outputs/HighwayEven1966-<S>.ttf`, for E(M)
  `highway61/outputs/Highway1966-EM-Constructed.ttf`), the earlier traced fonts of this project.

## Regenerating the fonts

Requirements: Python 3.10 or later (tested with 3.13) and

```
pip install fonttools shapely numpy scipy pillow opencv-python-headless
```

These versions were used: fontTools 4.67, shapely 2.2, numpy 2.5, scipy 1.18, Pillow 12.3, OpenCV 5.0. The
specimens use the DejaVu Sans font for their labels (`/usr/share/fonts/truetype/dejavu/`).

All scripts run from `highway61/work/`.  The fonts and the two headline specimens are written to the top level
of the repository, and other specimens and checks to `highway61/outputs/`:

```
cd highway61/work
for S in A B C D E EM F; do python3 build1961.py $S; done     # about 15–60 s per series; they can run in parallel
for S in A B C D E EM F; do python3 strokecheck.py $S; done   # every glyph must be a valid shape (no INVALID lines)
python3 spacingcheck.py                                        # every table pair within a unit of the standard
```

`build1961.py <S>` writes `Highway61Revisited-<S>.otf` and `.ttf` at the top level. Rebuilding from a fresh clone reproduces the
committed fonts exactly, apart from the timestamp in the `head` table.

The build reads these inputs, all in the repository:

- the construction modules: `dim*.py`, `series*1961.py`, `vgeom.py`, `lcE.py`, `lccondense.py`, `punct.py`, `strokeop.py`,
  and `spacing.py` with the tables in `spacingcheck.py`
- the lower-case parameters: `work/lc1966/fit.json` and `poly.json`
- the base fonts for glyph order and the punctuation's side bearings: `highway61/outputs/HighwayEven1966-{B,C,D,E,F}.ttf` and
  `highway61/outputs/Highway1966-EM-Constructed.ttf`

The scans (`work/p1961/`, `work/lc1966/p3?.jpg`) are needed only for refitting and checking.

Everything in `highway61/outputs/` is either one of these base fonts or the output of a current script:
`specimen1961.py` (`specimen-1961-<S>.png`), `lcspecimen.py` (`specimen-lc-E.png`), `emverify.py`
(`verify-EM-caps.png`) and `fhwacompare.py` (`compare-fhwa2000.png`). The base fonts come from the project's
earlier traced fonts and are not regenerated by the current process.

Specimens:

```
python3 lcseries.py          # specimen-lc-series.png (top level): all seven series
python3 specimen1961.py D    # highway61/outputs/specimen-1961-D.png: one series' capitals and figures
```

Refitting and retuning, only when changing the designs:

```
python3 lcfitall.py nhm      # refit the named lower-case letters to the 1966 scans (rewrites lc1966/fit.json)
python3 lcwidth.py B         # solve a series' lower-case width factor; copy the result into lccondense.F_LC
python3 lcspace.py B         # then its lower-case spacing factor; copy the result into spacing.LC_TIGHT
```

### Comparison with the FHWA 2000 fonts

The 2000 FHWA *Standard Highway Signs* alphabets publication embeds its fonts (Series B–F 2000 and Series E Modified
2000) as Type 1 font resources. Given the publication converted to PostScript (for example with `pdftops`):

```
python3 fhwa2otf.py Alphabets.ps   # extracts highway61/fhwa2000/*.pfa verbatim and converts them to .otf
python3 fhwaspecimen.py            # comparison-1966-vs-2000.png (top level): two-column specimen at equal cap height
python3 fhwacompare.py             # highway61/outputs/compare-fhwa2000.png: per-glyph overlays
python3 linespecimen.py            # specimen-lines-h61.png, specimen-lines-fhwa.png (top level): a line per font
```

The repository does not include these fonts (`highway61/fhwa2000/` is git-ignored). They are part of the official
publication, and `fhwa2otf.py` recreates them from it. Run it before `fhwaspecimen.py`, `fhwacompare.py`, `linespecimen.py`,
`lcwidth.py` or `lcspace.py`, which need them.
