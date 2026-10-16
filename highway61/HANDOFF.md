# Highway Gothic 1966 / 1961 font project — handoff

Exported 2026-10-15 from a claude.ai session so work can continue in Claude Code.

> **2026-10-16 reorganisation:** this directory was renamed from `highway1966/` to `highway61/`, and the
> Highway 61 Revisited fonts (`Highway61Revisited-*.otf/.ttf`) and the specimens `specimen-lc-series.png` and
> `comparison-1966-vs-2000.png` moved to the repository's top level, where the scripts now write them.  References
> below to `outputs/Highway61Revisited-*` and to those two specimens mean the top-level files.
>
> **2026-10-16 cleanup:** the files in `outputs/` left over from earlier stages (superseded fonts, specimens and
> overlays that no current script makes) were deleted, as were the committed copies of the FHWA 2000 fonts
> (`outputs/fhwa2000/`, `work/ps/`); `work/fhwa2otf.py` recreates those from the publication.  `outputs/` now holds
> only the base fonts the build reads and the outputs of current scripts.  Older files named below are in git history.

## Layout of this archive

```
work/                 every script and data file from the working directory (/home/claude/w)
  p1961/              the 1961 "Standard Alphabets for Highway Signs" page scans (from 1961.zip)
  ps/                 Type 1 fonts extracted from the 2000 FHWA PostScript file (.pfa)
  *.pkl               bitmaps and parameters (glyph bitmaps are numpy uint8 arrays, 1 px = 1 font unit)
outputs/              everything delivered so far (fonts, specimens, comparison sheets)
transcripts/          full conversation transcripts from the claude.ai sessions (read these for decisions/history)
```

Not included (regenerable or already yours): `work/hr/` (76 MB of PNG conversions of the hi-res JP2
scans, made by `hires.py` from the `Image063017043007_00NN.jp2` uploads) and the original uploaded scans.
Scripts that reference `/mnt/user-data/uploads/...` or `hr/...` need those paths pointed at your copies.

Font units: cap height 690, x-height 525, ascender 700, descender about -260. 1000 units/em.

## Part 1 — 1966-based fonts (finished, delivered)

Fonts in `outputs/`: `HighwayEven1966-{B,C,D,E,F}.ttf`, `Highway1966-EM-Constructed.ttf`
(plus older traced references `Highway1966-E-Reference.ttf`, `Highway1966-EM.ttf`).
Capitals/figures are traced from the 1966 sheets (`caps_*.pkl`); lowercase is geometrically constructed.

Lowercase construction (dispatch in `assemble.py`, `SRC` dict):
- `geomc2.py` o, e(old) brush ring · `geomc3.py` straight letters + Canvas · `geomc4.py` j (two-ellipse rings)
- `geomc5.py` bowls b d p q (inner ellipse + outward pen profile; `right_to` widens only the counter)
- `geomc6.py` c n m h r u · `geomc7.py` e · `geomc8.py` g · `geomc9.py` t f a · s is traced (`master_mk.pkl`)
- Series parameters in `narrowmk.SER` (width factor f, stroke T): B .542/82, C .78/102, D .89/110, E 1.0/123, E(M) 1.0/140, F 1.10/132
- Pipeline: regenerate letters into `lcE_<series>.pkl` → `buildfam.py <series>` → `fitspace.py <series>`
  (spacing matched to the 2000 FHWA alphabet totals, deltas in `spacing_delta.json`)
- Last changes: a's bowl outline drawn from the original (unwidened) counter so the notch at the stem shows;
  t terminal shift reduced to (TM-T)/3 (in output units) so Series B's tail is longer.

## Part 2 — 1961 dimensioned drawings (in progress)

Goal: capitals and figures built exactly from the 1961 dimension tables (lines + circular arcs), for
Series A (dropped by 1966), B, C, D, E, F.

Page files `p1961/stdalphabets1961-NNN.png`: 000 cover, 008 intro text, then 6 pages per series:
A 010 021 028 036 043 050 · B 058 066 072 079 087 095 · C 102 108 114 121 128 134 ·
D 141 146 152 160 166 173 · E 180 186 193 201 209 216 · F 223 227 235 242 249 255.
Page order within a series: 1 A–F, 2 G–L, 3 M–R, 4 S–X, 5 Y Z 1–4, 6 5–9 0
(check: some series' pages may group differently).

### Series E — complete (all 26 capitals + 10 figures)
- `dim1961.py` page E-5 (Y Z 1 2 3 4) + shared `Grid` (mask raster, 1 unit/px), `raster()` for polygons
- `dimE1.py` A–F · `dimE2.py` G–L · `dimE3.py` M–R (+ `strip`, `cross_solve` helpers) · `dimE4.py` S–X (+ `arc`, `internal_tangent`) · `dimE6.py` 5–9 0 (+ `oval` four-centre oval, `poly_mask`)
- `seriesE1961.py` registry `REG[char]() -> (bitmap, x_left, baseline_row)`; results cached in `seriesE1961.pkl`
- `dimcheck.py` overlays against the 1966 tracing (`caps_E.pkl`) and reports % area difference
- Outputs: `seriesE-1961-constructed.png`, `seriesE-1961-vs-1966-overlay.png`, `pilot-1961-vs-traced.png`

Each page file holds that page's table (8-inch and/or 18-inch row, inches) → font units (×690/8 or ×690/18).

### Conventions learned (apply to the other series)
- The intro (page 008): tangents, major arcs and control points first; connecting arcs of shorter radius
  fitted afterwards — their radii/centres are approximate ("found by experiment"). Stroke uniform unless indicated.
- Curves are circular arcs joined by tangency; most unlabelled centres follow from tangency to neighbours.
- Overshoot (letter J or T or N on a page, about 10.8 units) above cap height / below baseline for round tops/bottoms.
- Heights such as N, F, R, T, H on curved figures are measured from the overshoot lines, not baseline/cap height.
- Large arcs are often centred on the opposite overshoot line (5's E/C arcs, 6/9's L arcs, 2's spine L at the baseline).
- Recurring construction lines (e.g. x = D−O = 275.5 and x = P = 323.4 in Series E); terminals usually cut along the radius.
- Dimension letters are per page; the same letter means different things on different pages.
- Not every stroke is A: L's foot and T's bar use their own dimension.

### Decisions made with the user (Series E)
- 3: bar has a rectangular end at x=P extended as a rectangle to the construction line x=275.5 (no serif-like vertical).
  Its right-side notch: arcs' natural crossing kept (R fits when measured from the overshoot line; J ~13 units off — table error).
- 5: notch where the sloping stem foot meets the bowl's inner E−A arc; lower terminal tip on the left edge (H is ~7 units off).
- N: diagonal's lower edge runs straight to (D−F, 0) and trims the right stem's foot.
- K: arm's lower edge from the stem at height J to the top-right corner (O−N, cap), passing the crotch at cap−M.
- W: Q is the height of the middle counter's apex; middle strokes rise from the bottom flats to it.
- Q: tail rectangle F wide, from (W/2+L, BOT+K) to (W/2+O, BOT−P).
- 6/9 (9 = 6 rotated 180°): tail made a connecting arc tangent to the bowl's N arc at the left edge, centred level
  with the N centres, flat end on the cap line (6) / baseline (9) at M from the side → radius 492 instead of L=420.5.
  This removed the 9's notch; the 6 now sits ~15 units inside the 1966 tracing on its upper left (open question:
  alternative is radius L, tangent, tail 12–16 units wider than B).

### Exact outlines and the 1961 Series E font (done, Claude Code session 2026-10-15)
- `vgeom.py`: vector back end. `VGrid` has the raster `Grid`'s interface, but `g.X`/`g.Y` are linear expressions whose
  comparisons give half-planes, and disks/sectors/polygons are shapely regions combined with & | ~. `install()` patches
  dim1961/dimE1–6 to use it and returns {char: shapely geometry} (K and 9 re-expressed: K's sign-test band as a strip,
  9 as the 6 rotated). The construction code itself is unchanged and still produces the bitmaps (`seriesE1961.pkl`).
- `outline.py`: recovers lines and circular arcs from the polygon vertices (concyclic runs of short chords), draws arcs
  as cubic Béziers of ≤45° each. Max deviation from the exact construction 0.012 units.
- `tancheck.py`: lists joins that aren't tangent (arc pairs by their centre/radius gap, others by angle). Clean now.
- `build1961.py E` → `outputs/Highway61Revisited-E.otf` (font "Highway 61 Revisited E", after the Dylan album) (cubic, fractional coordinates: the exact version) and `.ttf`
  (quadratic, integer). Capitals and figures from the constructions, keeping the 1966 E glyphs' side bearings; lowercase
  and punctuation copied from `HighwayEven1966-E.ttf`. `specimen1961.py E` → `outputs/specimen-1961-E.png`.
- Tangency fixes, found once the outlines were exact (the raster hid them). Following the 1961 note that shorter
  connecting arcs were fitted afterwards, the shorter arc was refitted and the extreme point the table fixes was kept:
  - O/Q (`dimE3._O`): top/bottom radius H 274.9 → 274.3 (the table's H crossed the side arcs with a 2.3° kink).
  - 0 (`dimE6.fig0`): top/bottom radius V 269.5 → 273.1 (missed the side arcs by 2.5 units).
  - 6/9 (`dimE6.bowl6`): bowl side arcs N 194.1 → 192.7 (left), 198.6 (right); the off-centre L arcs (centred at
    x = B−V) can't be tangent to equal side arcs at x = 0 and x = B.
  Inner arcs stay concentric, so the stroke stays A.

### Series D — complete (Claude Code session 2026-10-15)
- `dimD1.py`..`dimD6.py` (pages 141 146 152 160 166 173), written directly against the vector grid (`vgeom.VGrid`);
  `seriesD1961.py` registry; `vgeom.series('D')`; `build1961.py D` → `outputs/Highway61Revisited-D.otf/.ttf`;
  `specimen1961.py D` → `outputs/specimen-1961-D.png`. Tables read for the 8, 12 and 18-inch rows; `dimtab.check`
  finds no row disagreements. Dimensions taken from the 18-inch row.
- `scanfit.py`: fits a construction onto its 1961 drawing (the pages draw outlines) and reports the distance from the
  construction to the drawn lines in font units, with an overlay. Every Series D glyph: p95 ≤ 3.3 units (most ≤ 1.2).
  Series E calibration: p95 ≈ 1.2. Use it for the remaining series; letter boxes on the page are in pixels.
- Constructions that differ from Series E, confirmed by the table and the scan:
  - B: outer bowls H and I both centred on x = D−F−H = D−I, tangent to the cap line/baseline; they cross at the crotch
    (G). Counter ends: arcs K, J level with those centres (which makes them tangent to the middle bar exactly), then
    the outer arcs' offsets H−A, I−A towards the top/bottom bars (tangent to them; they share K's and J's right
    extremes, so the joins are smooth and the stroke there is A). Corrected after review: K, J alone left corners.
  - C: upper arc N and lower arc O differ. G: big arc C = the letter height. O/0 (same dimensions): side arcs 690.
  - S: top arc C and upper bowl D share a centre height (touch on the left); lower bowl C and bottom arc E share a
    centre x (touch at the bottom). Lower tip at x=1.3.
  - 2: spine N centred on the baseline at x=N; bowl P touches the right side, externally tangent to the spine (centre
    height 487.3, drawing ≈489); M, L by tangency, terminal at the left edge O below the overshoot line.
  - 3: notch where bowls S and T cross = 376.7 vs table X 377.3; inner arc V level with T's centre, tangent to T−A and
    to the bar; fillet K above the bar (label read as K).
  - 4: diagonal and crossbar I (97) thick, not A. 1: flag arc I.
  - 6: four-centre oval bowl (B × O) with a tail arc K centred on its centre line; 9 = 6 rotated.
  - 7: stem ring centred on the baseline (computed centre y 0.9). 8: two four-centre ovals, the upper one's side arcs
    S centred above its top arc's centre, the lower one's side arcs R below its bottom arc's centre; notch
    (381.8, 378.2) vs table (382.1, 377.3).
  - W: the steep inner strokes poked 14 units past the feet; the letter is clipped to the outer strokes' edges.
- Decisions taken by analogy with Series E (to confirm): G's jaw at the lower arc's right edge (460.8; B = 463.6 is
  nominal); O/Q/0 top/bottom arcs refitted for tangency 220.4 → 223.4; 6/9 side arcs 198.9 → 195.9; 5's lower
  terminal tip on the left edge (y 23.4 vs the table's I above the overshoot line, 27.6).

### Series C — complete (Claude Code session 2026-10-16)
- `dimC1.py`..`dimC6.py` (pages 102 108 114 121 128 134), `seriesC1961.py`; `build1961.py C` →
  `outputs/Highway61Revisited-C.otf/.ttf`; `specimen1961.py C`. Every glyph fits its drawing with p95 ≤ 1.2 units
  (W: 5.7, its drawing is ~3% wider than its own table's M; the table is followed). Tables agree across rows except
  C-5 N (167–170, 18-inch used) and C-5 D, whose 18-inch entry prints 9 7/8 (9 27/32 used, as the other rows).
- `scanfit.edge_samples` snaps a model's boundary to the drawn lines so the drawing's own arcs can be fitted
  (used to settle the 3 and 8).
- Series C idioms: horizontal terminal cuts at an arc's centre height (C, G, S, 2, 3 lower), with a straight run
  up to a dimensioned height where one is given (C lower: N; S lower: F); two arcs on one centre line touching at a
  side extreme (S, 3, 8, 5); flat-topped bowls with corner arcs whose counters are the bowl offset by A (B, P, R:
  `dimC1._counter_end`).
- Specific constructions: B middle bar centred on the crotch height J; G bar top F above the bottom overshoot line
  (scan); K arm to (L, cap), leg to (B, 0); 2 spine M centred on the baseline at x=M, bowl R by tangency (top 1.6
  above the overshoot line); 3: upper N+L, lower R+Q, T from Q to the notch on L at height D (U then reads 6 units
  wider than drawn), fillets G at the bar; 5: D, F, E arcs, F centred at height G; 6/9: oval with the tail on the
  left side arc M; 7: ring centred on the baseline through (P, 0) (inner edge reaches the bar corner at 377.4 vs B
  377.3); 8: upper circle S + side arcs W, lower R sides + S bottom (refitted 177.8 → 180.7), counters: lower R−A /
  S−A / T, upper a V-stadium whose bottom is A above the lower counter (scan); notch (312.1, 377.6) vs table
  (312.7, 377.3); 0 = O.
- Tangency refits (by the Series E convention): O/Q/0 top arcs 174.9 → 175.3; 6/9 top/bottom arcs L 189.3 → 187.7;
  8 bottom arc 177.8 → 180.7. Unexplained label: 8's D (183.3; the upper centre is TOP−S = 522.9, not TOP−D 517.5).

### Series B — complete (Claude Code session 2026-10-16)
- `dimB1.py`..`dimB6.py` (pages 058 066 072 079 087 095), `seriesB1961.py`, `build1961.py B` →
  `outputs/Highway61Revisited-B.otf/.ttf`. Every glyph fits its drawing with p95 ≤ 2.4 units. Same idioms as Series C.
  Page B-6's 18-inch row is partly illegible: B-6 uses the 12-inch row (agrees with the 8-inch row).
- Differences and decisions: M's diagonals meet at height L (not the baseline); P/R bowl side arc C centred on the
  stem's left edge on the line K below the cap (corners therefore uneven); B counter corners H−A equal the table's I
  (confirms the offset reading); S spine: C and D−A overlapped by 1.6 (an 8.6° kink), C refitted 274.3 → 270.3;
  2: spine centre raised 4.6 above the baseline so the bowl sits on the overshoot line, stays tangent and the spine
  stays vertical at the corner (best fit of four variants); 3: N/Q and R/Q pairs on different centre lines (S, T),
  K joins lower Q to the notch, fillets H; 4: top flat F wider than the stem, diagonal and crossbar V; 1: flag drop I,
  depth H (from the scan; the labels are ambiguous); 7: right edge reaches the bar 5 units past B (clipped at B);
  O/Q/0 top arcs 134.2 → 136.9; 6/9 top arcs 140.2 → 140.6; 8 notch (246.8, 372.9) vs table (237.2, 366.6).
- `dimE6.oval` now samples arcs at 0.25 units (small arcs were partly faceted); all four series rebuilt.

### Series F — complete (Claude Code session 2026-10-16)
- `dimF1.py`..`dimF6.py` (pages 223 227 235 242 249 255), `seriesF1961.py`, `build1961.py F` →
  `outputs/Highway61Revisited-F.otf/.ttf`. Every glyph fits its drawing with p95 ≤ 4.7 units, except M (9.3: its
  drawing is ~1.4% wider than the table; the table is followed). F-2/F-3 use the 12-inch row and F-5 the 8-inch row
  (the scan cuts off the 18-inch row); the rows agree.
- Constructions: B as the Series D B (counter arcs H, G tangent to the bar, outer offsets E−A, F−A); C a ring of
  radius L = half the overshoot height; D, P, R semicircular bowls; O a circle; G jaw at the lower arc's extreme,
  clipped to that arc below its centre (the bar's foot lies lower); S as Series E (top D, bottom C, bowls I, H,
  internal-tangent spine); 2: top O on the axis x=Q, bowl K tangent to it, spine N at (N, 0), straight spine edges
  (internal tangents; the bowl and spine stand 10.8 apart); 3: O arcs top and bottom on x=D−Q, bowls T (I in) and K,
  bar on the notch height, counter ends tangent to O−A and the bar, A in from the bowls — table S 64.7 and F 97 can't
  be tangent to both, refitted to 57.1 and 85.2 (the drawing agrees: measured 56 and 88); 4: diagonal and crossbar U;
  5: bowl top E (centred low on the axis, top D above the overshoot line), right F, bottom C (centred high);
  6/9: oval B × K with K measured from the baseline (then the side arc refits to 194.6 against L 194.1), tail an arc J
  tangent to the bowl's left side and reaching the cap at B−M; 7: inner edge radius C through (B−O, 0) and the bar
  corner (centre 11.7 above the baseline); 8: upper oval S/W with the counter's bottom E (inner, the bowl's E+A),
  lower oval U/R with the counter's top T (inner): notch (543.8, 376.3) vs table (549.8, 377.3); 0: oval with sides J,
  top/bottom refitted 318.6 → 322.1.
- `scanfit.edge_samples` and the 4×4-opening trick (drawn outlines survive, thin dimension lines don't) were used to
  measure arcs where the labels were ambiguous.

### Series E(M) lower case — constructed from the 1966 grid sheets (Claude Code session 2026-10-16)
No dimensioned drawings exist for the lower case; the user supplied hi-res scans of pp. 33-35 (`work/lc1966/p3[345].jpg`)
and overlay drawings showing ellipses (bowls, shoulders, hooks) and straight construction lines (stems, slanted
ascender tops and descender bottoms, terminal cuts, diagonals).
- `lcgrid.py`: every printed grid line fitted individually (common shear, then per-line fit); scan points map to grid
  squares by interpolating between neighbouring lines.  One square = 1 inch at 8-inch letters = 86.25 units; x-height
  6 squares (517.5), ascender 8 (690), descender g/j 3, p/q/y about 2.2-2.4.
- `lcscan.py` / `lcview.py`: sub-pixel letter outlines in grid squares / font units (`lc1966/contours.pkl` cache).
  Printing gain measured from flat baselines and x-heights: 1.5 units per side (`lcmodel.SPREAD`).
- `lcgeom.py`: the curve primitive is the **quad oval** — four quarter-ellipses (axis-aligned) joined at their axis
  points, so tangent-continuous; an ellipse is the special case.  This matches the user's paired ellipses (e.g. e's
  left half from one ellipse, right half from another).  Every quarter used is recorded in `lcgeom.REG`.
- `lcE.py`: one construction per letter (o n a e g l i h m u r b d p q c t f j s), parameters in squares in
  `lc1966/fit.json`, adopted values = fit snapped to 1/32 square, angles to 0.5 degree (`lcE.dims()`).
  k v w x y z are straight-line polygons (`lcpoly.py` -> `lc1966/poly.json`: each edge fitted as a line, vertices
  = intersections, snapped to 1/32 square, baseline/x-height/ascender exact).
- `lcfitall.py <letters>`: fits a construction to the scan (two-sided distance, soft-L1, plus a per-letter x offset
  `ox` for grid registration); `evaluate()` reports median/p95 of the snapped values and writes overlays.
- `lcoutline.py`: exact outlines: each polygon edge is classified as straight or as lying on one of the recorded
  quarter-ellipses; runs become elliptical arcs drawn as cubic Beziers (<= 45 degrees of parameter each).
- `build1961.py E` now replaces a-z with these constructions (side bearings kept from Highway Even 1966 E);
  `lcspecimen.py` -> `outputs/specimen-lc-E.png` (text + every letter over its scan).
- Constructions in brief (all joins tangent unless a corner is drawn in the scan; `lcnotch.py` flags any corner
  that touches a short segment — the only ones left are the bowl/stem notches at the baseline of b and d):
  n/h/m/r/u: stems about 1 9/16 square, shoulder = two quarter-ellipses sharing the top point (right one tangent to the
  right stem's outer edge); counter = quarter-ellipses tangent to both stems' inner edges, sharing a top point
  (`counter_arch`); m's second arch reaches over the middle stem; u = n turned 180 degrees.
  a/d/g/q/b/p: bowl = two quad ovals sharing the left (b/p: right) extreme and top/bottom points, counter ellipse cut
  by the stem edge.  a: arch quad oval minus an underside quad oval tangent to the stem's inner edge; gap = underside
  outside the bowl.  e: outer quad oval cut by the bar's vertical end; counters = two ellipses; terminal 3:4 slope.
  c: outer quad oval, counter ellipse, opening = wedge between two terminal lines.  g: tail = outer/inner ellipses.
  t/f/j: one hook construction (stem, two quarter-ellipses sharing the bottom point, below the right one's vertical
  tangent the edge runs straight on; inner ellipse tangent to the stem; terminal cut beyond the inner ellipse's
  centre), mirrored for j and f.  s: two outer quad ovals (each gives one arc, one side, one spine edge); each counter
  is a quad oval shifted vertically until it touches the spine edge from outside (`rest`), so the counter's round end
  runs tangentially into the spine edge; past the contact the counter is bounded by its own far edge and opens past
  its terminal line; the counter and outer ovals all run out past the terminal lines so every cut ends in a corner.
  Ascender tops slope 0.4 (2:5), p/q descender bottoms 0.354.
- Fit quality (snapped, p95 of scan<->model distance, units): most letters 2-4; f, s about 5, t 5.5 (the tail's
  lower right is straighter than an ellipse), i 4.6 (dot / grid offset).  Not yet harmonised across letters
  (stem widths 1.53-1.57 squares, overshoots 0.06-0.11) — a candidate for a cleanup pass with shared values.
- Note the x-height mismatch with the 1966 tracings (525) — the constructions use 6 squares = 517.5.

### Series E(M) capitals and Series E lower case by stroke change (Claude Code session 2026-10-16)
The 1966 sheets recommend capitals and numerals for the E(M) lower case be "Series E uppercase, but with the stroke
width widened to approximately one-fifth of the letter height".  `strokeop.py` changes stroke weight on the exact
outlines: every segment moves sideways by d (lines along their normal, circular/elliptical arcs to the concentric
arc with radius / semi-axes changed by d), corners are re-intersected, segments that vanish or reverse are dropped,
then a piecewise-linear vertical map restores the metrics (arcs stay exact axis-aligned ellipses).
- E(M) capitals and figures (`build1961.py EM` -> `Highway61Revisited-EM.otf/.ttf`): Series E constructions offset
  by d = (690/5 - A)/2 = 9.70 (A = 3 3/32 in at 18 in = 118.59) and squeezed about mid-height by 690/(690+2d) back to
  baseline and cap line (vertical strokes 138, horizontal 134; letters about 19 units wider); an acute corner pushed
  through the baseline (M's middle point) is clipped back to it.  Lower case: the E(M) constructions; spacing (side
  bearings) from `Highway1966-EM-Constructed.ttf`.
- Checked against the 1983 E(M) tracings (`caps_EM.pkl`): `emverify.py` -> `outputs/verify-EM-caps.png`; median
  distance 1-10 units per glyph, p95 mostly 5-19 (mean 13.9).  Straight letters match closely; the 1983 drawings
  redraw some details (L's bar, J narrower, G and R bars, S, 6, 9), which is where the larger differences are.
  Choosing between offset-and-squeeze and offset-and-clip: squeeze fitted better (mean p95 12.8 against 15.2).
- Series E lower case (`build1961.py E`): the opposite change, offset by -9.70, then a vertical map in three bands
  (below the baseline shifted down by d, baseline to x-height stretched by X/(X-2d), above shifted up by d) so the
  baseline, x-height, ascender and descenders stay on their lines.  No pre-digital Series E lower case exists to
  check against.  Curves crossing the band lines get a tangent change of at most about 1 degree.
- `strokecheck.py E|EM`: every glyph's flattened outline must be a valid shape (all pass).

### Condensed (B, C, D) and expanded (F) lower case (Claude Code session 2026-10-16)
`lccondense.py` makes each series' lower case from the E(M) constructions, by the rule the 1961 capitals follow:
the width less the stroke scales by f = 0.475 / 0.65 / 0.825 / 1.175 (B C D F; measured from the 1961 capitals), and
the stroke steps with the series as the capitals' table strokes do (B 83.0, C 93.8, D 104.6, E 115.4, F 126.1 units;
E(M) 135).  Both vertical and horizontal strokes come out at the series stroke, and terminal cuts and slanted stem
ends keep their angles.
- `strokeop.condense2`: scale x by f; offset the outline horizontally by gx = (T - V f)/2 (vertical strokes to T);
  stretch y so that a vertical offset by -ey = -(V - T)/2 (horizontal strokes to T) lands the baseline and x-height
  back on their lines; dots stay round.  (An earlier "thin core" method — offset in, scale, offset out — collapsed
  every feature thinner than twice the core offset, e.g. the t's tail; it is still in strokeop.condense, unused.)
- Terminals: each letter is also built without its cuts (`lcE.build(c, nocut=True)`, the `NOCUT` switch); the cut
  lines are the cut letter's outline lines inside the uncut letter (`cuts`).  A short stub of the uncut letter is
  added beyond each cut (kept clear of other parts of the letter), the whole is condensed, and each cut is made
  again at its original angle through the mapped edge point at its middle; what is removed beyond it is the condensed
  image of the stub (condensed separately with a band of the letter on the near side), and any small ears left at
  the cut's two ends.  k (a straight-line letter) gets the same treatment for its slanted stem top; v w x y z are
  condensed directly (their diagonals change angle, as the capitals' do).
- Outlines are recovered exactly (lines and elliptical arcs, lcoutline) from the result; `strokecheck.py <S>`
  validates every glyph.  `build1961.py B|C|D|F` now uses these lower cases; specimen `outputs/specimen-lc-series.png`.
- Spacing (side bearings) is still that of the Highway Even 1966 fonts' lower case for each series.

### Dots and punctuation (Claude Code session 2026-10-16)
- i and j dots in B C D F scale with the stroke: radius = E(M) dot radius x series stroke / E(M) stroke (B 44.8,
  C 50.7, D 56.5, F 68.1), the dot's top kept on the ascender line (`strokeop.condense2(dot_scale=...)`).
- `punct.py`, used by `build1961.py` for every series: period = the font's i dot on the baseline; comma = the y's
  tail below the baseline, raised so its top is level with the period's top; quoteright and quotesingle = the comma
  raised so its top is on the cap line; quoteleft = the comma turned 180 degrees with its top on the cap line;
  exclam = the i turned upside down, its dot on the baseline; quotedbl,
  quotedblright, quotedblleft = two of those, half a cap stroke apart; hyphen = the H's crossbar between the stems.
  All are lines and circles (exact).  Side bearings from the base fonts' glyphs.

### Lighter b d p q bowls (Claude Code session 2026-10-16)
- The user found the bowls of b d p q heavier than their neighbours: the fitted counters left them 1.59-1.66 squares
  thick against the o's 1.56 (sides) and 1.50 (top/bottom).  `lcE.dims()` now enlarges each counter ellipse (centre x
  kept, snapped to 1/32) so the bowl is `lcE.BOWL_T` = (1.53125, 1.5) squares thick at its free side, top and bottom
  (`lcE.lighten`); fit.json is unchanged.  Every series inherits it through the E(M) constructions.

### FHWA 2000 comparison (Claude Code session 2026-10-16)
- `work/fhwa2otf.py <Alphabets.pdftops.ps>` writes the Type 1 fonts embedded in the 2000 Standard Highway Signs
  Alphabets PostScript (`%%BeginResource: font ...+SeriesB2000-Normal` etc.) verbatim to `fhwa2000/<name>.pfa` and
  converts them to `.otf` (outlines and advances unchanged, cmap from glyph names, 250-unit space added).  Series B, C,
  D, E, ModE and F 2000 have A-Z a-z 0-9 and ! " # $ & ( ) * + , - . / : = ? @ cent (no apostrophe or curly quotes,
  no space); HighwayD98 and HighwayME98 are fragments.  The fonts carry a URW 1998 copyright notice, so they are
  git-ignored; regenerate them from the PostScript.  Their caps are 285 units of a 1000 em, x-height 208.
- `work/fhwaspecimen.py` -> `outputs/comparison-1966-vs-2000.png`: two-column specimen (B C D | E E(M) F), sign copy,
  complete lower case, capitals, figures and punctuation in each font at equal cap height with natural spacing (the
  FHWA fonts get our space of the same series, in cap heights).
- `work/fhwacompare.py` -> `outputs/compare-fhwa2000.png`: each series (E(M) vs ModE) at equal cap height, two lines
  in each font and a per-glyph overlay (ours red, FHWA blue, both black).

### Lower-case widths matched to FHWA 2000 (Claude Code session 2026-10-16)
- At the user's request the lower case of each series is as wide in total (sum of the a-z advances, in cap heights)
  as the FHWA 2000 font's: `lccondense.F_LC` = B 0.55, C 0.82, D 0.919, E 1.035, F 1.124 (were 0.475, 0.65, 0.825,
  1.0, 1.175, the capitals' 1961 factors, which stay in `lccondense.F`).  Solved by `lcwidth.py <S>` (secant on the
  build's own advances).  Before: B -6.7%, C -11.7%, D -5.7%, E -1.8%, F +2.5% against FHWA; now all within 0.1%.
  E(M) was already within 0.3% and is unchanged.
- The Series E lower case now goes through `lccondense.series('E')` like the others (f = 1 reproduced the old
  offset-and-bands narrowing to within 0.3% in width), so it can be widened with its cuts kept at their angles.
- Ear removal after a re-cut (lccondense) now takes the two ends of the straight run of outline along the cut, not
  the outermost vertices anywhere on the cut line: at B 0.55 the inside of the t's foot crossed the t's terminal
  line and had a 15-unit notch taken out of it.  `series(S, f, only)` builds a subset of letters for tests.

### Lighter w (Claude Code session 2026-10-16)
- `build1961.W_LIGHT` = 0.15: in every series the w is offset inward by 0.075 x the series stroke and its baseline
  and x-height restored (vmap bands), so its strokes are about 0.82 of the v's across (FHWA Series B 2000 w: 0.80;
  their other series leave the w at v weight).  The offset also deepens the crotches.

### Series A — complete, with a lower case (Claude Code session 2026-10-16)
- `dimA1.py`..`dimA6.py` (pages 010 021 028 036 043 050), `seriesA1961.py`, `build1961.py A` ->
  `Highway61Revisited-A.otf/.ttf`.  Rows agree on every page.  Mostly the Series B constructions with this series'
  letters; every glyph fits its drawing with p95 <= 3.3 units.  Differences and decisions:
  - A-4 W: the middle strokes run full height as Series B's do (a box too narrow first made it look different).
  - A-5 2: bowl L touching the top overshoot line and the right side; spine M (centre A below the top overshoot line,
    internally tangent to L) running into N-A (N centred O below the baseline, on the right, reaching the bottom-left
    corner); N's centre moved 1 unit for exact tangency.  3: the notch is at (D-F, W above the bottom overshoot line),
    joined to T and U by arcs S whose concentric inner arcs S-A equal the table's H; P refitted 118.6 -> 122.8 so the
    upper terminal reaches the left side.  4: the crossbar runs H past the D width.
  - A-6 5: heights B and I measured from the bottom overshoot line (p95 3.3, from the baseline 8.2); 6/9: tail arc K
    centred M above the baseline, tangent to the bowl's left side, its tip on the top overshoot line; 7: the stem's
    centre solved from (Q, 0) on its left edge and (B, cap-A) on its right (it lies 488 below the baseline, not on it).
  - O/Q/0 end arcs refitted for tangency 118.6 -> 120.7; S spine C 345.0 -> 343.5.
- Capital width factor 0.450 (median core ratio to Series E, the method that gives B 0.475); lower case: stroke 61.5,
  `lccondense.F_LC['A']` 0.521 (B's ratio F_LC/F); spacing from Series B scaled (`build1961.SB_SCALE` 0.766,
  `SPACE_SCALE` 0.814, from a linear fit of the capitals' side bearings against the stroke over B..F).
- Since the lower case was widened to the FHWA totals, the horizontal offset gx = (T - V f)/2 is negative in A, C, D,
  E and F (the outline is eroded sideways), which left nicks and hooks of a few units at re-cuts and junctions.
  `lccondense._tidy` (an opening by 4 units, mitred) now cleans them when gx < 0; a closing as well was tried and
  rejected (it fills in beside acute terminal tips).  Stray pieces are dropped unless they are the i/j dot.  The stub
  beyond a cut must also keep 3 units clear of the rest of the letter (the e's stub had reached its bar).

### e: rounded into the bar in E(M) (Claude Code session 2026-10-16)
- The user found the E(M) e's right side still flat: the outer oval's right extreme was 30 units above the bar's lower
  edge, so the side ran straight down below it.  `lcfitall.TIES['e']` now ties `Oyr` to `yb` (the extreme on the bar's
  lower edge) and the e was refitted (p95 4.44 against 4.50 before).  Every series' e changes by 5-8 units.

### s as ellipses (Claude Code session 2026-10-16)
- Following the user's reading of the 1966 drawing, `lcE.letter_s_ell` replaces the two-quad-oval s
  (`letter_s_quads`, kept).  Upper family, all through the top point: U1 (big: the top stroke's upper edge to the
  upper terminal line), U2 (upper left side; its bottom arc the spine's lower edge), UM (its bottom arc the bottom
  stroke's upper edge).  Lower family through the bottom point: L1, L2, LM likewise.  The spine's middle is straight
  (the user's suggestion): its edges leave L2's top arc and U2's bottom arc along their tangents at x = s1, s2.  The
  counters are the bands inside LM above the spine and inside UM below it, each ended by a small ellipse (C1, C2:
  centre x and half-width fitted) that touches both band edges exactly (`_touching`).  Terminal lines cut U1/L1 and,
  beyond the terminal corners, U2/L2 outside LM/UM.  Fit p95 4.12 (snapped) against the quad-oval s's 4.77.
- The condensed and expanded s are built by mapping the construction, not by offsetting its outline
  (`lccondense._s_direct`): ink ellipses a*f+gx, b*k-ey; counter ellipses a*f-gx, b*k+ey; centres (x f, y k-ey);
  C1/C2 and s1/s2 scaled; each cut keeps its angle through the mapped middle of the original cut.  Offsetting the
  outline parted its non-axis joins (strokeop moves an arc by changing its semi-axes, exact only at axis points), so
  the construction re-solves its tangencies instead.

### Spacing by the standard tables (2026-10-16)
- `spacingcheck.py` transcribes the 1966/1983 spacing tables (Table I letter codes, Table II figure codes, Table VI
  factors 1.5/1.2/0.8/0.4 of the series stroke, E(M) at Series E's spaces; the 1983 lower-case table in mm at 375 mm
  loop height) and measures the fonts against them, including GPOS kerning.  Before: capitals off by 10-27 units on
  average (max ~96), E(M) lower case by 30 (max 105); side bearings alone could reach only ~7-11 mean.
- `spacing.py` solves per series: right side bearing per preceding character, one left side bearing per following
  class, least squares with a weak (0.05) pull to the old bearings; remainders become class kerning
  (`pos <glyph> @class <value>`, ~2000 pairs, within about +-48 units), added by feaLib in `build1961.setup`.
  Lower case for series other than E(M): the table scaled by `lccondense.stroke(S)/V_EM` (the user's choice).
- Result: all 776 capital/figure and 1352 lower-case table pairs within 1 unit in every series; outlines unchanged.
- Lower case re-tightened (user: the fonts must be drop-in replacements for the FHWA ones): `spacing.LC_TIGHT`
  multiplies the stroke-scaled table so the a-z advance sum equals FHWA 2000's (`lcspace.py`, imports build1961
  without writing): B 0.608, C 0.635, D 0.736, E 0.508, F 0.814; A uses B's; E(M) is the table unscaled (canon;
  its alphabet 23.25 cap heights against FHWA ModE's 23.95).  Totals now within 0.03 cap heights, kerned or not.
- Series E lower-case width factor `F_LC['E']` 1.035 -> 1.0 (user): its core now equals E(M)'s, as the capitals'
  does; the FHWA alphabet width is made up by spacing (LC_TIGHT E 0.447 -> 0.508).
  Uncovered pairs (lc before capital, letter/figure, punctuation, space) rely on side bearings only.

### Series A 8 waist, condensed s spine (2026-10-16)
- Series A 8 (`dimA6.fig8`): with the lower bowl's top arc W the waist between the counters came out 36 units, not
  the dimensioned A (64.7).  That arc is refitted (W 107.8 -> 117.8, still tangent to the side arcs S) so the lower
  counter's top is A below the upper's bottom; the notch moves from 375.1 to 362.9 (table U: 366.6).
- Condensed s spine (`lccondense.S_SPINE`): its tangent points s1/s2 moved towards the lobes' extremes by 9 (A) and
  12 (B) units; the narrowest width (twice the erosion that splits the s) goes from 0.77/0.76 to 0.89/0.87 of the
  stroke (E(M) 0.81, C-F 0.82-0.83).

### Next steps
0. Lower case: harmonise shared dimensions (stem, overshoot, bowl family a/b/d/g/p/q), improve t, f, s; the user
   reviewed the first full set; their notches/joins list (a f h n r s t u) is fixed (2026-10-16).
1. (done) Spacing from the 1966/1983 tables.  Possible follow-up: punctuation and figure/letter pairs, which the
   tables do not cover.
2. (done) Series A.
3. Consider refactoring the per-page code into a shared library of constructions parameterised by series.
