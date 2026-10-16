# Highway Gothic 1966 / 1961 font project — handoff

Exported 2026-10-15 from a claude.ai session so work can continue in Claude Code.

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

### Next steps
1. Spacing of the 1961 capitals/figures currently reuses the 1966 E side bearings; check it against the 1961 sheets if
   they give spacing.
2. Do Series A, F page by page with the same structure: read each page's table, write `dim<S><page>.py`,
   overlay against the 1966 tracings (`caps_B/C/D/F.pkl`); Series A has no tracing to compare against.
3. Consider refactoring the per-page code into a shared library of constructions parameterised by series.
