# Crown-Delineation Results: Non-Technical Report

Prepared: 2026-06-16

This report summarizes the crown-delineation experiments in the `results/`
folder — specifically the work in `crown-segmentation-results.md`. Crown
delineation is the step *after* tree-top detection: once a tree top is found,
the method draws the outline of that tree's crown and measures how wide it is.
The measured width is then compared against crown diameters surveyed in the
field. Deep-learning / machine-learning methods are out of scope for this
version; only the classical (non-ML) methods are covered here.

The short version is:

- After tree tops are found, several classical methods can estimate crown width
  to within roughly 2.4–2.7 m on average. Field crowns here average about 5 m
  wide, so the error is real but usable for many forest-mapping purposes.
- The best classical methods are very close to each other. *Having a sensible
  rule for where a crown stops growing matters more than which method you pick.*
- Better tree tops produce better crowns. Improving the detection step feeds
  straight through into more accurate crown widths.
- Crown width survives sparse data. Unlike understory tree counting, crown
  width stays accurate even when the LiDAR is thinned to 1 point/m^2.
- The "widest-axis" definition of crown diameter stays hard for every method,
  and the main failure mode is crowns being drawn too large.

## 1. What Was Measured

The benchmark uses NEON 2021 airborne LiDAR plus NEON field-surveyed trees at
three California forest sites, the same data as the tree-top detection work:

| Site | Plain-English forest type | Why it matters |
| --- | --- | --- |
| SJER | Open oak / foothill-pine savanna | Open crowns, easiest to separate |
| SOAP | Mixed conifer forest | Main anchor case; multi-layered canopy |
| TEAK | Dense red-fir conifer forest | Closed, packed crowns; hardest |

In total, 225 field trees across 40 plots were matched and scored (SJER 22,
SOAP 87, TEAK 116). Field crowns averaged about 5 m wide (the "equivalent
width" measure) to about 6 m (the "widest-axis" measure).

The key design choice: **every method starts from the same tree tops.** Tree
tops are detected once per plot, then handed to each crown method as the same
set of seeds. That way, any difference between methods is a difference in *how
they grow the crown*, not a difference in which trees they started from.

## 2. The Two Diameter Definitions (Important)

A crown is not a perfect circle, so "crown diameter" can mean two different
things, and they are **not interchangeable**:

- **Equivalent-circle width** — the diameter of a circle with the same area as
  the crown outline. This is the natural match to NEON's "ninety crown
  diameter" field column.
- **Widest-axis width** — the single longest straight line across the crown
  outline. This is the natural match to NEON's "maximum crown diameter" field
  column, but it is always biased a bit large for any irregular shape.

Across all methods, errors on the widest-axis measure run about 1.4–1.6× larger
than errors on the equivalent-circle measure for the *same* crowns. The
practical rule: compare like with like. Use equivalent-circle width against the
field's ninety-crown-diameter, and widest-axis against the field's
maximum-crown-diameter — never cross them.

## 3. How the Methods Differ

Each method takes the same canopy height model (a top-down height map of the
canopy surface) and the same tree tops, then decides which pixels belong to
which crown. The plain-English differences:

| Method | What it does, in plain English |
| --- | --- |
| lasR region growing | Grows each crown outward from its top, with a built-in rule to stop before it reaches its neighbor. Conservative — keeps crowns compact. |
| dalponte2016 | A similar grow-from-the-top method with a height-based stop rule. |
| silva2016 | Splits the canopy between tops, Voronoi-style (each pixel goes to its nearest top). |
| watershed (marker-free) | Treats the canopy like a landscape and floods the valleys; finds its own basins, which tends to over-grow crowns. |
| watershed (seeded) | The same flooding, but anchored so there is exactly one crown per tree top. |
| random walker | Spreads each crown label outward like a diffusion; with no stop rule it tiles the whole canopy and over-grows. |
| random walker + stop rule | The same method, but each crown is cut off below a fraction of its own top's height so it stops creeping into the gaps. |

There are also three **3-D point-cloud methods** (Li 2012, ptrees, AMS3D) that
work directly in the 3-D points instead of on the height map. These are
classical algorithms, not machine learning. They behave differently from the
height-map methods and are discussed in Section 6.

## 4. Main Result: the Classical Methods Are Close, and a Stop Rule Wins

Pooled over all 225 trees, on the **equivalent-circle width** (lower is
better):

| Method | Average error (RMSE) | Tendency |
| --- | ---: | --- |
| Random walker + stop rule | 2.42 m | Best; slight over-estimate |
| lasR region growing | 2.62 m | Compact crowns, low bias |
| dalponte2016 | 2.70 m | Close behind |
| silva2016 | 2.79 m | Slightly over-grows |
| Random walker (no stop rule) | 3.19 m | Over-grows |
| Watershed (marker-free) | 3.42 m | Over-grows the most |

The headline: the four leading methods sit within about 0.3 m of each other.
The single biggest lever is not *which* method you choose, but *whether the
method has a sensible stop rule*. The random walker went from worst (3.19 m) to
best (2.42 m) simply by adding a rule that truncates each crown below a fraction
of its own top's height — stopping it from tiling the empty gaps between trees.

The same ranking holds on the widest-axis width, just with larger numbers
(3.57 m for the best method, up to 5.61 m for the marker-free watershed).

## 5. Two Robustness Findings

### Better Tree Tops Give Better Crowns

When the crown methods were re-seeded from a stronger tree-top detector
(`multichm`, the best classical detector in conifer forests), crown-width
accuracy improved on both methods that accept external seeds — and they reached
about 1.5× more trees (225 → 330) by finding crowns the simpler detector
missed. Detection and crown drawing are separate steps, but detection quality
clearly carries through into crown accuracy. The practical takeaway: invest in
the detection step; it pays off downstream.

### Crown Width Survives Sparse Data

The LiDAR was thinned from native density down to 1 point/m^2. Crown-width
accuracy **did not degrade** — it stayed flat or even improved slightly as the
data got sparser (the coarser surface smooths the crown outline toward the
field width). This is a notable contrast with understory *tree counting*, which
collapses over the same range. Once the canopy surface is resolved well enough
to see the dominant crowns, measuring their width tolerates far sparser data
than counting hidden trees does. A crown-width product can run on cheaper,
lower-density data than a stem-census product.

## 6. The 3-D Point Methods: a Complementary Profile

The three 3-D point-cloud methods (Li 2012, ptrees, AMS3D) have a different and
in some ways complementary behavior:

| Method | Average error (equiv. width) | Notable behavior |
| --- | ---: | --- |
| AMS3D | 2.97 m | Slightly *under*-estimates width; reaches many more trees |
| ptrees | 3.03 m | Slightly *under*-estimates width; reaches many more trees |
| Li 2012 | 5.23 m | Over-segments badly; not recommended |

Two useful points:

- AMS3D and ptrees **under**-estimate crown width, the opposite of the
  height-map methods' tendency to over-grow. And they match far more trees
  (about 410–475 vs 225) because they reach down to sub-canopy crowns the
  height-map methods cannot see.
- On the harder **widest-axis** measure, AMS3D and ptrees actually beat the
  height-map methods, because tracing the actual 3-D points tracks the widest
  reach of a crown better than the smoothed surface outline does.
- Li 2012 splits single trees into too many pieces here and is not recommended
  for crown width.

## 7. What Stays Hard

No method achieved a genuinely good fit on the **widest-axis** measure — every
one tends to draw crowns somewhat too large, especially in the gaps between
neighboring trees. The marker-controlled (seeded) watershed reduces this
over-grow compared with the marker-free version, but still trails the
grow-from-the-top methods. So "crowns drawn too big" remains the dominant error
across the board, and the widest-axis definition is the harder of the two
targets.

Accuracy also degrades for suppressed and intermediate (sub-canopy) trees,
where a small real crown sits beneath a larger one and the method claims too
much canopy for it. Those classes are also thinly sampled here, so treat their
numbers as indicative.

## 8. Practical Recommendations

For an **equivalent-width** crown product:

- Use a height-map method with a good stop rule — the random walker with a
  height cutoff, or lasR region growing — seeded from the strongest available
  detector. These hold accuracy down to about 1 point/m^2.

For a **widest-axis** product, or when **sub-canopy reach** matters:

- Add a 3-D point method (AMS3D or ptrees) as a complementary arm. They reach
  more trees and track the widest axis better, at the cost of a small
  under-estimate on the equivalent-width measure.

General guidance:

1. Always report which diameter definition you are using, and compare against
   the matching field column.
2. Expect crown width to be robust to point density — you can map width on
   cheaper, sparser data than you need for counting understory trees.
3. Expect the main error to be over-grown crowns; a stop rule is the cheapest
   fix.
4. Crown delineation is a separate downstream problem from tree-top detection.
   Strong detection helps, but report the two accuracies separately.

## Source Result Files

This report summarizes:

- `crown-segmentation-results.md`

It deliberately excludes the deep-learning / machine-learning crown arm
(TreeisoNet) documented in that file, per the non-ML scope of this version.
