# CanvasXpressBio 0.99.5

* New `makeDemoSE()` builds the vignette's demonstration `SummarizedExperiment`
  in one call, carrying the derived differential-expression table and the
  per-cohort paired t-test p-values in its `metadata()`, so the vignette no
  longer constructs data inline.
* New `makeBeeswarmDemo()` returns a dense single-feature object for showing a
  sina / beeswarm overlay at scale.
* `cxboxplot()` gains `points` (`"all"`, `"sina"`, `"jitter"`, `"none"`) to pack
  the observations into a beeswarm that stays legible with many points.
* `cxsurvival()` gains `tooltip` to carry extra `colData` columns, so hovering a
  step or a censored observation shows that subject's metadata.
* Vignette: build the dataset with `makeDemoSE()`; document the enriched cell
  hover and the view-reset control on the heatmap; surface per-subject metadata,
  the censored-observation rug and the risk table on the Kaplan-Meier curve; and
  add a beeswarm-at-scale example.

# CanvasXpressBio 0.99.4

* Fix `cxvolcano()`: the result matrix was transposed, so every gene collapsed
  onto two plotted points. Genes are now the variables and the fold-change /
  p-value are the axis columns, giving one point per gene.
* `cxvolcano()` now colours points by regulation (up / down / not significant),
  draws the fold-change and p-value threshold lines, and gains `label`,
  `label_genes`, `label_fc` and `label_p` to label significant or named genes.
* `cxplot()` gains `overlays`, `var_overlays` and `cluster` to draw `colData` /
  `rowData` annotations as overlay tracks and to cluster rows and columns.
* New `cxsurvival()` draws a Kaplan-Meier curve from `colData`.
* New `cxboxplot()` draws a paired boxplot of one feature with per-subject
  connecting lines, faceting and colouring from `colData`.
* Rewrote the vignette around a single `SummarizedExperiment` shown as a
  heatmap, volcano, Kaplan-Meier curve and paired boxplot, and added a section
  motivating the package relative to existing Bioconductor visualization tools.

# CanvasXpressBio 0.99.3

* Version bump to re-run the Bioconductor staging build (no code changes; 0.99.2 built
  and passed checks but failed to publish).

# CanvasXpressBio 0.99.2

* Vignette charts now render in the built HTML: require `canvasXpress (>= 1.70.3)`, which
  ships a CanvasXpress library that BiocStyle's self-contained HTML post-processing no longer
  corrupts (earlier versions were mangled on inlining, leaving the widgets blank). Verified by
  rendering the vignette against the CRAN 1.70.3 library.

# CanvasXpressBio 0.99.1

* Add maintainer ORCID iD to `Authors@R` (BiocCheck note).

# CanvasXpressBio 0.99.0

* Initial submission to Bioconductor.
* `cxplot()` generic renders a `SummarizedExperiment` (and extending classes such as
  `SingleCellExperiment` and `DESeqDataSet`) as an interactive CanvasXpress heatmap, with
  `colData` as sample annotations and `rowData` as feature annotations. `n=` keeps only the
  most variable features.
* `cxvolcano()` draws a volcano plot from a differential-expression result table
  (`DESeqResults`, `limma::topTable()`, or a data frame).
