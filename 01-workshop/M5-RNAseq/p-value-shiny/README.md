# P-values and multiple testing

A standalone Shiny app for BIOL-3212. Runtime dependency: `shiny` only.

## Run locally

From R, install Shiny once if necessary, then launch:

```r
install.packages("shiny")
shiny::runApp("01-workshop/M5-RNAseq/p-value-shiny")
```

The path above assumes the working directory is the course repository root.

## Deploy to shinyapps.io

Install `rsconnect`, configure your shinyapps.io account using the account's
provided `rsconnect::setAccountInfo()` command in your own R session, and run:

```r
rsconnect::deployApp(
  appDir = "01-workshop/M5-RNAseq/p-value-shiny",
  appName = "biol3212-p-values"
)
```

Keep account tokens out of these files and version control. Deploy this app
directory alone; the surrounding course datasets are not needed.

## Teaching and simulation choices

- Each experiment contains 5,000 independent genes and 2–50 independent normal
  measurements per group (default 10), on the log2 scale. Inputs specify population variance; the generator passes its square root
  as the SD to `rnorm()`.
- Assigned genes use the treated distribution; all other genes use the control
  distribution in both groups. A true DEG has unequal population means.
  All genes and both groups share one adjustable population variance.
  Equal population means produce no true DEGs; these genes are labeled True Null.
- Tests are two-sided pooled Student t-tests, with t and differences oriented
  as treated minus control. Calculations follow
  `t.test(treated, control, var.equal = TRUE)`, with df = 2n - 2.
  Under the independent normal, equal-variance null model, this is an exact test.
  The simulation is a teaching model, not DESeq2.
- Explore a gene contains sample-mean and mean-difference densities for the
  selected truth class. Density estimates need at least two genes in that class.
  Selecting the known-variance z-test adds dashed theoretical null rejection
  boundaries at ±qnorm(0.975) * sqrt(2 * variance / n) on the difference plot.
  These always use the stored experiment settings and null center zero, even
  for the True DEG class. They represent raw alpha = 0.05, independently of
  the classification threshold and multiple-testing correction.
- A t/z selector above the tabs reanalyzes the same measurements without
  resimulating or picking a new example gene. Default: pooled Student t-test.
  The z-test uses the known shared population variance and normal reference
  distribution; it reports z rather than t and has no degrees-of-freedom column.
  All tables, p-values, corrections, and classification plots use the chosen test.
  Known variance is available here because this is a simulation; sample size
  still controls the standard error in both tests.
- Classification contains the decision controls, discovery and false-negative
  counts, two composition pie charts (TP/FP among called DEGs and TN/FN among called
  non-DEGs), and precision, recall, and specificity with definitions. Undefined
  precision or recall is displayed as N/A. DEG output is the third tab.
- A new single experiment is generated on launch and on Simulate experiment.
  Picking a gene and changing correction methods do not resimulate its data.
- The output includes all genes, sorted by raw p-value, in pages of 50; CSV
  download preserves full numeric precision. On-screen values are rounded.
- The app intentionally uses a new random draw each time rather than a fixed
  seed. Independent browser sessions have separate experiment state.

Suggested sequence: start at 0% DEGs, compare raw and adjusted decisions,
then introduce 10–20% DEGs and compare FWER, FDR, and missed discoveries.

Statistical references: R documentation for
[t.test](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/t.test.html)
and [p.adjust](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html).
