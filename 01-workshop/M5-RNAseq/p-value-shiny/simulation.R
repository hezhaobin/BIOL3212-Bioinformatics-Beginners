# Simulation and decision functions, independent of the Shiny interface.
# Each row is a gene; each column is an independent replicate.
# Inputs use population VARIANCE. rnorm() requires SD, hence sqrt() below.
simulate_genes <- function(proportion, control_mean, variance, treated_mean, n_genes = 5000, replicates = 10) {
  # Randomly choose which genes receive the specified treated distribution.
  assigned <- rep(FALSE, n_genes)
  n_assigned <- round(n_genes * proportion / 100)
  if (n_assigned > 0) assigned[sample.int(n_genes, n_assigned)] <- TRUE
  # The null being tested is equality of means, defined by unequal population means.
  truth <- assigned & (control_mean != treated_mean)
  # Draw all control values, then transform standard normals for treated values.
  # R recycles the gene-specific means and SDs across replicate columns.
  control <- matrix(rnorm(n_genes * replicates, control_mean, sqrt(variance)), ncol = replicates)
  treated <- matrix(rnorm(n_genes * replicates), ncol = replicates)
  treated <- treated * sqrt(variance) +
    ifelse(assigned, treated_mean, control_mean)
  # Estimate each gene's group means and unbiased sample variances (n - 1).
  mean_c <- rowMeans(control)
  mean_t <- rowMeans(treated)
  var_c <- rowSums((control - mean_c)^2) / (replicates - 1)
  var_t <- rowSums((treated - mean_t)^2) / (replicates - 1)
  # Pooled two-sided Student tests, equivalent to t.test(..., var.equal = TRUE).
  # With equal sample sizes, pooled variance is the average sample variance.
  # The population variance is shared, but each gene estimates it from its data.
  pooled_variance <- (var_c + var_t) / 2
  se2 <- pooled_variance * 2 / replicates
  statistic <- (mean_t - mean_c) / sqrt(se2)
  df <- 2L * replicates - 2L
  p <- 2 * pt(-abs(statistic), df)
  # Adjust the complete family of 5,000 p-values, preserving row order.
  results <- data.frame(
    Gene = sprintf("Gene_%04d", seq_len(n_genes)),
    True_DEG = truth, Control_mean = mean_c, Treated_mean = mean_t,
    Difference = mean_t - mean_c, t_statistic = statistic, df = df,
    Raw_p = p, Bonferroni_p = p.adjust(p, "bonferroni"),
    BH_adjusted_p = p.adjust(p, "BH")
  )
  list(control = control, treated = treated, results = results)
}

# Reanalyze the same simulated measurements; selecting a test never draws data.
# The t-test estimates each gene's variance. The z-test uses the population
# variance stored with this experiment, not any subsequently edited UI value.
analyze_genes <- function(experiment, variance, replicates, test = "t") {
  results <- experiment$results
  if (test == "z") {
    results$t_statistic <- NULL
    results$df <- NULL
    results$z_statistic <- results$Difference / sqrt(2 * variance / replicates)
    results$Raw_p <- 2 * pnorm(-abs(results$z_statistic))
    results$Bonferroni_p <- p.adjust(results$Raw_p, "bonferroni")
    results$BH_adjusted_p <- p.adjust(results$Raw_p, "BH")
    results <- results[, c("Gene", "True_DEG", "Control_mean", "Treated_mean",
                          "Difference", "z_statistic", "Raw_p", "Bonferroni_p", "BH_adjusted_p")]
  }
  results
}

# Compare the chosen statistic against the cutoff using the known population truth.
classify_genes <- function(results, method, threshold) {
  called <- results[[method]] <= threshold
  truth <- results$True_DEG
  c(TP = sum(called & truth), FP = sum(called & !truth),
    TN = sum(!called & !truth), FN = sum(!called & truth))
}

