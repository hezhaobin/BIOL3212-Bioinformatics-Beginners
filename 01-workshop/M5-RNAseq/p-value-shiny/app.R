# Shiny interface for the simulation functions in simulation.R.
# Run this directory with shiny::runApp(); only the shiny package is required.
library(shiny)
source("simulation.R", local = TRUE)

# Map student-facing correction labels to columns in the simulation results.
method_names <- c("Raw p-value" = "Raw_p", "Bonferroni-adjusted p-value" = "Bonferroni_p",
                  "BH-adjusted p-value (FDR control)" = "BH_adjusted_p")
# UI: population controls on the left, exploration and classification on the right.
ui <- fluidPage(
  tags$head(tags$style(HTML("body {background:#f7f9fc;color:#203047;font-size:14px} .form-control,.btn,.selectize-input,.selectize-dropdown {font-size:inherit} .well {background:white} .tab-content {background:white;padding:20px} .warning {background:#fff0cc;padding:12px;margin:12px 0} table {font-variant-numeric:tabular-nums} h2,h3 {color:#164e63}"))),
  titlePanel("5,000 genes: chance, discoveries, and false positives"),
  p("Simulate log2 expression with equal numbers of control and treated replicates per gene. The population truth is known here, so we can check our discoveries."),
  sidebarLayout(
    sidebarPanel(
      sliderInput("proportion", "Genes assigned to the treatment distribution (%)", 0, 100, 10, step = 1),
      numericInput("replicates", "Replicates per group", 6, min = 2, max = 50, step = 1),
      numericInput("control_mean", "Control mean (log2)", 2, step = 0.1),
      numericInput("variance", "Shared population variance", 1, min = 0.01, step = 0.1),
      numericInput("treated_mean", "Treated mean for assigned genes (log2)", 4, step = 0.1),
      p("All genes and both groups share this variance. For unassigned genes, both groups use the control mean."),
      uiOutput("population_warning"),
      actionButton("simulate", "Simulate experiment", class = "btn-primary"),
      p("Settings apply when you simulate. Changing the decision rule in Classification reuses the same data.")
    ),
    mainPanel(
      textOutput("experiment_description"),
      selectInput("test", "Statistical test (uses the same simulated data)",
                  c("Student t-test: estimated variance" = "t",
                    "Z-test: known population variance" = "z"), selected = "t"),
      p("The t-test estimates variance from each gene. The z-test uses the known simulation variance. Sample size affects the precision of both tests."),
      tabsetPanel(
        tabPanel("Explore a gene",
          uiOutput("gene_choices"), actionButton("pick", "Pick another"),
          plotOutput("gene_plot", height = "330px"), div(style = "overflow-x:auto", tableOutput("gene_stats")),
          plotOutput("group_mean_density", height = "280px"),
          plotOutput("difference_density", height = "280px"),
          p("These densities use all genes in the selected truth class. The difference is treated − control."),
          uiOutput("rejection_note"),
          p("A p-value is the probability, assuming equal population means and the test model, of a test statistic at least as extreme as observed (in either direction). It is not the probability that the null hypothesis is true."),
          plotOutput("p_hist", height = "260px"),
          p("Under a valid null model, raw p-values are approximately uniform. At a raw threshold of 0.05, about 5% of true nulls are rejected over repeated sampling; each experiment varies.")
        ),
        tabPanel("Classification",
      selectInput("method", "Decision statistic", method_names),
      numericInput("threshold", "Reject when statistic ≤ threshold", 0.05, min = 0, max = 1, step = 0.01),
      h3("This experiment: discoveries versus truth"),
      plotOutput("classification_counts", height = "300px"),
      plotOutput("classification_proportions", height = "340px"),
      tableOutput("classification_metrics"),
      p("Undefined ratios are shown as N/A."),
      textOutput("observed"),
      p("A negative decision means failure to reject the null. Power is TP / (TP + FN); type-II error rate is FN / (TP + FN)."),
      p("Bonferroni targets FWER; Benjamini–Hochberg (BH) targets FDR. A threshold is a long-run target, not a guaranteed outcome for one experiment. These simulations assume independent genes. With no true DEGs, every discovery is false, so FDR and FWER coincide."),
      tags$small("Teaching model: independent normal log2 expression, not a count-based RNA-seq analysis. Pooled Student t-tests assume normal measurements and a common population variance; df = 2n − 2.")
        ),
        tabPanel("DEG output",
          p("All 5,000 genes, sorted by raw p-value. Differences and test statistics use treated − control. Results use the selected t- or z-test. BH values are adjusted p-values used to control FDR."),
          numericInput("page", "Page (50 genes per page)", 1, min = 1, max = 100, step = 1),
          downloadButton("download", "Download all genes (CSV)"),
          div(style = "overflow-x:auto", tableOutput("results"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # Validate editable population inputs. Values represent variance, not SD.
  settings <- reactive({
    values <- c(input$proportion, input$control_mean, input$variance,
                input$treated_mean, input$replicates)
    validate(need(length(values) == 5 && all(is.finite(values)), "Enter all simulation settings."),
             need(input$variance > 0, "Variances must be positive."),
             need(input$proportion >= 0 && input$proportion <= 100, "Proportion must be between 0 and 100%."))
    validate(need(input$replicates >= 2 && input$replicates <= 50 &&
                    input$replicates == floor(input$replicates),
                  "Replicates must be a whole number from 2 to 50."))
    list(replicates = input$replicates, proportion = input$proportion, control_mean = input$control_mean,
         variance = input$variance, treated_mean = input$treated_mean)
  })
  # Decision inputs are separate from population settings, so changing the
  # correction or cutoff does not generate a new single experiment.
  cutoff <- reactive({
    validate(need(length(input$threshold) == 1 && is.finite(input$threshold) &&
                    input$threshold >= 0 && input$threshold <= 1, "Threshold must be between 0 and 1."))
    input$threshold
  })
  # Snapshot settings and data on launch or an explicit simulation click.
  experiment <- eventReactive(input$simulate, {
    s <- settings()
    list(settings = s, data = do.call(simulate_genes, s))
  }, ignoreNULL = FALSE)
  # Switching tests recomputes statistics and adjusted p-values without resampling.
  analysis <- reactive({
    req(input$test)
    e <- experiment()
    analyze_genes(e$data, e$settings$variance, e$settings$replicates, input$test)
  })
  # Raw two-sided alpha = 0.05 boundaries, always centered at the null difference 0.
  # They do not depend on truth class or the classification correction/cutoff.
  z_boundaries <- reactive({
    s <- experiment()$settings
    c(-1, 1) * qnorm(0.975) * sqrt(2 * s$variance / s$replicates)
  })
  output$rejection_note <- renderUI({
    if (input$test == "z")
      p("Dashed lines mark the theoretical null rejection boundaries for the known-variance z-test, raw two-sided alpha = 0.05. Differences outside the lines reject the null. These are not Bonferroni or BH boundaries, and remain centered on zero for true DEGs.")
    else
      p("The t-test has gene-specific boundaries on this difference scale because each gene estimates its own variance; fixed rejection lines are shown only for the known-variance z-test.")
  })
  # Equal population means imply no true DEGs.
  output$population_warning <- renderUI({
    s <- settings()
    if (s$control_mean == s$treated_mean)
      div(class = "warning", "Equal population means: there are no true DEGs. All positive calls will be false positives.")
  })
  output$experiment_description <- renderText({
    e <- experiment(); s <- e$settings
    paste0("Displayed experiment: ", s$replicates, " replicates per group; control mean ", s$control_mean, ", variance ", s$variance,
           "; assigned treated mean ", s$treated_mean, ", shared variance ", s$variance,
           "; ", sum(e$data$results$True_DEG), " true DEGs / 5,000. Population log2 difference: ",
           s$treated_mean - s$control_mean, " (", signif(2^(s$treated_mean - s$control_mean), 3),
           "× on the original scale, comparing geometric means).",
           if (!identical(s, settings())) " Settings changed: click Simulate experiment to apply." else "")
  })
  # Offer only truth classes present in this experiment (handles 0% and 100%).
  output$gene_choices <- renderUI({
    truth <- experiment()$data$results$True_DEG
    choices <- c(if (any(truth)) "True DEG", if (any(!truth)) "True Null")
    selectInput("gene_type", "Sample gene from", choices)
  })
  # Cache one random row until the experiment, truth class, or pick button changes.
  gene <- reactive({
    input$pick
    truth <- experiment()$data$results$True_DEG
    req(input$gene_type)
    candidates <- which(truth == (input$gene_type == "True DEG"))
    req(length(candidates) > 0)
    candidates[sample.int(length(candidates), 1)]
  })
  # Fixed horizontal offsets keep observations visible without jitter.
  # Scale all plot text by 25%; explicit legend sizes are relative to this base.
  output$gene_plot <- renderPlot({
    par(cex = 1.25)
    d <- experiment()$data; i <- gene()
    values <- c(d$control[i, ], d$treated[i, ])
    n <- ncol(d$control)
    offsets <- seq(-0.15, 0.15, length.out = n)
    plot(rep(1:2, each = n) + rep(offsets, 2), values,
         xlim = c(0.5, 2.5), xaxt = "n", xlab = "", ylab = "Expression (log2)",
         pch = 19, cex = 1.5 / 1.25, col = rep(c("#167d9a", "#c05c32"), each = n),
         main = paste(d$results$Gene[i], "—", input$gene_type))
    axis(1, 1:2, c("Control", "Treated"))
    segments(c(0.8, 1.8), c(mean(d$control[i, ]), mean(d$treated[i, ])),
             c(1.2, 2.2), c(mean(d$control[i, ]), mean(d$treated[i, ])), lwd = 4)
    legend("topright", "Horizontal bars: sample means", bty = "n", cex = 0.85)
  })
  # The density panels share the truth selector, but not the picked example gene.
  selected_results <- reactive({
    req(input$gene_type)
    r <- analysis()
    r <- r[r$True_DEG == (input$gene_type == "True DEG"), ]
    validate(need(nrow(r) >= 2, "At least two genes in this truth class are needed for a density."))
    r
  })
  output$group_mean_density <- renderPlot({
    par(cex = 1.25)
    r <- selected_results()
    dc <- density(r$Control_mean); dt <- density(r$Treated_mean)
    plot(dc, col = "#167d9a", lwd = 2, xlim = range(dc$x, dt$x),
         ylim = c(0, max(dc$y, dt$y)), xlab = "Sample mean expression (log2)",
         main = paste("Gene sample means —", input$gene_type))
    lines(dt, col = "#c05c32", lwd = 2)
    legend("topright", c("Control", "Treated"), col = c("#167d9a", "#c05c32"), lwd = 2, bty = "n")
  })
  output$difference_density <- renderPlot({
    par(cex = 1.25)
    r <- selected_results()
    curve <- density(r$Difference)
    limits <- range(curve$x)
    if (input$test == "z") limits <- range(limits, z_boundaries())
    plot(curve, lwd = 2, col = "#167d9a", xlim = limits,
         xlab = "Treated − control sample mean (log2)",
         main = paste("Differences in sample means —", input$gene_type))
    if (input$test == "z") {
      abline(v = z_boundaries(), lty = 2, lwd = 2, col = "#c05c32")
      legend("topright", "Null rejection boundaries (raw alpha = 0.05)",
             col = "#c05c32", lty = 2, lwd = 2, bty = "n", cex = 0.8)
    }
  })
  # Round only for display; the downloadable CSV retains numeric precision.
  format_results <- function(x) {
    x[] <- lapply(x, function(column) if (is.numeric(column)) format(signif(column, 4), trim = TRUE) else column)
    x
  }
  output$gene_stats <- renderTable(format_results(analysis()[gene(), ]), digits = 4)
  output$p_hist <- renderPlot({
    par(cex = 1.25)
    hist(analysis()$Raw_p, breaks = seq(0, 1, 0.05),
         col = "#93cbd7", border = "white", main = "Raw p-values across 5,000 genes",
         xlab = "Raw p-value", ylab = "Number of genes")
  })
  # Pagination changes the view only; sorting always covers all 5,000 genes.
  sorted <- reactive({
    r <- analysis()
    r[order(r$Raw_p), ]
  })
  output$results <- renderTable({
    req(input$page)
    page <- max(1, min(100, floor(input$page)))
    format_results(sorted()[seq.int((page - 1) * 50 + 1, page * 50), ])
  }, digits = 4, rownames = FALSE)
  output$download <- downloadHandler("gene-results.csv", function(file) write.csv(sorted(), file, row.names = FALSE))
  # Recompute the 2 x 2 table from the fixed data and the current decision rule.
  counts <- reactive(classify_genes(analysis(), input$method, cutoff()))
  # Classification views reuse the fixed experiment when the decision rule changes.
  output$classification_counts <- renderPlot({
    par(cex = 1.25)
    n <- counts()
    bars <- rbind(c(n["TP"], 0), c(n["FP"], 0), c(0, n["FN"]))
    colnames(bars) <- c("Called DEG", "Missed true DEGs")
    positions <- barplot(bars, col = c("#167d9a", "#c05c32", "#8b6bb1"),
                        ylim = c(0, max(1, colSums(bars)) * 1.2),
                        ylab = "Number of genes", main = "Discoveries and false negatives")
    text(positions, colSums(bars), labels = colSums(bars), pos = 3)
    legend("topright", c("True positives", "False positives", "False negatives"),
           fill = c("#167d9a", "#c05c32", "#8b6bb1"), bty = "n")
  })
  # Each pie conditions on the decision: positive calls or negative calls.
  # Empty decision sets have no defined composition and are labeled explicitly.
  output$classification_proportions <- renderPlot({
    par(cex = 1.25)
    n <- counts()
    old <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
    par(cex = 1.25)
    on.exit(par(old))
    draw_pie <- function(values, labels, colors, heading) {
      if (sum(values) == 0) {
        plot.new(); title(heading)
        text(0.5, 0.5, "No genes in this decision group", cex = 0.9)
      } else {
        # Keep labels for zero-sized categories in the legend, too.
        captions <- paste0(labels, ": ", values, " (",
                           round(100 * values / sum(values), 1), "%)")
        pie(values, labels = NA, col = colors, main = heading, radius = 0.72)
        legend("bottom", legend = captions, fill = colors, bty = "n", cex = 0.85)
      }
    }
    draw_pie(n[c("TP", "FP")], c("True positive", "False positive"),
             c("#167d9a", "#c05c32"), "Called DEGs")
    draw_pie(n[c("TN", "FN")], c("True negative", "False negative"),
             c("#548858", "#8b6bb1"), "Called non-DEGs")
  })
  output$classification_metrics <- renderTable({
    n <- counts()
    ratio <- function(a, b) if (b == 0) "N/A" else paste0(round(100 * a / b, 1), "%")
    data.frame(Metric = c("Precision", "Recall / power", "Specificity"),
               Definition = c("TP / (TP + FP): fraction of called DEGs that are true DEGs",
                              "TP / (TP + FN): fraction of true DEGs detected",
                              "TN / (TN + FP): fraction of true nulls not called DEG"),
               Value = c(ratio(n["TP"], n["TP"] + n["FP"]),
                         ratio(n["TP"], n["TP"] + n["FN"]),
                         ratio(n["TN"], n["TN"] + n["FP"])))
  })
  output$observed <- renderText({
    n <- counts(); discoveries <- n["TP"] + n["FP"]
    paste0("Discoveries: ", discoveries, ". Observed false discovery proportion: ",
           round(100 * n["FP"] / max(1, discoveries), 1), "%. At least one FP: ",
           if (n["FP"] > 0) "yes." else "no.")
  })
}

shinyApp(ui, server)
