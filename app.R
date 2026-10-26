# Run in RStudio with Run App, or: shiny::runApp("build/week08/mod09/clt-dashboard")
# Only dependency: install.packages("shiny")
library(shiny)

# Each preset supplies an iid generator and its exact population moments.
population <- function(kind, values = "0, 1, 5, 20", weights = "70, 20, 9, 1") {
  switch(kind,
    exponential = list(draw = function(k) rexp(k), mu = 1, sigma = 1,
      note = "Right-skewed: exponential with mean 1."),
    bimodal = list(draw = function(k) sample(c(-3, 3), k, TRUE) + rnorm(k, 0, .6),
      mu = 0, sigma = sqrt(9 + .6^2), note = "Two peaks: equal mixtures centered at -3 and 3."),
    ushape = list(draw = function(k) rbeta(k, .3, .3), mu = .5,
      sigma = sqrt(.3^2 / (.6^2 * 1.6)), note = "U-shaped: most observations are near 0 or 1."),
    rare = list(draw = function(k) 30 * rbinom(k, 1, .01) + rnorm(k),
      mu = .3, sigma = sqrt(1 + 30^2 * .01 * .99),
      note = "Rare outliers: 99% from N(0,1), 1% from N(30,1). Even n = 50 can look far from normal."),
    lognormal = list(draw = function(k) rlnorm(k, 0, 1.4), mu = exp(1.4^2 / 2),
      sigma = sqrt((exp(1.4^2) - 1) * exp(1.4^2)),
      note = "Strong right skew: lognormal with log-scale SD 1.4. Convergence can be slow."),
    bernoulli = list(draw = function(k) rbinom(k, 1, .02), mu = .02,
      sigma = sqrt(.02 * .98), note = "Rare successes: each observation is 1 with probability 0.02, otherwise 0."),
    normal = list(draw = function(k) rnorm(k), mu = 0, sigma = 1,
      note = "Normal population: the sample mean is exactly normal for every sample size."),
    custom = {
      parse_numbers <- function(x) suppressWarnings(as.numeric(strsplit(trimws(x), "[,;[:space:]]+")[[1]]))
      v <- parse_numbers(values); w <- parse_numbers(weights)
      if (length(v) < 2 || length(v) != length(w) || any(!is.finite(c(v, w))) ||
          any(w < 0) || sum(w) <= 0) stop("Enter at least two values and matching nonnegative weights with a positive total.")
      w <- w / sum(w); mu <- sum(v * w); sig <- sqrt(sum(w * (v - mu)^2))
      if (!is.finite(sig) || sig <= 0) stop("Give positive weight to at least two different values.")
      list(draw = function(k) sample(v, k, TRUE, prob = w), mu = mu, sigma = sig,
        note = "Custom discrete population. Weights are normalized to probabilities.")
    }, stop("Choose a population.")
  )
}

simulate_clt <- function(pop, sizes, B, seed) {
  set.seed(seed)
  preview <- pop$draw(40000L)
  # Each row is one independent sample. Its mean is one plotted observation.
  means <- lapply(sizes, function(n) rowMeans(matrix(pop$draw(B * n), nrow = B, ncol = n)))
  list(pop = pop, sizes = sizes, B = B, preview = preview, means = means)
}

# Histogram density uses ALL observations in the denominator, including off-screen ones.
# This avoids silently renormalizing a cropped heavy tail.
draw_hist <- function(x, limits, bins, title, xlabel, color, target = NULL, center = NULL) {
  edges <- seq(limits[1], limits[2], length.out = bins + 1L)
  bin <- findInterval(x, edges, rightmost.closed = TRUE)
  counts <- tabulate(bin[bin >= 1L & bin <= bins], nbins = bins)
  heights <- counts / (length(x) * diff(edges))
  grid <- seq(limits[1], limits[2], length.out = 500L)
  normal_y <- if (is.null(target)) 0 else dnorm(grid, target[1], target[2])
  plot(NA, xlim = limits, ylim = c(0, max(heights, normal_y, .001) * 1.12),
    xlab = xlabel, ylab = "Density", main = title, bty = "l", yaxs = "i")
  rect(head(edges, -1), 0, tail(edges, -1), heights, col = color, border = "white")
  if (!is.null(target)) lines(grid, normal_y, col = "#2E6E6B", lwd = 2.5)
  if (!is.null(center)) abline(v = center, col = "#80273F", lty = 2, lwd = 1.5)
  hidden <- sum(x < limits[1] | x > limits[2])
  if (hidden > 0) mtext(sprintf("Outside view: %s / %s observations", hidden, length(x)), side = 3, cex = .65, col = "#555555")
}

plot_clt <- function(sim, standardized = FALSE, overlay = TRUE, zoom = TRUE, bins = 45) {
  old <- par(mfrow = c(2, 2), mar = c(4, 4, 3.6, 1.2), oma = c(0, 0, .5, 0),
    family = "sans", col.axis = "#333333", col.lab = "#333333", cex.main = 1.15)
  on.exit(par(old))
  p <- sim$pop
  raw_limits <- if (zoom) as.numeric(quantile(sim$preview, c(.001, .999))) else range(c(sim$preview, unlist(sim$means)))
  if (diff(raw_limits) == 0) raw_limits <- range(c(sim$preview, unlist(sim$means)))
  if (diff(raw_limits) == 0) raw_limits <- raw_limits + c(-1, 1)
  raw_limits <- raw_limits + c(-1, 1) * diff(raw_limits) * .04
  draw_hist(sim$preview, raw_limits, bins, "Population draws", "Individual observation X", "#D6C0C8", center = p$mu)
  for (j in seq_along(sim$sizes)) {
    n <- sim$sizes[j]; x <- sim$means[[j]]; se <- p$sigma / sqrt(n)
    if (standardized) {
      x <- (x - p$mu) / se
      limits <- if (zoom) c(-4, 4) else range(c(-4, 4, unlist(Map(function(m, n) (m-p$mu)/(p$sigma/sqrt(n)), sim$means, sim$sizes))))
      target <- c(0, 1); center <- 0; label <- "Z = (sample mean - population mean) / SE"
    } else {
      limits <- raw_limits; target <- c(p$mu, se); center <- p$mu; label <- "Sample mean"
    }
    draw_hist(x, limits, bins, paste("Sample means: n =", n), label,
      c("#D6C0C8", "#B9D4D1", "#DFD2B5")[j], if (overlay) target else NULL, center)
  }
}

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body {font-family: 'Helvetica Neue', Arial, sans-serif; color:#262626; background:white;}
    .container-fluid {max-width:1500px; margin:auto; padding:24px;}
    h1 {color:#80273F; font-size:32px; font-weight:600; margin-top:0;}
    h3 {font-size:18px; font-weight:600;}
    .well {background:#f7f2f4; border:0; border-radius:6px; box-shadow:none;}
    .btn-primary {background:#80273F; border-color:#80273F;}
    .btn-primary:hover,.btn-primary:focus {background:#642035; border-color:#642035;}
    .help-block {color:#555; font-size:13px;}
    .note {padding:12px 16px; background:#eef4f3; margin-bottom:12px; border-radius:5px;}
    .shiny-output-error-validation {color:#80273F; padding:12px;}
    :focus-visible {outline:3px solid #2E6E6B; outline-offset:2px;}
    @media(max-width:768px){.container-fluid{padding:12px;} h1{font-size:26px;}}
  "))),
  h1("From random samples to a sampling distribution"),
  p("Keep the population fixed. Change the number of observations in each sample."),
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput("kind", "Population", choices = c("Two peaks" = "bimodal", "Right-skewed" = "exponential",
        "U-shaped" = "ushape", "Rare outliers" = "rare", "Strong right skew" = "lognormal",
        "Rare 0/1 outcomes" = "bernoulli", "Normal" = "normal", "Make your own" = "custom")),
      conditionalPanel("input.kind == 'custom'",
        textInput("values", "Possible values", "0, 1, 5, 20"),
        textInput("weights", "Weights", "70, 20, 9, 1"),
        helpText("Separate entries with commas. Weights need not sum to 1.")),
      h3("Observations per sample"),
      fluidRow(column(4, numericInput("n1", "First", 1, 1, 500, 1)),
        column(4, numericInput("n2", "Second", 5, 1, 500, 1)),
        column(4, numericInput("n3", "Third", 50, 1, 500, 1))),
      sliderInput("B", "Number of repeated samples", min = 500, max = 10000, value = 3000, step = 500, sep = ","),
      checkboxInput("standardized", "Standardize the sample means", FALSE),
      checkboxInput("overlay", "Show the normal comparison", TRUE),
      checkboxInput("zoom", "Zoom in (report observations outside view)", TRUE),
      sliderInput("bins", "Histogram bins", 15, 80, 45),
      numericInput("seed", "Random seed", 1151, 1, 1000000, 1),
      actionButton("draw", "Draw new samples", class = "btn-primary"),
      br(), br(), downloadButton("png", "Save plots"), downloadButton("csv", "Save means")
    ),
    mainPanel(width = 9,
      uiOutput("description"),
      plotOutput("plots", height = "650px"),
      p("Bars: simulated observations. Teal curve: normal comparison. Dashed line: population mean, or 0 after standardizing."),
      tableOutput("summary"),
      tags$details(tags$summary("What changes when you change n?"),
        p("Each sample contains n independent observations from the same population. Its mean is one observation in a sample-mean histogram."),
        p("The CLT applies here because all offered populations have finite, positive variance. For sufficiently large n, the standardized sample mean is approximately normal. The population itself does not change."),
        p("More repeated samples smooth the histogram. Increasing n changes the sampling distribution. These are different operations."),
        p("Raw sample-mean panels share the population's horizontal scale. Standardized panels share a Z scale. Zooming changes only the view; every observation remains in the calculations."),
        p("Discrete populations remain discrete at finite n; histogram bin choices can emphasize gaps. The population panel is a fixed-seed simulation of 40,000 observations, not an exact density."))
    )
  )
)

server <- function(input, output, session) {
  pop <- reactive({
    result <- tryCatch(population(input$kind, input$values, input$weights), error = function(e) e)
    validate(need(!inherits(result, "error"), if (inherits(result, "error")) conditionMessage(result) else ""))
    result
  })
  simulation <- reactive({
    sizes <- c(input$n1, input$n2, input$n3)
    validate(need(length(sizes) == 3 && all(is.finite(sizes)) && all(sizes >= 1 & sizes <= 500) && all(sizes == floor(sizes)), "Use three whole sample sizes between 1 and 500."))
    validate(need(length(input$seed) == 1 && is.finite(input$seed) && input$seed >= 1 && input$seed <= 1000000, "Use a seed between 1 and 1,000,000."))
    p <- pop()
    # Population display stays fixed between redraws; change only the sampled means.
    sim <- simulate_clt(p, sizes, input$B, as.integer(input$seed + input$draw))
    set.seed(as.integer(input$seed)); sim$preview <- p$draw(40000L)
    sim
  })
  output$description <- renderUI({p <- pop(); div(class = "note", p$note,
    br(), sprintf("Population mean = %.3f. Population SD = %.3f.", p$mu, p$sigma))})
  output$plots <- renderPlot({plot_clt(simulation(), input$standardized, input$overlay, input$zoom, input$bins)}, res = 110)
  output$summary <- renderTable({
    s <- simulation()
    data.frame(`Sample size` = s$sizes, `Repeated samples` = s$B,
      `Average of sample means` = vapply(s$means, mean, numeric(1)),
      `SD of sample means` = vapply(s$means, sd, numeric(1)),
      `Theoretical SE` = s$pop$sigma / sqrt(s$sizes), check.names = FALSE)
  }, digits = 3, striped = FALSE, bordered = FALSE)
  output$png <- downloadHandler(filename = function() "clt-plots.png", content = function(file) {
    png(file, width = 1800, height = 1250, res = 160)
    on.exit(dev.off()); plot_clt(simulation(), input$standardized, input$overlay, input$zoom, input$bins)
  })
  output$csv <- downloadHandler(filename = function() "clt-sample-means.csv", content = function(file) {
    s <- simulation()
    result <- do.call(rbind, lapply(seq_along(s$sizes), function(j) data.frame(
      panel = j, n = s$sizes[j], repetition = seq_len(s$B), sample_mean = s$means[[j]],
      standardized_mean = (s$means[[j]] - s$pop$mu)/(s$pop$sigma/sqrt(s$sizes[j])))))
    write.csv(result, file, row.names = FALSE)
  })
}
shinyApp(ui, server)
