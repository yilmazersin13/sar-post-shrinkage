# Simulation study (Section 5): Experiment 1 (adaptive-Lasso selection, 500 replications)
# and Experiment 2 (fixed restriction, 1,000 replications), Tables 1-2 and Figures 1-4.
# Run from the repository folder: source("simulation.R")
library(ggplot2)
source("sar_functions.R")

# Fixed design: x1 (p1 = 4 target columns) and x2 (p2 candidate columns) with
# Cov(x1, x2) = c / sqrt(p2/p1) in blocks; uniform margins, centred and scaled (divisor n).
make_block_design <- function(n, p1 = 4, p2 = 8, coupling = .6, seed = 720000) {
  stopifnot(p2 %% p1 == 0, coupling >= 0, coupling < 1)
  groups <- p2 / p1
  B <- coupling * matrix(rep(diag(p1), groups), p1, p2) / sqrt(groups)
  Sigma <- rbind(cbind(diag(p1), B), cbind(t(B), diag(p2)))
  set.seed(seed)
  Z <- matrix(runif(n * (p1 + p2), -sqrt(3), sqrt(3)), n, p1 + p2)
  X <- Z %*% chol(Sigma)
  X <- sweep(X, 2, colMeans(X), "-")
  X <- sweep(X, 2, sqrt(colMeans(X^2)), "/")
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  attr(X, "population_covariance") <- Sigma
  X
}

# Departure direction for beta2: two nonzero entries, the rest exactly zero.
sparse_direction <- function(p2, active = 2) {
  stopifnot(active > 0, active <= p2)
  c(rep(1 / sqrt(active), active), rep(0, p2 - active))
}

# Losses for beta1 and for average direct/indirect impacts, one row per method.
make_record <- function(estimates, beta, p1, rho_true, rho_common, rho_native,
                        weights, setting, replication, block) {
  truth <- outer(beta[seq_len(p1)], impact_factors(rho_true, weights))
  do.call(rbind, lapply(names(estimates), function(method) {
    error <- estimates[[method]][seq_len(p1)] - beta[seq_len(p1)]
    common <- outer(estimates[[method]][seq_len(p1)], impact_factors(rho_common, weights))
    native <- outer(estimates[[method]][seq_len(p1)], impact_factors(rho_native[method], weights))
    row <- data.frame(setting, replication = replication, block = block, method = method,
      loss = sum(error^2), rho_error = rho_common - rho_true,
      rho_native_error = rho_native[method] - rho_true,
      ADE_loss = sum((common[, "ADE"] - truth[, "ADE"])^2),
      ATE_loss = sum((common[, "ATE"] - truth[, "ATE"])^2),
      AIE_loss = sum((common[, "AIE"] - truth[, "AIE"])^2),
      ADE_native_loss = sum((native[, "ADE"] - truth[, "ADE"])^2),
      AIE_native_loss = sum((native[, "AIE"] - truth[, "AIE"])^2))
    for (j in seq_len(p1)) row[[paste0("error", j)]] <- error[j]
    rownames(row) <- NULL
    row
  }))
}

selection_rates <- function(selected, beta, p1) {
  strong <- seq_len(p1)
  weak <- which(abs(beta) > 0 & !seq_along(beta) %in% strong)
  zeros <- which(beta == 0)
  c(strong_recall = mean(strong %in% selected),
    all_strong_retained = as.numeric(all(strong %in% selected)),
    weak_recall = if (length(weak)) mean(weak %in% selected) else NA_real_,
    FPR = if (length(zeros)) mean(zeros %in% selected) else NA_real_,
    exact_support = as.numeric(setequal(selected, which(beta != 0))))
}

# Experiment 1 (adaptive-Lasso selection): beta = (1.5, -1, 2, 0.8, Delta* nu)'.
# Seed of replication r is seed + r. Failed fits are recorded, not replaced.
finite_cell <- function(setting, weights, X, R, nlambda, seed) {
  records <- list(); diagnostics <- list()
  p1 <- 4; p2 <- setting$p2
  beta <- c(1.5, -1, 2, .8, setting$delta_star * sparse_direction(p2))
  for (replication in seq_len(R)) {
    set.seed(seed + replication)
    y <- spatial_solve(weights, as.vector(X %*% beta) + rnorm(setting$n), setting$rho)
    diagnostic <- data.frame(setting, replication = replication, seed = seed + replication,
      fixed_ok = FALSE, selected_ok = FALSE, failure = "",
      T_fixed = NA_real_, T_selected = NA_real_, q = NA_integer_, fallback = NA,
      strong_recall = NA_real_, all_strong_retained = NA_real_, weak_recall = NA_real_,
      FPR = NA_real_, exact_support = NA_real_,
      AL_lambda = NA_real_, RIDGE_lambda = NA_real_, AL_edge = NA, RIDGE_edge = NA,
      AL_kkt = NA_real_, RIDGE_kkt = NA_real_, rho_boundary = NA)
    tryCatch({
      full <- sar_qmle(y, X, weights)
      C <- sar_information(X, full$beta, full$rho, full$sigma2, weights)$C
      restricted <- sar_qmle(y, X[, 1:p1, drop = FALSE], weights)
      fixed <- post_shrinkage(full, restricted, 1:p1, C, setting$n)
      native <- c(OF = full$rho, UF = restricted$rho, PT = full$rho, S = full$rho, PS = full$rho)
      records[[length(records) + 1]] <- make_record(fixed$estimates, beta, p1, setting$rho,
        full$rho, native, weights, setting, replication, "fixed")
      diagnostic$fixed_ok <- TRUE; diagnostic$T_fixed <- fixed$T
      diagnostic$rho_boundary <- full$boundary
      al <- tune_penalty(y, X, weights, full, "AL", nlambda)
      ridge <- tune_penalty(y, X, weights, full, "ridge", nlambda)
      selected_refit <- sar_qmle(y, X[, al$selected, drop = FALSE], weights)
      selected <- post_shrinkage(full, selected_refit, al$selected, C, setting$n)
      selected$estimates$AL <- al$beta; selected$estimates$RIDGE <- ridge$beta
      native <- c(OF = full$rho, UF = selected_refit$rho, PT = full$rho, S = full$rho,
                  PS = full$rho, AL = al$rho, RIDGE = ridge$rho)
      records[[length(records) + 1]] <- make_record(selected$estimates, beta, p1,
        setting$rho, full$rho, native, weights, setting, replication, "selected")
      diagnostic$selected_ok <- TRUE; diagnostic$T_selected <- selected$T
      diagnostic$q <- selected$q; diagnostic$fallback <- selected$fallback
      rates <- selection_rates(al$selected, beta, p1)
      for (name in names(rates)) diagnostic[[name]] <- rates[name]
      diagnostic$AL_lambda <- al$lambda; diagnostic$RIDGE_lambda <- ridge$lambda
      diagnostic$AL_edge <- al$edge; diagnostic$RIDGE_edge <- ridge$edge
      diagnostic$AL_kkt <- al$kkt; diagnostic$RIDGE_kkt <- ridge$kkt
    }, error = function(e) diagnostic$failure <<- conditionMessage(e))
    diagnostics[[replication]] <- diagnostic
  }
  list(raw = if (length(records)) do.call(rbind, records) else NULL,
       diagnostics = do.call(rbind, diagnostics), X = X, beta = beta)
}

# Experiment 2 (fixed restriction): beta2 = kappa / sqrt(n), kappa scaled so that
# kappa' C22.1 kappa = Delta, with C from the null coefficients and true nuisance values.
theory_cell <- function(setting, weights, X, R, seed) {
  p1 <- 4; p2 <- setting$p2
  beta <- c(1.5, -1, 2, .8, rep(0, p2))
  C0 <- sar_information(X, beta, setting$rho, 1, weights)$C
  baseline <- adr_formula(C0, rep(0, p2), p1)
  nu <- sparse_direction(p2)
  kappa_local <- sqrt(setting$Delta / as.numeric(crossprod(nu, baseline$C22 %*% nu))) * nu
  beta[(p1 + 1):(p1 + p2)] <- kappa_local / sqrt(setting$n)
  theory <- adr_formula(C0, kappa_local, p1)
  expected <- data.frame(setting[rep(1, length(theory$risk)), , drop = FALSE],
    method = names(theory$risk), ADR = as.vector(theory$risk),
    theoretical_bias2 = unname(vapply(theory$bias, function(z) sum(z^2), numeric(1))),
    row.names = NULL)
  raw <- list(); diagnostics <- list()
  for (replication in seq_len(R)) {
    set.seed(seed + replication)
    y <- spatial_solve(weights, as.vector(X %*% beta) + rnorm(setting$n), setting$rho)
    diagnostic <- data.frame(setting, replication = replication, seed = seed + replication,
      valid = FALSE, T = NA_real_, rejected = NA, failure = "")
    tryCatch({
      full <- sar_qmle(y, X, weights)
      restricted <- sar_qmle(y, X[, 1:p1], weights)
      C <- sar_information(X, full$beta, full$rho, full$sigma2, weights)$C
      post <- post_shrinkage(full, restricted, 1:p1, C, setting$n)
      native <- c(OF = full$rho, UF = restricted$rho, PT = full$rho, S = full$rho, PS = full$rho)
      raw[[length(raw) + 1]] <- make_record(post$estimates, beta, p1, setting$rho,
        full$rho, native, weights, setting, replication, "theory")
      diagnostic$valid <- TRUE; diagnostic$T <- post$T
      diagnostic$rejected <- post$T > qchisq(.95, p2)
    }, error = function(e) diagnostic$failure <<- conditionMessage(e))
    diagnostics[[replication]] <- diagnostic
  }
  list(raw = if (length(raw)) do.call(rbind, raw) else NULL,
    diagnostics = do.call(rbind, diagnostics), expected = expected,
    C0 = C0, kappa = kappa_local, direction = nu, X = X, beta = beta)
}

row_mcse <- function(z) {
  R <- ncol(z)
  if (R < 2) return(rep(NA_real_, nrow(z)))
  sqrt(rowSums(sweep(z, 1, rowMeans(z), "-")^2) / (R * (R - 1)))
}

# Sitewise errors of S(rho)^-1 X1 b against the truth, with Monte Carlo SEs.
spatial_summary <- function(cell, weights, block = "selected") {
  raw <- cell$raw[cell$raw$block == block, ]
  methods <- unique(raw$method)
  if (!length(methods)) return(NULL)
  target <- cell$beta[1:4]
  X1 <- cell$X[, 1:4]
  truth <- spatial_solve(weights, as.vector(X1 %*% target), raw$rho[1])
  errors <- list()
  for (method in methods) {
    z <- raw[raw$method == method, ]; z <- z[order(z$replication), ]
    errors[[method]] <- vapply(seq_len(nrow(z)), function(r) {
      b <- target + as.numeric(z[r, paste0("error", 1:4)])
      spatial_solve(weights, as.vector(X1 %*% b), z$rho[r] + z$rho_error[r]) - truth
    }, numeric(nrow(X1)))
  }
  of_loss <- errors$OF^2
  of_mse <- rowMeans(of_loss)
  do.call(rbind, lapply(methods, function(method) {
    e <- errors[[method]]; loss <- e^2; mse <- rowMeans(loss)
    ratio <- mse / of_mse
    paired <- loss - of_loss
    influence <- -100 * sweep(loss - sweep(of_loss, 1, ratio, "*"), 1, of_mse, "/")
    data.frame(id = raw$id[1], rho = raw$rho[1], delta_star = raw$delta_star[1],
      coupling = raw$coupling[1], n = nrow(X1), p2 = raw$p2[1], block = block,
      method = method, site = seq_len(nrow(X1)), x = weights$xy[, 2], y = weights$xy[, 1],
      valid = ncol(e), truth = truth, bias = rowMeans(e), MSE = mse, RMSE = sqrt(mse),
      MSE_MCSE = row_mcse(loss), paired_MSE_difference = rowMeans(paired),
      paired_difference_MCSE = row_mcse(paired), gain_percent = 100 * (1 - ratio),
      gain_MCSE = row_mcse(influence))
  }))
}

# LaTeX tables from saved summaries.
mse_table_tex <- function(chosen, methods = c("OF", "UF", "PT", "S", "PS", "AL", "RIDGE")) {
  settings <- unique(chosen[, c("id", "experiment", "n", "coupling", "delta_star")])
  settings <- settings[order(settings$experiment == "small_n", settings$n,
                            settings$coupling, settings$delta_star), ]
  rows <- c("\\begin{tabular}{@{}ccc*{7}{r}@{}}", "\\toprule",
            "$n$ & $c$ & $\\Delta^\\ast$ & OF & UF & PT & S & PS & AL & RIDGE \\\\",
            "\\midrule")
  for (i in seq_len(nrow(settings))) {
    a <- settings[i, ]
    z <- chosen[chosen$id == a$id, ]
    z <- z[match(methods, z$method), ]
    stopifnot(!anyNA(z$method))
    new_n <- i == 1L || a$n != settings$n[i - 1L] ||
      a$experiment != settings$experiment[i - 1L]
    new_c <- new_n || a$coupling != settings$coupling[i - 1L]
    same_n <- settings$n == a$n & settings$experiment == a$experiment
    same_c <- same_n & settings$coupling == a$coupling
    if (i > 1L && new_n) rows <- c(rows, "\\midrule")
    else if (i > 1L && new_c) rows <- c(rows, "\\cmidrule(l){2-10}")
    # Two actual rows per setting permit precise multirow spans.
    # Include the vertical space occupied by internal booktabs rules.
    rules <- length(unique(settings$coupling[same_n])) - 1L
    shift <- if (rules > 0L)
      sprintf("[-\\dimexpr(\\aboverulesep+\\belowrulesep+\\cmidrulewidth)*%d/2\\relax]", rules) else ""
    n_cell <- if (new_n) sprintf("\\multirow[c]{%d}{*}%s{%d}", 2L * sum(same_n), shift, a$n) else ""
    c_cell <- if (new_c) sprintf("\\multirow[c]{%d}{*}{%.1f}", 2L * sum(same_c), a$coupling) else ""
    d_cell <- sprintf("\\multirow[c]{2}{*}{%.2f}", a$delta_star)
    rows <- c(rows,
      paste0(paste(c(n_cell, c_cell, d_cell, sprintf("%.3f", 100 * z$MSE)), collapse = " & "), " \\\\"),
      paste0(paste(c("", "", "", sprintf("{\\scriptsize(%.3f)}", 100 * z$MSE_MCSE)), collapse = " & "), " \\\\"))
  }
  c(rows, "\\bottomrule", "\\end{tabular}")
}

wald_table_tex <- function(cal) {
  # The common replication count appears once in the caption.
  stopifnot(length(unique(cal$valid)) == 1L, all(cal$valid == 1000L))
  cal <- cal[order(cal$n, cal$coupling), ]
  for (n in unique(cal$n)) stopifnot(length(unique(cal$rejection[cal$n == n])) == 1L)
  cal <- cal[!duplicated(cal$n), ]
  rows <- c("\\begin{tabular}{@{}crrr@{}}", "\\toprule",
            "$n$ & Rejection (\\%) & MCSE (pp) & 95\\% MC interval (\\%) \\\\", "\\midrule")
  for (i in seq_len(nrow(cal))) {
    z <- cal[i, ]
    rows <- c(rows, sprintf("%d & %.1f & %.2f & [%.2f, %.2f] \\\\",
                           z$n, 100*z$rejection, 100*z$MCSE, 100*z$MC_lower, 100*z$MC_upper))
  }
  c(rows, "\\bottomrule", "\\end{tabular}")
}

selection_table_tex <- function(sel_range) {
  # Constant strong-variable retention appears once in the caption.
  stopifnot(all(sel_range$all_strong_min == 1))
  sel_range <- sel_range[order(sel_range$n, sel_range$coupling), ]
  rows <- c("\\begin{tabular}{@{}ccrrrr@{}}", "\\toprule",
            "$n$ & $c$ & Exact support & Weak recall & FPR & Fallback (max.) \\\\", "\\midrule")
  for (i in seq_len(nrow(sel_range))) {
    z <- sel_range[i, ]
    new_n <- i == 1L || z$n != sel_range$n[i - 1L]
    if (i > 1L && new_n) rows <- c(rows, "\\midrule")
    n_cell <- if (new_n) sprintf("\\multirow[c]{%d}{*}{%d}", sum(sel_range$n == z$n), z$n) else ""
    rows <- c(rows, sprintf("%s & %.1f & %.3f--%.3f & %.3f--%.3f & %.3f--%.3f & %.3f \\\\",
      n_cell, z$coupling, z$exact_min, z$exact_max, z$weak_min, z$weak_max,
      z$FPR_min, z$FPR_max, z$fallback_max))
  }
  c(rows, "\\bottomrule", "\\end{tabular}")
}

# Settings (as in the paper). Important: changing seeds or replications changes the results.
mode <- "full" # "full" for the paper; "small" for a quick check (R = 3)
R <- if (mode == "full") 500 else 3
Rtheory <- if (mode == "full") 1000 else 3
nlambda <- 16
outdir <- file.path("output", "simulation")
dir.create(file.path(outdir, "cells"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "figures"), showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE)

core_grid <- expand.grid(n = c(250, 500), rho = c(.2, .5, .8), coupling = c(.3, .6),
  p2 = 8, delta_star = c(0, .1, .25, .5, .75, 1), KEEP.OUT.ATTRS = FALSE)
core_grid$experiment <- "core"
small_n <- expand.grid(n = 100, rho = .5, coupling = .9, p2 = 8,
  delta_star = c(0, .75), KEEP.OUT.ATTRS = FALSE)
small_n$experiment <- "small_n"
finite_settings <- rbind(core_grid, small_n)
finite_settings$W <- "queen"
finite_settings$direction <- "sparse_two"
finite_settings$id <- seq_len(nrow(finite_settings))
finite_settings$design_seed <- 710000L + match(finite_settings$n, c(100, 250, 500)) * 10000L
finite_settings$response_seed <- 810000L + match(finite_settings$n, c(100, 250, 500)) * 100000L

theory_settings <- expand.grid(n = c(100, 250, 500), rho = .5, coupling = c(.3, .6), p2 = 8,
  Delta = c(0, 1, 4, 9), KEEP.OUT.ATTRS = FALSE)
theory_settings$id <- seq_len(nrow(theory_settings))
theory_settings$experiment <- "theory"
theory_settings$design_seed <- 710000L + match(theory_settings$n, c(100, 250, 500)) * 10000L
theory_settings$response_seed <- 1810000L + match(theory_settings$n, c(100, 250, 500)) * 100000L

# Experiment 1: 74 settings, R replications each. X is generated once per n and c.
# Run time for mode = "full": about 2.5 hours. Saved cells are reused.
weights_cache <- list(); design_cache <- list()
get_inputs <- function(setting) {
  wk <- as.character(setting$n); dk <- paste(setting$n, setting$p2, setting$coupling, sep = "_")
  if (is.null(weights_cache[[wk]])) weights_cache[[wk]] <<- make_weights(setting$n)
  if (is.null(design_cache[[dk]])) design_cache[[dk]] <<- make_block_design(setting$n, 4,
    setting$p2, setting$coupling, setting$design_seed)
  list(w = weights_cache[[wk]], X = design_cache[[dk]])
}
finite_results <- vector("list", nrow(finite_settings)); maps <- list()
for (i in seq_len(nrow(finite_settings))) {
  setting <- finite_settings[i, ]; inp <- get_inputs(setting)
  file <- file.path(outdir, "cells", sprintf("finite_%03d.rds", i))
  if (file.exists(file)) result <- readRDS(file) else {
    result <- finite_cell(setting, inp$w, inp$X, R, nlambda, setting$response_seed)
    saveRDS(result, file)
  }
  finite_results[[i]] <- result
  if (setting$n == 250 && setting$coupling == .6 && setting$rho %in% c(.2, .8) &&
      setting$delta_star %in% c(0, .75))
    maps[[length(maps) + 1]] <- spatial_summary(result, inp$w, "selected")
  cat(sprintf("Finite %d/%d: valid %d/%d\n", i, nrow(finite_settings),
              sum(result$diagnostics$selected_ok), nrow(result$diagnostics)))
}

# Experiment 2: 24 settings, Rtheory replications each
theory_results <- vector("list", nrow(theory_settings))
for (i in seq_len(nrow(theory_settings))) {
  setting <- theory_settings[i, ]; inp <- get_inputs(setting)
  file <- file.path(outdir, "cells", sprintf("theory_%03d.rds", i))
  if (file.exists(file)) result <- readRDS(file) else {
    result <- theory_cell(setting, inp$w, inp$X, Rtheory, setting$response_seed)
    saveRDS(result, file)
  }
  theory_results[[i]] <- result
  cat(sprintf("Theory %d/%d: valid %d/%d\n", i, nrow(theory_settings),
              sum(result$diagnostics$valid), nrow(result$diagnostics)))
}
finite <- do.call(rbind, lapply(finite_results, function(z) z$raw))
theory <- do.call(rbind, lapply(theory_results, function(z) z$raw))
fd <- do.call(rbind, lapply(finite_results, function(z) z$diagnostics))
td <- do.call(rbind, lapply(theory_results, function(z) z$diagnostics))
expected <- do.call(rbind, lapply(theory_results, function(z) z$expected))
sm <- do.call(rbind, maps)

# Failed fits are counted, not replaced
cat("Experiment 1: requested", nrow(fd), "; failed", sum(!fd$selected_ok), "\n")
cat("Experiment 2: requested", nrow(td), "; failed", sum(!td$valid), "\n")

# Risk summaries: MSE of beta1, relative efficiency RE = MSE(OF) / MSE, Monte Carlo SEs
summarize_risk <- function(raw) {
  groups <- split(raw, interaction(raw$id, raw$block, raw$method, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    baseline <- raw[raw$id == z$id[1] & raw$block == z$block[1] & raw$method == "OF", ]
    of_loss <- baseline$loss[match(z$replication, baseline$replication)]
    errors <- as.matrix(z[, paste0("error", 1:4)])
    MSE <- mean(z$loss); RE <- mean(of_loss) / MSE
    keep <- intersect(c("id", "n", "rho", "coupling", "p2", "delta_star", "Delta",
                        "experiment", "block", "method"), names(z))
    answer <- z[1, keep, drop = FALSE]
    answer$valid <- nrow(z)
    answer$bias2 <- sum(colMeans(errors)^2)
    answer$variance <- mean(rowSums(sweep(errors, 2, colMeans(errors), "-")^2))
    answer$MSE <- MSE; answer$MSE_MCSE <- sd(z$loss) / sqrt(nrow(z))
    answer$gain_percent <- 100 * (1 - MSE / mean(of_loss))
    answer$RE <- RE
    answer$scaled_risk <- z$n[1] * MSE
    answer$scaled_MCSE <- z$n[1] * answer$MSE_MCSE
    for (metric in c("ADE", "ATE", "AIE"))
      answer[[paste0(metric, "_RMSE")]] <- sqrt(mean(z[[paste0(metric, "_loss")]]))
    answer
  }))
}
fs <- summarize_risk(finite); ts <- summarize_risk(theory)
ts <- merge(ts, expected[, c("id", "method", "ADR")], by = c("id", "method"), sort = FALSE)
write.csv(fs, file.path(outdir, "finite_summary.csv"), row.names = FALSE)
write.csv(ts, file.path(outdir, "theory_summary.csv"), row.names = FALSE)

safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
selection <- do.call(rbind, lapply(split(fd, fd$id), function(z) {
  answer <- z[1, c("id", "n", "rho", "coupling", "p2", "delta_star", "experiment")]
  for (name in c("strong_recall", "all_strong_retained", "weak_recall", "FPR", "exact_support"))
    answer[[name]] <- safe_mean(z[[name]])
  answer$q_below_three <- safe_mean(z$q < 3)
  answer$q_at_most_four <- safe_mean(z$q <= 4)
  answer$S_PS_fallback <- safe_mean(z$fallback)
  answer
}))
# Wald rejection rate at Delta = 0, with Wilson Monte Carlo intervals
calibration <- do.call(rbind, lapply(split(td, td$id), function(z) {
  answer <- z[1, c("id", "n", "coupling", "p2", "Delta")]
  answer$valid <- sum(z$valid); answer$rejection <- safe_mean(z$rejected)
  answer$MCSE <- sqrt(answer$rejection * (1 - answer$rejection) / answer$valid)
  h <- qnorm(.975); m <- answer$valid; r <- answer$rejection
  centre <- (r + h^2 / (2 * m)) / (1 + h^2 / m)
  half <- h * sqrt(r * (1 - r) / m + h^2 / (4 * m^2)) / (1 + h^2 / m)
  answer$MC_lower <- centre - half; answer$MC_upper <- centre + half
  answer
}))

# Tables: Table 1 (coefficient MSE), Table 2 (Wald), appendix selection table
methods <- c("OF", "UF", "PT", "S", "PS", "AL", "RIDGE")
chosen <- fs[fs$block == "selected" & fs$rho == .5 &
  ((fs$experiment == "core" & fs$delta_star %in% c(0, .25, .75)) | fs$experiment == "small_n"), ]
writeLines(mse_table_tex(chosen, methods), file.path(outdir, "tables", "table1_coefficient_mse.tex"))
null <- calibration[calibration$Delta == 0, ]
writeLines(wald_table_tex(null), file.path(outdir, "tables", "table2_wald.tex"))
sel <- selection[selection$experiment == "core", ]
sel_range <- do.call(rbind, lapply(split(sel, interaction(sel$n, sel$coupling, drop = TRUE)), function(z)
  data.frame(n = z$n[1], coupling = z$coupling[1], all_strong_min = min(z$all_strong_retained),
    exact_min = min(z$exact_support), exact_max = max(z$exact_support),
    weak_min = min(z$weak_recall, na.rm = TRUE), weak_max = max(z$weak_recall, na.rm = TRUE),
    FPR_min = min(z$FPR), FPR_max = max(z$FPR), fallback_max = max(z$S_PS_fallback))))
writeLines(selection_table_tex(sel_range), file.path(outdir, "tables", "selection.tex"))
print(null[, c("n", "coupling", "valid", "rejection", "MC_lower", "MC_upper")], row.names = FALSE)

# Figures
colours <- c(OF = "#666666", UF = "#0072B2", PT = "#E69F00", S = "#CC79A7", PS = "#009E73",
             AL = "#D55E00", RIDGE = "#56B4E9")
types <- c(OF = "solid", UF = "longdash", PT = "dashed", S = "dotdash", PS = "solid",
           AL = "twodash", RIDGE = "dotted")
theme_set(theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "#F3F4F5", colour = "#D5D8DB"),
  legend.position = "bottom", legend.title = element_blank(), legend.key.width = grid::unit(1.1, "cm")))
colour_scale <- scale_colour_manual(values = colours, breaks = methods, drop = FALSE)
line_scale <- scale_linetype_manual(values = types, breaks = methods, drop = FALSE)
save_figure <- function(p, name, width, height)
  ggsave(file.path(outdir, "figures", paste0(name, ".png")), p, width = width, height = height,
         dpi = 320, bg = "white")
facet_nc <- facet_grid(coupling ~ n, labeller = label_bquote(rows = c == .(coupling), cols = n == .(n)))

# Figure 1: relative efficiency against Delta*, rho = 0.5
main <- fs[fs$experiment == "core" & fs$block == "selected" & fs$rho == .5 & fs$method != "OF", ]
main$method <- factor(main$method, levels = methods)
p1 <- ggplot(main, aes(delta_star, RE, colour = method, linetype = method)) +
  geom_hline(yintercept = 1, colour = "#777777") + geom_line(linewidth = .75) +
  geom_point(size = 1.8, show.legend = FALSE) + facet_nc +
  labs(x = expression(Delta^"*"), y = "Relative efficiency") +
  scale_colour_manual(values = colours, breaks = methods[-1]) +
  scale_linetype_manual(values = types, breaks = methods[-1])
save_figure(p1, "fig1_relative_efficiency", 9.5, 6.4)

# Figure 2: empirical scaled risk (points, 95% MC intervals) and ADR curves (Theorem 4)
curves <- list()
for (i in which(theory_settings$Delta == 0)) {
  cell <- theory_results[[i]]; nu <- cell$direction
  base <- adr_formula(cell$C0, rep(0, 8), 4)
  for (Delta in seq(0, 9, length.out = 80)) {
    kappa <- sqrt(Delta / as.numeric(crossprod(nu, base$C22 %*% nu))) * nu
    a <- adr_formula(cell$C0, kappa, 4)
    curves[[length(curves) + 1]] <- data.frame(n = theory_settings$n[i],
      coupling = theory_settings$coupling[i], Delta = Delta, method = names(a$risk), ADR = a$risk)
  }
}
curves <- do.call(rbind, curves)
p2 <- ggplot(curves, aes(Delta, ADR, colour = method, linetype = method)) + geom_line(linewidth = .75) +
  geom_errorbar(data = ts, aes(y = scaled_risk, ymin = pmax(0, scaled_risk - 1.96 * scaled_MCSE),
    ymax = scaled_risk + 1.96 * scaled_MCSE), width = .1, linewidth = .35, show.legend = FALSE) +
  geom_point(data = ts, aes(y = scaled_risk), size = 1.5, show.legend = FALSE) +
  facet_grid(coupling ~ n, scales = "free_y", labeller = label_bquote(rows = c == .(coupling), cols = n == .(n))) +
  scale_x_continuous(breaks = c(0, 1, 4, 9)) + labs(x = expression(Delta), y = "Scaled risk") +
  colour_scale + line_scale
save_figure(p2, "fig2_theory_agreement", 10.5, 6.8)

# Figure 3: coefficient, direct-effect and indirect-effect RMSE (n = 250, c = 0.6), rho_OF for all
sp <- fs[fs$experiment == "core" & fs$block == "selected" & fs$n == 250 & fs$coupling == .6 &
         fs$delta_star %in% c(0, .75), ]
spatial <- do.call(rbind, lapply(c("MSE", "ADE_RMSE", "AIE_RMSE"), function(metric) {
  z <- sp
  z$metric <- switch(metric, MSE = "Coefficient RMSE", ADE_RMSE = "Direct-effect RMSE",
                     AIE_RMSE = "Indirect-effect RMSE")
  z$error <- if (metric == "MSE") sqrt(z[[metric]]) else z[[metric]]
  z
}))
spatial$method <- factor(spatial$method, levels = methods)
p3 <- ggplot(spatial, aes(rho, error, colour = method, linetype = method)) +
  geom_line(linewidth = .7) + geom_point(size = 1.8, show.legend = FALSE) +
  facet_grid(metric ~ delta_star, scales = "free_y", labeller = label_bquote(cols = Delta^"*" == .(delta_star))) +
  scale_x_continuous(breaks = c(.2, .5, .8)) + labs(x = expression(rho), y = NULL) + colour_scale + line_scale
save_figure(p3, "fig3_spatial_effects", 9.6, 8)

# Figure 4: sitewise MSE reduction relative to OF for UF, S and PS
sm$cell <- factor(paste0("rho = ", sm$rho, ", d* = ", sm$delta_star))
gain_map <- function(z, file, height) {
  lim <- 10 * ceiling(max(abs(z$gain_percent), na.rm = TRUE) / 10)
  p <- ggplot(z, aes(x, y, fill = gain_percent)) + geom_tile(colour = "white", linewidth = .2) +
    facet_grid(method ~ cell) + coord_equal() +
    scale_fill_gradient2(low = "#B2182B", mid = "#F0F0F0", high = "#2166AC", midpoint = 0,
                         limits = c(-lim, lim), name = "Gain over OF (%)") + labs(x = NULL, y = NULL) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  save_figure(p, file, 11, height)
}
focus <- sm[sm$method %in% c("UF", "S", "PS"), ]
focus$method <- factor(focus$method, levels = c("UF", "S", "PS"))
gain_map(focus, "fig4_spatial_gain", 4.8)
all_maps <- sm[sm$method != "OF", ]; all_maps$method <- factor(all_maps$method, levels = methods[-1])
gain_map(all_maps, "figS_spatial_gain_all", 7.4)

# Appendix figures: selection rates and Wald QQ at Delta = 0
rates <- do.call(rbind, lapply(c("strong_recall", "exact_support", "q_below_three"), function(metric) {
  z <- selection[selection$experiment == "core" & selection$rho == .5, ]
  z$metric <- switch(metric, strong_recall = "Strong recall", exact_support = "Exact support",
                     q_below_three = "q < 3")
  z$rate <- z[[metric]]; z
}))
p5 <- ggplot(rates, aes(delta_star, rate, colour = metric, linetype = metric)) +
  geom_line(linewidth = .75) + geom_point(size = 1.7) + facet_nc +
  scale_y_continuous(limits = c(0, 1)) + labs(x = expression(Delta^"*"), y = "Proportion")
save_figure(p5, "figS_selection", 9.5, 6.4)
qq <- do.call(rbind, lapply(split(td[td$Delta == 0 & td$valid, ], td$id[td$Delta == 0 & td$valid]), function(z) {
  z$reference <- qchisq(ppoints(nrow(z)), z$p2[1]); z$empirical <- sort(z$T); z
}))
p6 <- ggplot(qq, aes(reference, empirical)) + geom_abline(slope = 1, intercept = 0, colour = "#777777") +
  geom_point(colour = "#0072B2", size = .9, alpha = .65) + facet_nc +
  labs(x = expression(chi[8]^2 ~ "quantile"), y = "Wald quantile")
save_figure(p6, "figS_wald_qq", 11.5, 6.5)
cat("Outputs saved in", outdir, "\n")
