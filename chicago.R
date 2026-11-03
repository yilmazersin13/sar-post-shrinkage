# Chicago income application (Section 6.1): full-sample fit, geographic validation, Table 3,

library(Matrix); library(ggplot2); library(sf)
source("sar_functions.R")
source("real_data_functions.R")

# rule = "original": PT/S/PS use their own coefficients with rho_OF
# rule = "fixed": prediction rule shrink_theta()
rule <- "original"

use_basemap <- TRUE
tile_alpha <- 0.5
outdir <- file.path("output", "chicago")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
methods <- c("OF", "UF", "PT", "S", "PS", "AL", "RIDGE")

fit_all <- function(md) {
  full <- marginal_qmle(md)
  al <- marginal_tune(md, full, "AL"); ridge <- marginal_tune(md, full, "RIDGE")
  restricted <- marginal_qmle(md, al$retained)
  info <- marginal_information(md, full$beta, full$rho, full$sigma2)
  post <- post_shrinkage(full, restricted, al$retained, info$C, md$n)
  list(full = full, restricted = restricted, AL = al, RIDGE = ridge, post = post,
       retained = al$retained)
}

# Tract maps, with an OpenStreetMap basemap 
map_tracts <- function(geo, values, fill_scale, ncol, file, width, height) {
  df <- do.call(rbind, lapply(names(values), function(nm) {
    g <- geo[, "geometry"]; g$panel <- nm; g$value <- values[[nm]]; g
  }))
  df$panel <- factor(df$panel, levels = names(values))
  base <- function(tiles) {
    g <- ggplot(df)
    if (tiles) g <- g + ggspatial::annotation_map_tile(type = "osm", zoomin = 0,
                                                       progress = "none", alpha = tile_alpha)
    g + geom_sf(aes(fill = value), colour = "white", linewidth = 0.05, alpha = if (tiles) 0.85 else 1) +
      fill_scale + facet_wrap(~panel, ncol = ncol) + theme_bw(base_size = 10) +
      theme(panel.grid = element_blank(), axis.text = element_text(size = 7),
            strip.background = element_rect(fill = "grey96", colour = "grey80"),
            strip.text = element_text(face = "bold"), legend.position = "bottom",
            legend.key.width = grid::unit(1.2, "cm"))
  }
  tiles <- use_basemap && requireNamespace("ggspatial", quietly = TRUE)
  for (ext in c("pdf", "png")) {
    out <- file.path(outdir, paste0(file, ".", ext))
    res <- if (tiles) try(ggsave(out, base(TRUE), width = width, height = height, dpi = 400,
                                 bg = "white"), silent = TRUE) else NULL
    if (!tiles || inherits(res, "try-error")) {
      if (tiles) message("Basemap failed: ", conditionMessage(attr(res, "condition")))
      tiles <- FALSE
      ggsave(out, base(FALSE), width = width, height = height, dpi = 400, bg = "white")
    }
  }
}
diverging <- function(values, name) {
  lim <- max(abs(unlist(values)))
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
                       limits = c(-lim, lim), name = name)
}

# Data: 791 census tracts, log per-capita income, queen contiguity
s <- st_transform(st_read(file.path("data", "Chi-SDOH.geojson"), quiet = TRUE), 26971)
stopifnot(all(st_is_valid(s)), !any(st_is_empty(s)))
d <- st_drop_geometry(s)
ids <- sprintf("%.0f", d$geoid10)
pct <- c("EP_NOHSDP", "Unemp14", "EP_AGE17", "Ovr6514P", "EP_DISABL", "EP_SNGPNT", "EP_CROWD",
         "EP_LIMENG", "EP_NOVEH", "EP_MUNIT", "EP_GROUPQ", "CTA14P")
fields <- c("EP_PCI", "Pop2014", pct, "CmTm14", "PDENS14", "MEANMI_14")
# no row is removed; stop if a value is invalid
stopifnot(!anyDuplicated(ids), complete.cases(d[, fields]), all(d$EP_PCI > 0), all(d$PDENS14 > 0),
          all(as.matrix(d[, pct]) >= 0 & as.matrix(d[, pct]) <= 100))
y <- log(d$EP_PCI); names(y) <- ids; n <- length(y)
w <- queen_weights(s, ids)
stopifnot(n == 791L)

# Full-sample fit (all responses observed)
X <- chicago_design(d)$X
fs <- fit_all(marginal_data(y, X, w, seq_len(n)))
res <- as.vector(y - fs$full$rho * (w$W %*% y) - X %*% fs$full$beta)
z <- (res - mean(res)) / sqrt(mean((res - mean(res))^2))
full_sample <- c(n = n, q = fs$post$q, T = fs$post$T, rho_OF = fs$full$rho,
  skewness = mean(z^3), excess_kurtosis = mean(z^4) - 3,
  Moran_I = n / sum(w$W) * sum(res * (w$W %*% res)) / sum(res^2))
print(round(full_sample, 4))
cat("AL excluded:", colnames(X)[-fs$retained], "\n")

# Geographic validation: five south-to-north bands of tract centroids, GEOID breaks ties
northing <- st_coordinates(st_centroid(st_geometry(s)))[, 2]
fold <- integer(n); fold[order(northing, ids)] <- as.integer(floor((seq_len(n) - 1) * 5 / n) + 1L)
fold_file <- file.path(outdir, "fold_fits.rds")
fits <- if (file.exists(fold_file)) readRDS(fold_file) else list()
pred <- list(); info_folds <- list()
for (k in 1:5) {
  train <- which(fold != k); held <- which(fold == k)
  # Important: centring, scaling, selection, tuning and refitting use training data only
  md <- marginal_data(y[train], chicago_design(d, train = train)$X, w, train)
  if (length(fits) < k) { fits[[k]] <- fit_all(md); saveRDS(fits, fold_file) }
  f <- fits[[k]]
  est <- c(f$post$estimates, list(AL = f$AL$beta, RIDGE = f$RIDGE$beta))
  rho <- c(OF = f$full$rho, UF = f$restricted$rho, PT = f$full$rho, S = f$full$rho,
           PS = f$full$rho, AL = f$AL$rho, RIDGE = f$RIDGE$rho)
  weight <- shrinkage_weights(f$post)
  if (rule == "fixed") for (m in c("PT", "S", "PS")) {
    th <- shrink_theta(f$full, f$restricted, f$retained, weight[[m]])
    est[[m]] <- th$beta; rho[[m]] <- th$rho
  }
  pred[[k]] <- data.frame(geoid = ids[held], fold = k, observed = y[held],
    sapply(methods, function(m) marginal_predict(md, est[[m]], rho[[m]])))
  info_folds[[k]] <- data.frame(fold = k, q = f$post$q, T = f$post$T, S_weight = weight[["S"]],
                                rho_OF = f$full$rho, rho_UF = f$restricted$rho)
  cat(sprintf("Fold %d done.\n", k))
}
info_folds <- do.call(rbind, info_folds)
p <- do.call(rbind, pred); p <- p[match(ids, p$geoid), ]

# Pooled RMSE and MAE on log income
ref <- mean((p$OF - p$observed)^2)
metrics <- do.call(rbind, lapply(methods, function(m) {
  e <- p[[m]] - p$observed
  data.frame(rule = rule, method = m, RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
             MSE_change_vs_OF = 100 * (mean(e^2) / ref - 1))
}))
print(info_folds, digits = 4, row.names = FALSE)
print(metrics, digits = 4, row.names = FALSE)
write.csv(metrics, file.path(outdir, paste0("metrics_", rule, ".csv")), row.names = FALSE)

# Table 3
r4 <- round(metrics$RMSE, 4); m4 <- round(metrics$MAE, 4)
cell <- function(x, best) paste0(sprintf("%.4f", x), ifelse(x == min(best), "$^{*}$", ""))
writeLines(c("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrr@{}}", "\\toprule",
  "Estimator & RMSE & MAE \\\\", "\\midrule",
  paste0(methods, " & ", cell(r4, r4), " & ", cell(m4, m4), " \\\\"),
  "\\bottomrule", "\\end{tabular*}"), file.path(outdir, paste0("table3_chicago_", rule, ".tex")))

# Figure 5
geo <- st_transform(s, 4326)
map_tracts(geo, list("Per-capita income" = d$EP_PCI / 1000),
           scale_fill_viridis_c(option = "C", name = "Income ($000)"), 1, "income_map", 6, 6.5)
qq <- ggplot(data.frame(normal = qnorm(ppoints(n)), residual = sort(z)), aes(normal, residual)) +
  geom_point(size = 1, colour = "#28648c", alpha = .65) +
  geom_abline(intercept = 0, slope = 1, colour = "#b24b43") +
  labs(title = "Residual QQ", x = "Normal quantile", y = "Standardized residual") + theme_bw(base_size = 11)
for (ext in c("pdf", "png"))
  ggsave(file.path(outdir, paste0("residual_qq.", ext)), qq, width = 6, height = 4.8, dpi = 300, bg = "white")

# Validation maps for PS
e <- list(OF = p$observed - p$OF, UF = p$observed - p$UF, PS = p$observed - p$PS)
map_tracts(geo, list(Observed = p$observed, OF = p$OF, UF = p$UF, PS = p$PS),
           scale_fill_viridis_c(option = "C", name = "Log income"), 4, "ps_income", 12, 5.5)
map_tracts(geo, e, diverging(e, "Observed - predicted"), 3, "ps_errors", 10, 5.5)
g2 <- list("PS vs OF" = e$OF^2 - e$PS^2, "PS vs UF" = e$UF^2 - e$PS^2)
map_tracts(geo, g2, diverging(g2, "Squared-error reduction"), 2, "ps_gains", 7.5, 5.8)
cat("Outputs saved in", outdir, "\n")
