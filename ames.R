# Ames housing application (Section 6.2): full-sample fit, geographic validation, Table 4,

library(Matrix); library(ggplot2)
source("sar_functions.R")
source("real_data_functions.R")


use_basemap <- TRUE
tile_alpha <- 0.5
jitter_m <- 15
outdir <- file.path("output", "ames")
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


map_points <- function(df, fill_scale, ncol, file, width, height) {
  base <- function(tiles) {
    if (tiles) {
      g <- ggplot(sf::st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)) +
        ggspatial::annotation_map_tile(type = "osm", zoomin = 0, progress = "none", alpha = tile_alpha) +
        geom_sf(aes(fill = value), shape = 21, colour = "grey20", stroke = 0.15, size = 2.4) +
        ggspatial::annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.6)
    } else {
      g <- ggplot(df, aes(longitude, latitude, fill = value)) +
        geom_point(shape = 21, colour = "grey20", stroke = 0.15, size = 2.4) +
        coord_fixed(ratio = 1 / cos(mean(df$latitude) * pi / 180))
    }
    g + fill_scale + facet_wrap(~panel, ncol = ncol) + labs(x = NULL, y = NULL) +
      theme_bw(base_size = 10) +
      theme(panel.grid = element_blank(), axis.text = element_text(size = 7),
            strip.background = element_rect(fill = "grey96", colour = "grey80"),
            strip.text = element_text(face = "bold"), legend.position = "bottom",
            legend.key.width = grid::unit(1.2, "cm"))
  }
  tiles <- use_basemap && requireNamespace("sf", quietly = TRUE) &&
    requireNamespace("ggspatial", quietly = TRUE)
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
panel_data <- function(values, lon, lat, by_size = FALSE) {
  out <- do.call(rbind, lapply(names(values), function(nm)
    data.frame(longitude = lon, latitude = lat, panel = nm, value = values[[nm]])))
  out$panel <- factor(out$panel, levels = names(values))
  # Large absolute values are drawn last
  if (by_size) out <- out[order(abs(out$value)), ]
  out
}
diverging <- function(x, name) {
  lim <- max(abs(x$value))
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
                       limits = c(-lim, lim), name = name)
}

# Data: all 648 sales from 2009, log sale price
load(file.path("data", "ames.rda"))
stopifnot(unname(tools::md5sum(file.path("data", "ames.rda"))) == "d8728a5f02457b4c0da32c0573559dc9")
d <- as.data.frame(ames); d$source_row <- seq_len(nrow(d))
d <- d[as.character(d$Year_Sold) == "2009", ]
check_ames_fields(d)
y <- log(d$Sale_Price); names(y) <- d$source_row; n <- length(y)
xy <- as.matrix(d[, c("Longitude", "Latitude")])
w <- make_observed_weights(xy, 4L); dimnames(w$W) <- list(names(y), names(y))
price <- d$Sale_Price
d$Sale_Price <- NULL #  the design never sees the response
stopifnot(n == 648L)

# Full-sample fit (all responses observed)
design <- ames_broad_features(d); X <- design$X
fs <- fit_all(marginal_data(y, X, w, seq_len(n)))
res <- as.vector(y - fs$full$rho * (w$W %*% y) - X %*% fs$full$beta)
z <- (res - mean(res)) / sqrt(mean((res - mean(res))^2))
q <- fs$post$q; Tn <- fs$post$T
full_sample <- c(n = n, columns = ncol(X) - 1, AL_columns = length(fs$retained) - 1, q = q,
  T = Tn, cutoff = qchisq(.95, q), S_weight = (q - 2) / Tn, rho_OF = fs$full$rho,
  skewness = mean(z^3), excess_kurtosis = mean(z^4) - 3,
  max_leverage = max(rowSums(qr.Q(qr(X))^2)),
  Moran_I = n / sum(w$W) * sum(res * (w$W %*% res)) / sum(res^2))
print(round(full_sample, 4))

# Geographic validation: five latitude bands; equal latitudes stay together
fold <- findInterval(d$Latitude, quantile(d$Latitude, (1:4) / 5, type = 1, names = FALSE)) + 1L
fold_file <- file.path(outdir, "fold_fits.rds")
fits <- if (file.exists(fold_file)) readRDS(fold_file) else list()
pred <- list(original = list(), fixed = list()); info_folds <- list()
for (k in 1:5) {
  train <- which(fold != k); held <- which(fold == k)
  trained <- ames_broad_features(d[train, ])
  applied <- apply_ames_validation_recipe(d, trained$recipe)
  md <- marginal_data(y[train], applied$design$X, w, train)
  if (length(fits) < k) { fits[[k]] <- fit_all(md); saveRDS(fits, fold_file) }
  f <- fits[[k]]
  est <- c(f$post$estimates, list(AL = f$AL$beta, RIDGE = f$RIDGE$beta))
  rho <- c(OF = f$full$rho, UF = f$restricted$rho, PT = f$full$rho, S = f$full$rho,
           PS = f$full$rho, AL = f$AL$rho, RIDGE = f$RIDGE$rho)
  weight <- shrinkage_weights(f$post)
  fixed <- est; rho_fixed <- rho
  for (m in c("PT", "S", "PS")) {
    th <- shrink_theta(f$full, f$restricted, f$retained, weight[[m]])
    fixed[[m]] <- th$beta; rho_fixed[[m]] <- th$rho
  }
  id <- data.frame(source_row = d$source_row[held], fold = k,
                   longitude = xy[held, 1], latitude = xy[held, 2])
  pred$original[[k]] <- data.frame(id, observed = y[held],
    sapply(methods, function(m) marginal_predict(md, est[[m]], rho[[m]])))
  pred$fixed[[k]] <- data.frame(id, observed = y[held],
    sapply(methods, function(m) marginal_predict(md, fixed[[m]], rho_fixed[[m]])))
  info_folds[[k]] <- data.frame(fold = k, columns = ncol(md$X) - 1, q = f$post$q, T = f$post$T,
    S_weight = weight[["S"]], rho_OF = f$full$rho, rho_UF = f$restricted$rho,
    unseen_homes = length(unique(applied$unseen$source_row)))
  cat(sprintf("Fold %d done.\n", k))
}
info_folds <- do.call(rbind, info_folds)

# Pooled RMSE and MAE on log price
metrics <- do.call(rbind, lapply(names(pred), function(rule) {
  p <- do.call(rbind, pred[[rule]]); ref <- mean((p$OF - p$observed)^2)
  do.call(rbind, lapply(methods, function(m) {
    e <- p[[m]] - p$observed
    data.frame(rule = rule, method = m, RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
               MSE_change_vs_OF = 100 * (mean(e^2) / ref - 1))
  }))
}))
print(info_folds, digits = 4, row.names = FALSE)
print(metrics, digits = 4, row.names = FALSE)
write.csv(metrics, file.path(outdir, "metrics.csv"), row.names = FALSE)

# Table 4: fixed prediction rule; * marks the smallest value in each column
tab <- metrics[metrics$rule == "fixed", ]
r4 <- round(tab$RMSE, 4); m4 <- round(tab$MAE, 4)
cell <- function(x, best) paste0(sprintf("%.4f", x), ifelse(x == min(best), "$^{*}$", ""))
writeLines(c("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrr@{}}", "\\toprule",
  "Estimator & RMSE & MAE \\\\", "\\midrule",
  paste0(methods, " & ", cell(r4, r4), " & ", cell(m4, m4), " \\\\"),
  "\\bottomrule", "\\end{tabular*}"), file.path(outdir, "table4_ames.tex"))

# Figure 6: sale prices and OF residual quantiles
set.seed(2009) # display jitter only
lat0 <- mean(xy[, 2]) * pi / 180
lon_j <- xy[, 1] + runif(n, -1, 1) * jitter_m / (111320 * cos(lat0))
lat_j <- xy[, 2] + runif(n, -1, 1) * jitter_m / 110540
map_points(panel_data(list("Sale price" = price / 1000), lon_j, lat_j),
           scale_fill_viridis_c(option = "C", name = "Price ($000)"), 1, "price_map", 6, 5.2)
qq <- ggplot(data.frame(normal = qnorm(ppoints(n)), residual = sort(z)), aes(normal, residual)) +
  geom_point(size = 1, colour = "#28648c", alpha = .65) +
  geom_abline(intercept = 0, slope = 1, colour = "#b24b43") +
  labs(title = "Residual QQ", x = "Normal quantile", y = "Standardized residual") + theme_bw(base_size = 11)
for (ext in c("pdf", "png"))
  ggsave(file.path(outdir, paste0("residual_qq.", ext)), qq, width = 6, height = 4.8, dpi = 300, bg = "white")

# Validation maps for PS under the fixed rule
p <- do.call(rbind, pred$fixed)
p <- p[match(d$source_row, p$source_row), ] # same row order as lon_j, lat_j
e <- data.frame(OF = p$observed - p$OF, UF = p$observed - p$UF, PS = p$observed - p$PS)
map_points(panel_data(list(Observed = p$observed, OF = p$OF, UF = p$UF, PS = p$PS), lon_j, lat_j),
           scale_fill_viridis_c(option = "C", name = "Log price"), 2, "ps_prices", 8, 8.5)
err <- panel_data(as.list(e), lon_j, lat_j, by_size = TRUE)
map_points(err, diverging(err, "Observed - predicted"), 3, "ps_errors", 11, 4.8)
gain <- panel_data(list("PS vs OF" = e$OF^2 - e$PS^2, "PS vs UF" = e$UF^2 - e$PS^2), lon_j, lat_j, TRUE)
map_points(gain, diverging(gain, "Squared-error reduction"), 2, "ps_gains", 8, 4.8)
gain_abs <- panel_data(list("PS vs OF" = abs(e$OF) - abs(e$PS), "PS vs UF" = abs(e$UF) - abs(e$PS)),
                       lon_j, lat_j, TRUE)
map_points(gain_abs, diverging(gain_abs, "Absolute-error reduction"), 2, "ps_gains_abs", 8, 4.8)
cat("Outputs saved in", outdir, "\n")
