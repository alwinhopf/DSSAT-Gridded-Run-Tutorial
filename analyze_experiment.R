# =============================================================================
# analyze_experiment.R — post-hoc analysis of a weather x soil sweep.
#
# Turns the combined.csv produced by run_experiment.R into the figures and
# tables that make input differences legible:
#
#   * summary_by_combo.csv          mean / sd / CV / range per input combination
#   * fig_yield_maps.png            small-multiple maps, one panel per combination
#   * fig_sensitivity_cv_map.png    per-point CV across combinations (where inputs matter most)
#   * fig_variance_decomposition.png  ANOVA share: weather vs soil vs period vs residual
#   * pairwise_rmsd.csv / .png       RMSD between every pair of combinations (how interchangeable)
#   * rank_stability.csv / .png      Spearman corr. of the spatial yield pattern between combinations
#
# USAGE
#   Rscript analyze_experiment.R results/EXPERIMENT_<name>_combined.csv
#   Rscript analyze_experiment.R <combined.csv> final_grain_kg_ha --treatment 1 --out analysis/run1
#
# The response variable defaults to final_grain_kg_ha; pass any numeric column as
# the 2nd argument. By default points are averaged over treatments and years
# (per combination); restrict with --treatment N.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
flag <- function(name, default = NA) { i <- which(args == name); if (length(i) && i < length(args)) args[i + 1] else default }
pos  <- args[!grepl("^--", args) & !(args %in% c(flag("--treatment"), flag("--out")))]
CSV  <- pos[1]
RV   <- if (length(pos) >= 2) pos[2] else "final_grain_kg_ha"
TRT  <- suppressWarnings(as.integer(flag("--treatment")))
if (is.na(CSV) || !file.exists(CSV)) stop("Usage: Rscript analyze_experiment.R <combined.csv> [response_var] [--treatment N] [--out DIR]")
OUT  <- flag("--out", file.path("analysis", sub("_combined\\.csv$", "", basename(CSV))))
# Optional boundary outline drawn under the maps (state/country borders). Defaults
# to the bundled US states; pass --boundary <file.shp> for another region, or
# --boundary none to disable.
BOUNDARY <- flag("--boundary", file.path("shapefile", "tl_2024_us_state.shp"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

have_gg <- requireNamespace("ggplot2", quietly = TRUE)
if (have_gg) suppressMessages(library(ggplot2))

# Build a boundary underlay (geom_sf) + matching coord_sf cropped to the data,
# or NULL to fall back to coord_quickmap(). Cheap & optional: needs sf + the file.
load_boundary <- function(df) {
  if (!have_gg || is.na(BOUNDARY) || identical(tolower(BOUNDARY), "none") ||
      !file.exists(BOUNDARY) || !requireNamespace("sf", quietly = TRUE)) return(NULL)
  b <- tryCatch(sf::st_transform(sf::read_sf(BOUNDARY), 4326), error = function(e) NULL)
  if (is.null(b)) return(NULL)
  mx <- 0.5
  xl <- c(min(df$longitude) - mx, max(df$longitude) + mx)
  yl <- c(min(df$latitude)  - mx, max(df$latitude)  + mx)
  bb <- sf::st_bbox(c(xmin = xl[1], ymin = yl[1], xmax = xl[2], ymax = yl[2]), crs = 4326)
  b  <- tryCatch(suppressWarnings(sf::st_crop(b, bb)), error = function(e) b)
  list(layer = geom_sf(data = b, fill = NA, colour = "grey55", linewidth = 0.2, inherit.aes = FALSE),
       coord = coord_sf(xlim = xl, ylim = yl, expand = FALSE))
}

d <- read.csv(CSV, stringsAsFactors = FALSE)
need <- c("latitude","longitude","weather_source","soil_source", RV)
if (!all(need %in% names(d))) stop("combined.csv missing columns: ", paste(setdiff(need, names(d)), collapse = ", "))
d[[RV]] <- suppressWarnings(as.numeric(d[[RV]]))
d <- d[is.finite(d[[RV]]), ]
if (!is.null(TRT) && !is.na(TRT) && "treatment" %in% names(d)) d <- d[d$treatment == TRT, ]
if ("point_id" %in% names(d)) d$point_id <- as.character(d$point_id) else
  d$point_id <- as.character(interaction(round(d$latitude,5), round(d$longitude,5), drop = TRUE))

# Combination label: weather | soil (| period when more than one period present)
multi_period <- "period" %in% names(d) && length(unique(d$period)) > 1
d$combo <- if (multi_period) paste(d$weather_source, d$soil_source, d$period, sep = " | ") else
                              paste(d$weather_source, d$soil_source, sep = " | ")
cat(sprintf("Loaded %d rows | response = %s | %d combinations | %d points%s\n",
            nrow(d), RV, length(unique(d$combo)), length(unique(d$point_id)),
            if (!is.na(TRT)) sprintf(" | treatment %d only", TRT) else ""))

# Per-point, per-combination mean of the response (averaged over treatments/years)
agg <- aggregate(stats::setNames(list(d[[RV]]), "y"),
                 by = list(point_id = d$point_id, combo = d$combo,
                           weather_source = d$weather_source, soil_source = d$soil_source,
                           latitude = d$latitude, longitude = d$longitude), FUN = mean)
if (multi_period) agg$period <- d$period[match(paste(agg$point_id, agg$combo), paste(d$point_id, d$combo))]

# ---- 1. Summary per combination -------------------------------------------
cv <- function(x) if (length(x) && mean(x) != 0) 100 * sd(x) / mean(x) else NA
summ <- do.call(rbind, lapply(split(agg, agg$combo), function(g) data.frame(
  combo = g$combo[1], weather_source = g$weather_source[1], soil_source = g$soil_source[1],
  n_points = nrow(g), mean = mean(g$y), sd = sd(g$y), cv_pct = cv(g$y),
  min = min(g$y), max = max(g$y), stringsAsFactors = FALSE)))
summ <- summ[order(-summ$mean), ]; rownames(summ) <- NULL
write.csv(summ, file.path(OUT, "summary_by_combo.csv"), row.names = FALSE)
cat("summary_by_combo.csv\n"); print(summ, digits = 5)

# ---- Wide matrix: points x combinations (for pairwise comparisons) ---------
combos <- sort(unique(agg$combo))
pts    <- sort(unique(agg$point_id))
W <- matrix(NA_real_, length(pts), length(combos), dimnames = list(pts, combos))
W[cbind(match(agg$point_id, pts), match(agg$combo, combos))] <- agg$y
keep <- rowSums(is.na(W)) == 0           # points present in every combination
Wc <- W[keep, , drop = FALSE]

# ---- 2. Per-point sensitivity (CV across combinations) --------------------
sens <- data.frame(point_id = rownames(Wc),
                   latitude  = agg$latitude[match(rownames(Wc), agg$point_id)],
                   longitude = agg$longitude[match(rownames(Wc), agg$point_id)],
                   mean = apply(Wc, 1, mean), cv_pct = apply(Wc, 1, cv),
                   range = apply(Wc, 1, function(x) max(x) - min(x)), stringsAsFactors = FALSE)
write.csv(sens, file.path(OUT, "sensitivity_by_point.csv"), row.names = FALSE)
cat(sprintf("sensitivity_by_point.csv  (median across-input CV = %.1f%%)\n", median(sens$cv_pct, na.rm = TRUE)))

# ---- 3. Pairwise RMSD + Spearman rank stability between combinations -------
nC <- length(combos); rmsd <- cor_s <- matrix(NA_real_, nC, nC, dimnames = list(combos, combos))
for (i in seq_len(nC)) for (j in seq_len(nC)) {
  rmsd[i, j]  <- sqrt(mean((Wc[, i] - Wc[, j])^2))
  cor_s[i, j] <- if (nrow(Wc) > 2) suppressWarnings(cor(Wc[, i], Wc[, j], method = "spearman")) else NA
}
write.csv(data.frame(combo = combos, rmsd, check.names = FALSE), file.path(OUT, "pairwise_rmsd.csv"), row.names = FALSE)
write.csv(data.frame(combo = combos, cor_s, check.names = FALSE), file.path(OUT, "rank_stability.csv"), row.names = FALSE)
cat("pairwise_rmsd.csv, rank_stability.csv\n")

# ---- 4. Variance decomposition (ANOVA) ------------------------------------
terms <- c()
if (length(unique(agg$weather_source)) > 1) terms <- c(terms, "weather_source")
if (length(unique(agg$soil_source))    > 1) terms <- c(terms, "soil_source")
if (multi_period && length(unique(agg$period)) > 1) terms <- c(terms, "period")
vd <- NULL
if (length(terms) >= 1 && nrow(agg) > length(terms) + 1) {
  fit <- tryCatch(aov(as.formula(paste("y ~", paste(terms, collapse = " + "))), data = agg), error = function(e) NULL)
  if (!is.null(fit)) {
    a <- summary(fit)[[1]]; ss <- a[["Sum Sq"]]
    vd <- data.frame(factor = sub("Residuals", "residual", trimws(rownames(a))),
                     variance_pct = round(100 * ss / sum(ss), 2), stringsAsFactors = FALSE)
    write.csv(vd, file.path(OUT, "variance_decomposition.csv"), row.names = FALSE)
    cat("variance_decomposition.csv\n"); print(vd, row.names = FALSE)
  }
}

# ---- Figures --------------------------------------------------------------
if (have_gg) {
  B <- load_boundary(agg)                       # boundary underlay, or NULL
  under <- if (!is.null(B)) B$layer else NULL    # drawn before the points
  geo   <- if (!is.null(B)) B$coord else coord_quickmap()
  if (!is.null(B)) cat(sprintf("map boundary: %s\n", BOUNDARY))
  # 4a. Faceted yield maps
  pm <- ggplot(agg, aes(longitude, latitude, colour = y)) + under +
    geom_point(size = 1.6) + facet_wrap(~combo) + geo + scale_colour_viridis_c(option = "D") +
    labs(title = paste("Mean", RV, "by input combination"), colour = RV, x = NULL, y = NULL) +
    theme_minimal(base_size = 10)
  ggsave(file.path(OUT, "fig_yield_maps.png"), pm, width = 11, height = 7.5, dpi = 150)
  # 4b. Sensitivity (CV) map
  ps <- ggplot(sens, aes(longitude, latitude, colour = cv_pct)) + under +
    geom_point(size = 2) + geo + scale_colour_viridis_c(option = "B") +
    labs(title = paste0("Input sensitivity: across-combination CV of ", RV),
         subtitle = "higher = the choice of weather/soil source moves the result more here",
         colour = "CV (%)", x = NULL, y = NULL) + theme_minimal(base_size = 11)
  ggsave(file.path(OUT, "fig_sensitivity_cv_map.png"), ps, width = 8, height = 6, dpi = 150)
  # 4c. Variance decomposition bar
  if (!is.null(vd)) {
    pv <- ggplot(vd, aes(x = "", y = variance_pct, fill = factor)) + geom_col(width = 0.6) +
      coord_flip() + labs(title = paste("Variance in", RV, "explained by each factor"),
        x = NULL, y = "share of variance (%)", fill = NULL) + theme_minimal(base_size = 12)
    ggsave(file.path(OUT, "fig_variance_decomposition.png"), pv, width = 8, height = 2.6, dpi = 150)
  }
  # 4d. Pairwise RMSD + rank-stability heatmaps
  heat <- function(M, title, lab, pal) {
    df <- expand.grid(a = rownames(M), b = colnames(M)); df$v <- as.vector(M)
    ggplot(df, aes(a, b, fill = v)) + geom_tile() +
      geom_text(aes(label = ifelse(is.na(v), "", formatC(v, format = "g", digits = 2))), size = 2.6) +
      scale_fill_gradientn(colours = pal, name = lab) + labs(title = title, x = NULL, y = NULL) +
      theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  }
  ggsave(file.path(OUT, "fig_pairwise_rmsd.png"),
         heat(rmsd, paste("Pairwise RMSD of", RV, "between combinations"), "RMSD", c("#ffffcc","#fd8d3c","#bd0026")),
         width = 8, height = 6.5, dpi = 150)
  ggsave(file.path(OUT, "fig_rank_stability.png"),
         heat(round(cor_s, 2), "Spatial rank stability (Spearman) between combinations", "rho", c("#b2182b","#f7f7f7","#2166ac")),
         width = 8, height = 6.5, dpi = 150)
  cat("figures: fig_yield_maps.png, fig_sensitivity_cv_map.png, fig_pairwise_rmsd.png, fig_rank_stability.png",
      if (!is.null(vd)) ", fig_variance_decomposition.png" else "", "\n", sep = "")
} else cat("(ggplot2 not installed — figures skipped)\n")

cat(sprintf("\nDone. Outputs in %s/\n", OUT))
