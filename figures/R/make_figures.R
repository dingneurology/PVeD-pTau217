#!/usr/bin/env Rscript
# Figures 2 and 3, Supplementary Figure S1, and the data panels of Figure 1 (sections 3 and 4).
# Reads the CSV files in ./source_data (written by figures/export_fig_sources.py, i.e. after run_all.py)
# and writes PDF/PNG/TIFF files to ./output.
#
# Packages: ggplot2, dplyr, readr, ggseg (>= 2.0 tested; 1.6 layout fallback included), scales.
# Layout uses base R grid; patchwork is not required.
# Usage, from the repository root:   Rscript figures/R/make_figures.R
#   or in R/RStudio:                  source("figures/R/make_figures.R")
#   Figure 3 only:                    options(fig3_only = TRUE); source("figures/R/make_figures.R")

needed <- c("ggplot2", "dplyr", "readr", "ggseg", "scales")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Please install missing packages first: install.packages(c(",
       paste0('"', missing, '"', collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(dplyr); library(readr); library(ggseg); library(scales)
})

# locate this script's folder
get_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f) == 1) return(dirname(normalizePath(f)))
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of)) return(dirname(normalizePath(of)))
  getwd()
}
base_dir <- get_dir()
src <- file.path(base_dir, "source_data")
out_dir <- file.path(base_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) readr::read_csv(file.path(src, f), show_col_types = FALSE)

font_family <- if (.Platform$OS.type == "windows") "sans" else "Arial"
pal <- c(neutral_dark = "#272727", neutral_mid = "#767676", signal_blue = "#2F6EA6",
         accent_red = "#C9483A", accent_orange = "#D98B2B", signal_teal = "#3BA99C")
profile_levels <- c("Low pTau217 / High PVeD", "Low pTau217 / Low PVeD", "High pTau217 / High PVeD", "High pTau217 / Low PVeD")
profile_cols <- c("Low pTau217 / High PVeD" = "#72A7D8", "Low pTau217 / Low PVeD" = "#9C9C9C",
                  "High pTau217 / High PVeD" = "#E2A85E", "High pTau217 / Low PVeD" = "#C9483A")
profile_short <- c("Low pTau217 / High PVeD" = "pTau217−/\nPVeD+", "Low pTau217 / Low PVeD" = "pTau217−/\nPVeD−",
                   "High pTau217 / High PVeD" = "pTau217+/\nPVeD+", "High pTau217 / Low PVeD" = "pTau217+/\nPVeD−")

theme_pub <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
          axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.text = element_text(size = base_size - 0.45, colour = "black"),
          legend.text = element_text(size = base_size - 0.6),
          strip.text = element_text(size = base_size - 0.1, face = "bold"),
          strip.background = element_rect(linewidth = 0.35),
          plot.title = element_text(size = base_size + 0.6, face = "bold", hjust = 0),
          plot.caption = element_text(size = base_size - 0.9, hjust = 1, colour = pal[["neutral_dark"]]),
          panel.grid = element_blank(), plot.margin = margin(4, 4, 4, 4))
}
theme_set(theme_pub())

# arrange a list of ggplots on a grid; cells = list(list(plot, row(s), col(s)), ...)
draw_grid <- function(cells, nrow, ncol, heights = rep(1, nrow), widths = rep(1, ncol)) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow, ncol, heights = unit(heights, "null"), widths = unit(widths, "null"))))
  for (cl in cells) {
    vp <- viewport(layout.pos.row = cl[[2]], layout.pos.col = cl[[3]])
    if (is.function(cl[[1]])) { pushViewport(vp); cl[[1]](); popViewport() } else print(cl[[1]], vp = vp)
  }
  popViewport()
}
save_fig <- function(draw, name, w_mm, h_mm, tiff_dpi = 600) {
  w <- w_mm / 25.4; h <- h_mm / 25.4
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  pdf_dev(file.path(out_dir, paste0(name, ".pdf")), width = w, height = h); draw(); dev.off()
  grDevices::png(file.path(out_dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 300, bg = "white"); draw(); dev.off()
  grDevices::tiff(file.path(out_dir, paste0(name, ".tiff")), width = w, height = h, units = "in", res = tiff_dpi,
                  compression = "lzw", bg = "white"); draw(); dev.off()
  message("Saved ", name, " (pdf/png/tiff) in ", out_dir)
}

# ============================== Figure 3 ===================================
forest <- rd("Fig3_forest.csv")
brain <- rd("Fig3_brainmap.csv")
pred <- rd("Fig3_prediction.csv")
sens <- rd("Fig3_sensitivity.csv")

dk_atlas <- if (is.function(dk)) dk() else dk   # ggseg >= 2.0 exports dk() as a function
# Regional maps: each hemisphere/view is drawn as its own small plot so that the four views can be
# placed close together (ggseg's combined layout leaves wide gaps between hemispheres).
brain_fill <- scale_fill_gradient(low = "#0F4C81", high = "#D6E6F2", limits = c(-0.12, 0), oob = squish,
                                  breaks = c(-0.12, -0.09, -0.06, -0.03, 0), na.value = "#E3E3E3",
                                  name = "Interaction β")
geom_bbox_w <- function(p) {
  pp <- ggplot_build(p)$layout$panel_params[[1]]
  xr <- if (!is.null(pp$x_range)) pp$x_range else if (!is.null(pp$x.range)) pp$x.range else pp$x$continuous_range
  diff(range(xr))
}
full_w <- NULL
one_view <- function(d, h, v) {
  if (is.null(full_w)) full_w <<- geom_bbox_w(ggplot() + geom_brain(data = d, atlas = dk_atlas, aes(fill = brain_beta)))
  for (arg in c("view", "side")) {
    a <- list(data = d, atlas = dk_atlas, mapping = aes(fill = brain_beta), colour = "grey25",
              linewidth = 0.12, hemi = h)
    a[[arg]] <- v
    p <- suppressWarnings(ggplot() + do.call(geom_brain, a))
    w <- tryCatch(suppressWarnings(geom_bbox_w(p)), error = function(e) NA)
    if (!is.na(w) && w < 0.45 * full_w) {
      return(p + brain_fill + coord_fixed(expand = FALSE) + theme_void(base_family = font_family) +
               theme(legend.position = "none", plot.margin = margin(0, 1, 0, 1)))
    }
  }
  stop("Could not draw a single hemisphere/view with this ggseg version.")
}
brain_views <- list(c("left", "lateral"), c("left", "medial"), c("right", "medial"), c("right", "lateral"))
brain_legend <- function() {
  d <- brain %>% filter(cohort == "ADNI") %>% select(region, brain_beta)
  p <- ggplot() + geom_brain(data = d, atlas = dk_atlas, aes(fill = brain_beta)) + brain_fill +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0, barheight = unit(26, "mm"),
                                  barwidth = unit(3, "mm"), ticks.colour = "black", frame.colour = "black")) +
    theme_void(base_size = 6.5, base_family = font_family) +
    theme(legend.title = element_text(size = 6.4, face = "bold"), legend.text = element_text(size = 6))
  g <- ggplotGrob(p)
  k <- grep("guide-box", g$layout$name)
  k <- k[!vapply(g$grobs[k], function(x) inherits(x, "zeroGrob"), logical(1))][1]
  g$grobs[[k]]
}
brain_band <- function() {
  ttl <- c("ADNI" = "A   ADNI: regional tau PET (flortaucipir)", "HABS-HD" = "E   HABS-HD: regional tau PET (PI-2620)")
  vlab <- c("Left lateral", "Left medial", "Right medial", "Right lateral")
  leg <- brain_legend()
  pushViewport(viewport(layout = grid.layout(5, 6,
    heights = unit(c(5, 1, 5, 1, 4), c("mm", "null", "mm", "null", "mm")),
    widths = unit(c(2, 1, 1, 1, 1, 24), c("mm", "null", "null", "null", "null", "mm")))))
  for (k in seq_along(ttl)) {
    coh <- names(ttl)[k]; r <- 2 * k - 1
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = 2:5))
    grid.text(ttl[[k]], x = unit(0, "npc"), y = unit(0.35, "npc"), just = "left",
              gp = gpar(fontsize = 7.1, fontface = "bold", fontfamily = font_family))
    popViewport()
    d <- brain %>% filter(cohort == coh) %>% select(region, brain_beta)
    for (j in seq_along(brain_views)) {
      print(one_view(d, brain_views[[j]][1], brain_views[[j]][2]), vp = viewport(layout.pos.row = r + 1, layout.pos.col = j + 1))
    }
  }
  for (j in seq_along(vlab)) {
    pushViewport(viewport(layout.pos.row = 5, layout.pos.col = j + 1))
    grid.text(vlab[j], y = unit(0.55, "npc"), gp = gpar(fontsize = 5.8, col = "grey30", fontfamily = font_family))
    popViewport()
  }
  pushViewport(viewport(layout.pos.row = 1:5, layout.pos.col = 6)); grid.draw(leg); popViewport()
  pushViewport(viewport(layout.pos.row = 4:5, layout.pos.col = 6))
  grid.rect(x = unit(3, "mm"), y = unit(0.42, "npc"), width = unit(3, "mm"), height = unit(2.4, "mm"),
            gp = gpar(fill = "#E3E3E3", col = "grey25", lwd = 0.4))
  grid.text("Not analysed", x = unit(5.6, "mm"), y = unit(0.42, "npc"), just = "left",
            gp = gpar(fontsize = 5.8, fontfamily = font_family))
  popViewport()
  popViewport()
}

forest_panel <- function(coh, letter) {
  d <- forest %>% filter(cohort == coh) %>% arrange(order) %>%
    mutate(star = ifelse(sig, "*", ""),
           lab = sprintf("%.3f (%.3f, %.3f)%s", beta, lo, hi, star),
           outcome_short = factor(outcome_short, levels = rev(unique(outcome_short))))
  xr <- max(abs(c(d$lo, d$hi))) * 1.08
  ggplot(d, aes(beta, outcome_short, colour = domain_group, alpha = sig)) +
    geom_vline(xintercept = 0, linewidth = 0.35, linetype = "22", colour = pal[["neutral_mid"]]) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.45) +
    geom_point(size = 1.8) +
    geom_text(aes(x = xr * 1.05, label = lab), hjust = 0, size = 1.9, alpha = 1, show.legend = FALSE) +
    annotate("text", x = xr * 1.05, y = Inf, label = "β (95% CI)", hjust = 0, vjust = -0.3, size = 2.0,
             fontface = "bold", colour = pal[["neutral_dark"]]) +
    scale_colour_manual(values = c("AD PET pathology" = pal[["accent_red"]], "Other outcomes" = pal[["neutral_mid"]])) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.45), guide = "none") +
    coord_cartesian(xlim = c(-xr, xr), clip = "off") +
    labs(x = expression(beta ~ "for pTau217" %*% "PVeD interaction"), y = NULL,
         title = paste0(letter, "  ", coh, ": all outcomes"), caption = "* FDR-adjusted P < 0.05") +
    theme(legend.position = "none", plot.margin = margin(10, 62, 4, 4))
}

pred_panel <- function(coh, letter) {
  d <- pred %>% filter(cohort == coh) %>%
    mutate(pved_level = factor(pved_level, levels = c("Higher PVeD (+1 SD)", "Lower PVeD (-1 SD)")),
           outcome = factor(outcome, levels = unique(outcome)))
  cols <- c("Higher PVeD (+1 SD)" = pal[["signal_blue"]], "Lower PVeD (-1 SD)" = pal[["accent_red"]])
  ggplot(d, aes(z_ptau217, fit, colour = pved_level, fill = pved_level)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ outcome, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
    labs(x = expression("Standardized " * log[10] * " plasma pTau217"), y = "Adjusted standardized outcome",
         title = paste0(letter, "  ", coh, ": model-predicted outcomes"), colour = NULL, fill = NULL) +
    theme(legend.position = "bottom")
}

sens_panel <- function(coh, letter) {
  lv <- c("Primary (all available)", "±24 months", "±12 months", "±6 months", "+ interval covariates",
          "M2 harmonized covariates", "M3 + CDR-SB")
  d <- sens %>% filter(cohort == coh) %>%
    mutate(model = factor(model, levels = rev(lv[lv %in% model])),
           outcome_short = factor(outcome_short, levels = unique(outcome_short)),
           star = ifelse(fdr < 0.05, "*", ""))
  cols <- setNames(c(pal[["signal_blue"]], pal[["accent_orange"]], pal[["accent_red"]]), levels(d$outcome_short))
  pd <- position_dodge(width = 0.6)
  ggplot(d, aes(beta, model, colour = outcome_short)) +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "22", colour = pal[["neutral_mid"]]) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.35, position = pd) +
    geom_point(size = 1.3, position = pd) +
    geom_text(aes(x = hi, label = star), hjust = -0.4, vjust = 0.75, size = 2.4, position = pd, show.legend = FALSE) +
    scale_colour_manual(values = cols) +
    labs(x = expression(beta ~ "for pTau217" %*% "PVeD interaction"), y = NULL,
         title = paste0(letter, "  ", coh, ": temporal alignment and adjustment"), colour = NULL,
         caption = "* FDR-adjusted P < 0.05") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.key.width = unit(4.5, "mm"), legend.margin = margin(0, 0, 0, 0))
}

# Layout: regional maps for both cohorts across the top (A, E); below, ADNI in the left column (B-D)
# and HABS-HD in the right column (F-H), so the two cohorts can be compared row by row.
fig3_draw <- function() {
  draw_grid(list(
    list(brain_band, 1, 1:2),
    list(forest_panel("ADNI", "B"), 2, 1), list(forest_panel("HABS-HD", "F"), 2, 2),
    list(pred_panel("ADNI", "C"), 3, 1),   list(pred_panel("HABS-HD", "G"), 3, 2),
    list(sens_panel("ADNI", "D"), 4, 1),   list(sens_panel("HABS-HD", "H"), 4, 2)),
    nrow = 4, ncol = 2, heights = c(62, 70, 60, 64))
}
save_fig(fig3_draw, "Figure_3_revised", 183, 258, tiff_dpi = 300)
if (isTRUE(getOption("fig3_only", FALSE))) stop_quietly <- TRUE else stop_quietly <- FALSE

if (!stop_quietly) {
# ======================= Supplementary Figure S1 ===========================
sc <- rd("FigS1_scatter.csv") %>% mutate(profile = factor(profile, levels = profile_levels))
pc <- rd("FigS1_profile_counts.csv")
ep <- rd("FigS1_endpoints.csv")

scatter_panel <- function(coh, letter) {
  d <- sc %>% filter(cohort == coh)
  med <- pc %>% filter(cohort == coh, profile == "__medians__")
  ggplot(d, aes(log10_ptau217, pved_total, colour = profile)) +
    geom_point(size = 0.35, alpha = 0.55) +
    geom_vline(xintercept = med$ptau_median, linetype = "22", linewidth = 0.35) +
    geom_hline(yintercept = med$pved_median, linetype = "22", linewidth = 0.35) +
    scale_colour_manual(values = profile_cols, guide = "none") +
    labs(x = expression(log[10] ~ "(plasma pTau217)"), y = "PVeD total", title = paste0(letter, "  ", coh))
}
count_panel <- function(coh, letter) {
  d <- pc %>% filter(cohort == coh, profile != "__medians__") %>% mutate(profile = factor(profile, levels = profile_levels))
  ggplot(d, aes(profile, n, fill = profile)) +
    geom_col(width = 0.7) + geom_text(aes(label = n), vjust = -0.4, size = 2.0) +
    scale_fill_manual(values = profile_cols, guide = "none") +
    scale_x_discrete(labels = profile_short) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Participants, n", title = paste0(letter, "  ", coh, ": profile sizes"))
}
endpoint_panel <- function(coh, letter) {
  d <- ep %>% filter(cohort == coh) %>% arrange(order) %>% mutate(label = factor(label, levels = label))
  ggplot(d, aes(label, n)) +
    geom_col(width = 0.65, fill = pal[["signal_blue"]]) + geom_text(aes(label = n), vjust = -0.4, size = 2.0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Participants, n", title = paste0(letter, "  ", coh, ": endpoint-specific samples"))
}
figS1_draw <- function() {
  draw_grid(list(
    list(scatter_panel("ADNI", "A"), 1, 1), list(count_panel("ADNI", "B"), 1, 2), list(endpoint_panel("ADNI", "C"), 2, 1:2),
    list(scatter_panel("HABS-HD", "D"), 3, 1), list(count_panel("HABS-HD", "E"), 3, 2), list(endpoint_panel("HABS-HD", "F"), 4, 1:2)),
    nrow = 4, ncol = 2, heights = c(1.15, 0.85, 1.15, 0.85))
}
save_fig(figS1_draw, "Supplementary_Figure_S1_revised", 183, 247)


# ============================== Figure 2 ===================================
f2 <- rd("Fig2_participants.csv") %>% mutate(profile = factor(profile, levels = profile_levels))
dx_cols <- c("CN" = "#9DBFD3", "CU" = "#9DBFD3", "MCI" = "#D9B96E", "AD dementia" = "#B95A4A", "Dementia" = "#B95A4A")
dx_panel <- function(coh, letter) {
  lv <- if (coh == "ADNI") c("CN", "MCI", "AD dementia") else c("CU", "MCI", "Dementia")
  d <- f2 %>% filter(cohort == coh, !is.na(dx)) %>% count(profile, dx) %>% group_by(profile) %>%
    mutate(pct = 100 * n / sum(n)) %>% ungroup() %>% mutate(dx = factor(dx, levels = rev(lv)))
  ggplot(d, aes(profile, pct, fill = dx)) + geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = ifelse(pct >= 6, sprintf("%.0f%%", pct), "")), position = position_stack(vjust = 0.5), size = 1.9) +
    scale_fill_manual(values = dx_cols, breaks = lv) + scale_x_discrete(labels = profile_short) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Clinical diagnosis, %", fill = NULL, title = paste0(letter, "  ", coh, ": clinical composition")) +
    theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"))
}
box_panel <- function(coh, letter, vars, labs_, ttl) {
  d <- bind_rows(lapply(seq_along(vars), function(k) f2 %>% filter(cohort == coh) %>%
                          transmute(profile, value = .data[[vars[k]]], measure = labs_[k]))) %>%
    filter(!is.na(value)) %>% mutate(measure = factor(measure, levels = labs_))
  ggplot(d, aes(profile, value, colour = profile)) +
    geom_jitter(width = 0.18, size = 0.3, alpha = 0.35) +
    geom_boxplot(outlier.shape = NA, fill = NA, width = 0.55, linewidth = 0.35, colour = "black") +
    facet_wrap(~ measure, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = profile_cols, guide = "none") + scale_x_discrete(labels = profile_short) +
    labs(x = NULL, y = NULL, title = paste0(letter, "  ", coh, ": ", ttl)) +
    theme(axis.text.x = element_text(size = 4.6))
}
fig2_draw <- function() {
  tau_lab_a <- "Tau PET, meta-temporal SUVR"; tau_lab_h <- "Tau PET, MsTL stage SUVR"
  draw_grid(list(
    list(dx_panel("ADNI", "A"), 1, 1),
    list(box_panel("ADNI", "B", c("amyloid_centiloid", "tau_suvr"), c("Amyloid PET, Centiloid", tau_lab_a), "amyloid and tau PET"), 1, 2),
    list(box_panel("ADNI", "C", c("log10_gfap", "log10_nfl"), c("log10 GFAP", "log10 NfL"), "plasma GFAP and NfL"), 2, 1),
    list(box_panel("ADNI", "D", c("pacc"), c("Baseline PACC"), "baseline PACC"), 2, 2),
    list(dx_panel("HABS-HD", "E"), 3, 1),
    list(box_panel("HABS-HD", "F", c("amyloid_centiloid", "tau_suvr"), c("Amyloid PET, Centiloid", tau_lab_h), "amyloid and tau PET"), 3, 2),
    list(box_panel("HABS-HD", "G", c("log10_gfap", "log10_nfl"), c("log10 GFAP", "log10 NfL"), "plasma GFAP and NfL"), 4, 1),
    list(box_panel("HABS-HD", "H", c("pacc"), c("Baseline PACC-like composite"), "baseline PACC"), 4, 2)),
    nrow = 4, ncol = 2)
}
save_fig(fig2_draw, "Figure_2_revised", 183, 240, tiff_dpi = 300)

# ============ Figure 1 replacement data panels (sections 3 and 4) ============
sl <- rd("Fig1_profile_slopes.csv")
fig1_s3_draw <- function() {
  draw_grid(list(list(forest_panel("ADNI", "") + labs(title = NULL) + theme(plot.margin = margin(14, 70, 4, 4)), 1, 1),
                 list(pred_panel("ADNI", "") + labs(title = NULL) + guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)), 1, 2)), nrow = 1, ncol = 2, widths = c(1.35, 1))
}
save_fig(fig1_s3_draw, "Figure_1_section3_panels_revised", 140, 76, tiff_dpi = 300)
slope_panel <- function(coh) {
  d <- sl %>% filter(cohort == coh, term %in% profile_levels) %>% mutate(term = factor(term, levels = profile_levels))
  ctr <- sl %>% filter(cohort == coh, grepl("^Contrast", term))
  ggplot(d, aes(term, est, colour = term)) +
    geom_hline(yintercept = 0, linetype = "22", linewidth = 0.3, colour = pal[["neutral_mid"]]) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5) + geom_point(size = 1.6) +
    scale_colour_manual(values = profile_cols, guide = "none") + scale_x_discrete(labels = profile_short) +
    labs(x = NULL, y = "Annual PACC change, SD/year", title = paste0(coh, " (PACC)"),
         caption = sprintf("Key contrast pTau217+/PVeD− vs pTau217+/PVeD+:\n%.3f (%.3f, %.3f); P = %.3f", ctr$est, ctr$lo, ctr$hi, ctr$p)) +
    theme(plot.caption = element_text(hjust = 0.5, size = 5.4))
}
fig1_s4_draw <- function() {
  draw_grid(list(list(slope_panel("ADNI"), 1, 1), list(slope_panel("HABS-HD"), 1, 2)), nrow = 1, ncol = 2)
}
save_fig(fig1_s4_draw, "Figure_1_section4_panels_revised", 112, 58, tiff_dpi = 300)

}
message("Done.")
