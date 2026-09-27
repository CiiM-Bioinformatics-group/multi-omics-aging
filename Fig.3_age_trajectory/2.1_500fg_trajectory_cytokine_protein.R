## Fig 1B (500FG): cytokine vs. protein aging-trajectory schematic --
## LOESS trajectory curves (with pointwise 95% CI ribbons) and each
## layer's transition point from an age-binned, weighted segmented
## regression.

library(dplyr)
library(ggplot2)
library(segmented)
source("header.R")

## ---- paths (override via env vars; defaults are placeholders) -------------
cfg <- list(
  merge_file  = Sys.getenv("TRAJECTORY_MERGE_FILE", "data/500fg_heatmap_merge_data_gender.csv"),
  output_stem = Sys.getenv("TRAJECTORY_B_OUTPUT_STEM", "output/fig2b_trajectory_500fg_nature_final")
)
dir.create(dirname(cfg$output_stem), recursive = TRUE, showWarnings = FALSE)

## ---- load, split into cytokine / protein blocks by column name range ------
merge_data <- read.csv(cfg$merge_file, row.names = 1)
cn <- colnames(merge_data)
cytokine <- merge_data[, match("IFNy_C.conidiaHK_WB_48h", cn):match("IL22_Bacteroides_PBMC_7days", cn)]
protein  <- merge_data[, (match("IL22_Bacteroides_PBMC_7days", cn) + 1):ncol(merge_data)]
Age <- merge_data$Age

## ---- 1. continuous LOESS trajectory (per-sample z-scored row-mean score) --
cyt_score  <- rowMeans(scale(cytokine), na.rm = TRUE)
prot_score <- rowMeans(scale(protein),  na.rm = TRUE)
df1 <- data.frame(Age, cyt_score, prot_score)
df1 <- df1[complete.cases(df1), ]

age_seq  <- seq(min(df1$Age), max(df1$Age), length.out = 200)
fit_cyt  <- loess(cyt_score  ~ Age, data = df1, span = 0.75)
fit_prot <- loess(prot_score ~ Age, data = df1, span = 0.75)
pred_cyt  <- predict(fit_cyt,  newdata = data.frame(Age = age_seq))
pred_prot <- predict(fit_prot, newdata = data.frame(Age = age_seq))

## Normalize to 0-1 for the schematic y-axis.
pred_cyt_norm  <- (pred_cyt  - min(pred_cyt))  / diff(range(pred_cyt))
pred_prot_norm <- (pred_prot - min(pred_prot)) / diff(range(pred_prot))

## ---- 2. transition point: age-binned, weighted segmented regression -------
cytokine_num <- as.data.frame(lapply(cytokine, as.numeric))
protein_num  <- as.data.frame(lapply(protein,  as.numeric))
cyt_z  <- scale(as.matrix(cytokine_num))
prot_z <- scale(as.matrix(protein_num))

bin_width <- 5
breaks <- seq(floor(min(Age, na.rm = TRUE) / bin_width) * bin_width,
              ceiling(max(Age, na.rm = TRUE) / bin_width) * bin_width,
              by = bin_width)
df_age <- data.frame(Age = Age) %>%
  mutate(bin = cut(Age, breaks = breaks, include.lowest = TRUE, right = FALSE))

## Mean expression per feature within each age bin, averaged across
## features into one layer score per bin (bins with < 5 samples dropped).
get_bin_traj <- function(zmat, df_age) {
  split_idx <- split(seq_len(nrow(df_age)), df_age$bin)
  split_idx <- split_idx[sapply(split_idx, length) >= 5]
  centers <- sapply(names(split_idx), function(b) {
    lr <- as.numeric(gsub("\\[|\\)|\\]|\\(", "", unlist(strsplit(b, ","))))
    mean(lr)
  })
  n_bin <- sapply(split_idx, length)
  layer_score <- sapply(split_idx, function(ix) {
    mean(colMeans(zmat[ix, , drop = FALSE], na.rm = TRUE), na.rm = TRUE)
  })
  data.frame(Age = as.numeric(centers), score = as.numeric(layer_score), n = as.numeric(n_bin))
}
traj_cyt  <- get_bin_traj(cyt_z,  df_age)
traj_prot <- get_bin_traj(prot_z, df_age)

## Breakpoint from the bin-level (sample-size-weighted) trajectory;
## psi_init values are the visually-estimated priors from the source
## analysis.
get_tau <- function(traj, psi_init) {
  fit_seg <- segmented(lm(score ~ Age, data = traj, weights = n), seg.Z = ~ Age, psi = psi_init)
  if ("Est." %in% colnames(fit_seg$psi)) as.numeric(fit_seg$psi[1, "Est."]) else as.numeric(fit_seg$psi[1, 2])
}
tau_cyt  <- get_tau(traj_cyt,  psi_init = 45)
tau_prot <- get_tau(traj_prot, psi_init = 60)

## ---- 3. pointwise 95% CIs on the LOESS curves, then build + export --------
p <- base_params()

## Re-derives each curve's CI from its loess fit, checking that the fit
## still reproduces the already-normalized curve from step 1 (guards
## against silently drifting from the fit above).
predict_norm_ci <- function(model, reference, label) {
  pr <- predict(model, newdata = data.frame(Age = age_seq), se = TRUE)
  r <- range(pr$fit, na.rm = TRUE)
  d <- diff(r)
  if (!is.finite(d) || d <= 0) stop("The fitted trajectory has no valid range.")
  value <- (pr$fit - r[1]) / d
  if (!isTRUE(all.equal(as.numeric(value), as.numeric(reference), tolerance = 1e-7)))
    stop("The model does not match the current normalized curve. Restore the original fit.")
  delta <- qt(0.975, df = pr$df) * pr$se.fit / d
  data.frame(Age = age_seq, value = value, lo = value - delta, hi = value + delta, type = label)
}

df_plot <- rbind(
  predict_norm_ci(fit_cyt, pred_cyt_norm, "Cytokine responses (Immunosenescence)"),
  predict_norm_ci(fit_prot, pred_prot_norm, "Protein (Inflammaging)")
)
df_plot$type <- factor(df_plot$type, levels = c("Cytokine responses (Immunosenescence)", "Protein (Inflammaging)"))

TRAJ_COLORS <- c("Cytokine responses (Immunosenescence)" = "#e63946", "Protein (Inflammaging)" = "#457b9d")
y_mid <- mean(range(c(df_plot$lo, df_plot$hi, df_plot$value), na.rm = TRUE))  # centers the transition labels

p1 <- ggplot(df_plot, aes(x = Age, y = value, colour = type)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = type), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_vline(xintercept = c(tau_cyt, tau_prot), linetype = 2, linewidth = 0.3, color = "gray40") +
  annotate("text", x = tau_cyt + 1, y = y_mid,
           label = paste0("Cytokine transition (", round(tau_cyt, 1), " y)"),
           size = p$FONT_PT / .pt, color = "black", angle = 90, hjust = 0.5, vjust = 0.5) +
  annotate("text", x = tau_prot + 1, y = y_mid,
           label = paste0("Protein transition (", round(tau_prot, 1), " y)"),
           size = p$FONT_PT / .pt, color = "black", angle = 90, hjust = 0.5, vjust = 0.5) +
  scale_color_manual(values = TRAJ_COLORS, name = NULL) +
  scale_fill_manual(values = TRAJ_COLORS, guide = "none") +
  scale_x_continuous(breaks = AGE_BREAKS, expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
  coord_cartesian(clip = "off") +
  labs(x = "Age (years)", y = "Normalized composite score") +
  common_ggplot_theme(p$FONT_PT) +
  theme(
    legend.position = "bottom",
    legend.box.spacing = grid::unit(0.3, "mm"),
    legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0)
  ) +
  guides(colour = guide_legend(ncol = 2))

## Exported at TARGET_WIDTH_MM (the shared journal-grid width), content
## padded by PAGE_MARGIN_MM on each side; height is tuned directly since
## a single trajectory panel has no natural "fit" height.
PAGE_MARGIN_MM <- 1.5
B_WIDTH_MM  <- TARGET_WIDTH_MM - 2 * PAGE_MARGIN_MM
B_HEIGHT_MM <- 40

page_w_mm <- B_WIDTH_MM  + 2 * PAGE_MARGIN_MM   # == TARGET_WIDTH_MM
page_h_mm <- B_HEIGHT_MM + 2 * PAGE_MARGIN_MM
w_in <- page_w_mm / 25.4
h_in <- page_h_mm / 25.4

draw_panel_B <- function() {
  grid::grid.newpage()
  print(p1, newpage = FALSE, vp = grid::viewport(
    x = unit(PAGE_MARGIN_MM, "mm"), y = unit(PAGE_MARGIN_MM, "mm"),
    width = unit(B_WIDTH_MM, "mm"), height = unit(B_HEIGHT_MM, "mm"),
    just = c("left", "bottom")
  ))
}

open_pdf_device(paste0(cfg$output_stem, ".pdf"), width = w_in, height = h_in)
draw_panel_B()
grDevices::dev.off()

grDevices::png(paste0(cfg$output_stem, ".png"), width = w_in * PNG_DPI, height = h_in * PNG_DPI, res = PNG_DPI)
draw_panel_B()
grDevices::dev.off()

message(sprintf("%s.pdf / .png saved at %.1f x %.1f mm (content %.1f x %.1f mm + %.1f mm margin).",
                 cfg$output_stem, page_w_mm, page_h_mm, B_WIDTH_MM, B_HEIGHT_MM, PAGE_MARGIN_MM))
