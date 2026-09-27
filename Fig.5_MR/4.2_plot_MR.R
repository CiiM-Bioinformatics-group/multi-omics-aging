# 5.2 Combined MR scatter plot (4 EAA outcomes) from the outputs of 5.1
#     Run from the repository root: Rscript MR/5.2_plot_MR.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

mr_dir <- "results/MR/ferulic_4_sulfate_p1e-06_kb10000_r20.5"
combined_dir <- file.path(mr_dir, "MR_combined_results")
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

outcome_labels <- c(PhenoAge_AA = "PhenoAge AA", Hannum_AA = "Hannum AA", GrimAge_AA = "GrimAge AA", IEAA = "IEAA")
tsmr_cols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
               "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")

# ---- Load per-outcome SNP effects and MR estimates ----
snp_list <- list()
line_list <- list()
for (d in list.dirs(mr_dir, recursive = FALSE)) {
  if (!startsWith(basename(d), "MR_output_")) next
  outcome <- sub("^MR_output_", "", basename(d))
  dat    <- read.csv(file.path(d, "harmonised_data_with_F.csv"), stringsAsFactors = FALSE)
  mr_res <- read.csv(file.path(d, "MR_all_methods.csv"), stringsAsFactors = FALSE)
  if ("mr_keep" %in% names(dat)) dat <- dat[dat$mr_keep %in% c(TRUE, "TRUE", 1), , drop = FALSE]

  # Orient SNPs to positive exposure effect (as in TwoSampleMR::mr_scatter_plot)
  flip <- dat$beta.exposure < 0
  dat$beta.exposure[flip] <- -dat$beta.exposure[flip]
  dat$beta.outcome[flip]  <- -dat$beta.outcome[flip]

  # Line intercepts: 0, except MR-Egger (weighted regression intercept)
  mr_res$a <- 0
  egger <- mr_res$method %in% c("MR Egger", "MR Egger (bootstrap)")
  if (any(egger)) {
    fit <- lm(beta.outcome ~ beta.exposure, data = dat, weights = 1 / dat$se.outcome^2)
    mr_res$a[egger] <- unname(coef(fit)[1])
  }
  snp_list[[outcome]]  <- mutate(dat, outcome_name = outcome)
  line_list[[outcome]] <- mutate(mr_res, outcome_name = outcome)
}

present <- intersect(names(outcome_labels), names(snp_list))
lab <- function(x) factor(outcome_labels[x], levels = unname(outcome_labels[present]))
snp_df  <- bind_rows(snp_list[present]) %>% mutate(outcome_lab = lab(outcome_name))
line_df <- bind_rows(line_list[present]) %>% mutate(outcome_lab = lab(outcome_name))

method_levels <- unique(unlist(lapply(line_list[present], `[[`, "method")))
line_df$method <- factor(line_df$method, levels = method_levels)
pal <- setNames(tsmr_cols[seq_along(method_levels)], method_levels)

# ---- Plot ----
p_scatter <- ggplot(snp_df, aes(x = beta.exposure, y = beta.outcome)) +
  geom_errorbar(aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome),
                colour = "grey", width = 0, linewidth = 0.25) +
  geom_errorbarh(aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure),
                 colour = "grey", width = 0, linewidth = 0.25) +
  geom_point(size = 0.7, stroke = 0.2) +
  geom_abline(data = line_df, aes(intercept = a, slope = b, colour = method), linewidth = 0.45, show.legend = TRUE) +
  scale_colour_manual(values = pal, drop = FALSE) +
  facet_wrap(~ outcome_lab, scales = "free", ncol = 2) +
  labs(colour = "MR Test", x = "SNP effect on ferulic 4-sulfate", y = "SNP effect on outcome") +
  theme_bw(base_size = 5) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.2),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(size = 5, colour = "black"),
    axis.title        = element_text(size = 5, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.2, colour = "black"),
    axis.ticks.length = unit(0.05, "cm"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    strip.clip        = "off",
    strip.text        = element_text(size = 5, colour = "black", margin = margin(1.2, 1.2, 1.2, 1.2)),
    strip.background  = element_rect(fill = "grey96", colour = "black", linewidth = 0.4),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = 5, colour = "black"),
    legend.text       = element_text(size = 5, colour = "black"),
    legend.key.size   = unit(0.28, "cm"),
    legend.spacing.x  = unit(0.12, "cm"),
    legend.margin     = margin(t = 1, b = 0),
    panel.spacing     = unit(0.6, "lines"),
    plot.margin       = margin(3, 4, 3, 3)
  ) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE))

ggsave(file.path(combined_dir, "combined_scatter_plot.pdf"), p_scatter,
       width = 180, height = 180 * 0.8 + 18, units = "mm", dpi = 300, bg = "white")
