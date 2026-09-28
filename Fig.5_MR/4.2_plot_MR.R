# 5.2 Forward MR figures: IVW forest plot and combined scatter plot (outputs of 5.1)
#     Run from the repository root: Rscript MR/5.2_plot_MR.R
library(data.table)
library(dplyr)
library(ggplot2)

mr_dir <- "results/MR"
combined_dir <- file.path(mr_dir, "MR_combined_results")
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

outcome_order  <- c("PhenoAge_AA", "Hannum_AA", "GrimAge_AA", "IEAA")
outcome_labels <- c(PhenoAge_AA = "PhenoAge AA", Hannum_AA = "Hannum AA", GrimAge_AA = "GrimAge AA", IEAA = "IEAA")

# ---- Load forward results ----
present <- outcome_order[dir.exists(file.path(mr_dir, paste0("MR_output_", outcome_order)))]
load_csv <- function(o, f) read.csv(file.path(mr_dir, paste0("MR_output_", o), f), stringsAsFactors = FALSE)
harm    <- setNames(lapply(present, load_csv, f = "harmonised_data_with_F.csv"), present)
methods <- setNames(lapply(present, load_csv, f = "MR_all_methods.csv"), present)
summ    <- bind_rows(lapply(present, load_csv, f = "MR_summary.csv"))

# =====================================================================
# 1. Forest plot: IVW estimate per outcome with nSNPs, Beta (95% CI), P
# =====================================================================
fp <- summ %>%
  mutate(lo = MR_beta - 1.96 * MR_se, hi = MR_beta + 1.96 * MR_se,
         outcome = factor(outcome, levels = outcome_order),
         beta_ci = sprintf("%.3f (%.3f, %.3f)", MR_beta, lo, hi),
         p_lab = ifelse(MR_pval < 0.001, formatC(MR_pval, format = "e", digits = 2), sprintf("%.3f", MR_pval))) %>%
  arrange(outcome) %>%
  mutate(y = row_number())                         # PhenoAge at the bottom, IEAA on top

x_lo <- min(-1, floor(min(fp$lo) * 2) / 2)          # axis in 0.5 steps, always including 0
x_hi <- max(0, ceiling(max(fp$hi) * 2) / 2)
span <- x_hi - x_lo
fig_w <- 90                                         # figure width (mm)
fig_h <- 50                                         # figure height (mm)
m_left <- 72; m_right <- 104                        # margins (pt) holding the text columns
panel_pt <- fig_w / 25.4 * 72 - m_left - m_right    # panel width (pt)
pt2x <- function(pt) pt / panel_pt * span           # convert pt offsets to x units
x_out  <- x_lo - pt2x(70)                           # text column positions (outside the panel)
x_nsnp <- x_lo - pt2x(24)
x_beta <- x_hi + pt2x(5)
x_p    <- x_hi + pt2x(74)
y_head <- nrow(fp) + 0.8
txt <- 1.8                                          # text size (mm), ~5 pt

p_forest <- ggplot(fp, aes(x = MR_beta, y = y)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.35, colour = "black") +
  geom_point(size = 1.4, colour = "#3C78C3") +
  geom_text(aes(x = x_out,  label = outcome), hjust = 0, size = txt) +
  geom_text(aes(x = x_nsnp, label = n_SNPs_used), hjust = 0, size = txt) +
  geom_text(aes(x = x_beta, label = beta_ci), hjust = 0, size = txt) +
  geom_text(aes(x = x_p,    label = p_lab), hjust = 0, size = txt) +
  annotate("text", x = c(x_out, x_nsnp, x_beta, x_p), y = y_head,
           label = c("Outcome", "nSNPs", "Beta (95% CI)", "p-value"),
           hjust = 0, size = txt, fontface = "bold") +
  scale_x_continuous(breaks = seq(x_lo, x_hi, by = 0.5), labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(limits = c(0.5, y_head + 0.3), expand = c(0, 0)) +
  coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
  labs(x = "MR effect (IVW; beta scale)", y = NULL) +
  theme_classic(base_size = 5) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.25), axis.ticks.x = element_line(linewidth = 0.25),
    axis.text.x = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size = 5, colour = "black"),
    plot.margin = margin(4, m_right, 4, m_left)
  )
ggsave(file.path(combined_dir, "forest_plot_IVW.pdf"), p_forest,
       width = fig_w, height = fig_h, units = "mm", bg = "white")
print(p_forest)

# =====================================================================
# 2. Scatter plot: SNP effects with MR method lines (180 mm wide)
# =====================================================================
tsmr_cols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
               "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")

snp_list <- list(); line_list <- list()
for (o in present) {
  dat <- harm[[o]]
  mr_res <- methods[[o]]
  if ("mr_keep" %in% names(dat)) dat <- dat[dat$mr_keep %in% c(TRUE, "TRUE", 1), , drop = FALSE]
  flip <- dat$beta.exposure < 0                       # orient to positive exposure effect
  dat$beta.exposure[flip] <- -dat$beta.exposure[flip]
  dat$beta.outcome[flip]  <- -dat$beta.outcome[flip]
  mr_res$a <- 0                                       # intercept 0, except MR-Egger
  egger <- mr_res$method %in% c("MR Egger", "MR Egger (bootstrap)")
  if (any(egger)) {
    fit <- lm(beta.outcome ~ beta.exposure, data = dat, weights = 1 / dat$se.outcome^2)
    mr_res$a[egger] <- unname(coef(fit)[1])
  }
  snp_list[[o]]  <- mutate(dat, outcome_name = o)
  line_list[[o]] <- mutate(mr_res, outcome_name = o)
}
lab <- function(x) factor(outcome_labels[x], levels = unname(outcome_labels[present]))
snp_df  <- bind_rows(snp_list)  %>% mutate(outcome_lab = lab(outcome_name))
line_df <- bind_rows(line_list) %>% mutate(outcome_lab = lab(outcome_name))
method_levels <- unique(unlist(lapply(line_list, `[[`, "method")))
line_df$method <- factor(line_df$method, levels = method_levels)
pal <- setNames(tsmr_cols[seq_along(method_levels)], method_levels)

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
print(p_scatter)
