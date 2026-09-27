# 5.1 Two-sample MR: ferulic 4-sulfate -> epigenetic age acceleration (EAA), plus reverse MR
#     Forward: metabolite instruments P < 1e-6 (n = 2,392), LD clumping r2 < 0.5, 10 Mb window, F >= 10
#     Reverse: EAA instruments P < 5e-8 (n = 34,710), same LD clumping and F filter
#     LD reference: 500FG genotypes (European ancestry)
#     Run from the repository root: Rscript MR/5.1_ferulic4sulfate_MR_epiclock.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(ieugwasr)
  library(ggplot2)
})

# ---- Parameters ----
P_EXP    <- 1e-6    # metabolite instrument threshold
N_EXP    <- 2392    # sample size reported by the metabolite GWAS (summary statistics source)
P_EAA    <- 5e-8    # EAA instrument threshold (reverse MR)
N_EAA    <- 34710   # sample size reported by the EAA GWAS (McCartney et al. 2021, Genome Biol)
CLUMP_KB <- 10000
CLUMP_R2 <- 0.5
CLUMP_P  <- 1
F_MIN    <- 10

# ---- Paths ----
exposure_gwas <- "data/MR/ferulic_4_sulfate/phenocode-X100005389_Plasma_NHW.tsv.gz"
outcomes <- c(                                     # EAA GWAS summary statistics
  PhenoAge_AA = "data/MR/AA_qtl_filter.RDS",
  Hannum_AA   = "data/MR/Hannum_qtl_filter.RDS",
  IEAA        = "data/MR/IEAA_qtl_filter.RDS",
  GrimAge_AA  = "data/MR/GrimAge_qtl_filter.RDS"
)
ref_bfile   <- "data/genotype/500FG/allchr"        # LD reference (PLINK bfile; IDs "chr:pos:ref:alt;rsID")
plink_bin   <- Sys.which("plink")                # or the full path to a PLINK 1.9 binary
output_dir  <- "results/MR/ferulic_4_sulfate_p1e-06_kb10000_r20.5"
combined_dir <- file.path(output_dir, "MR_combined_results")
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

has_mrpresso <- requireNamespace("MRPRESSO", quietly = TRUE)
has_raps     <- requireNamespace("mr.raps", quietly = TRUE)
has_grip     <- "mr_grip" %in% mr_method_list()$obj

# ---- Helpers ----
safe1 <- function(x, col) {
  if (is.null(x) || !nrow(x) || !(col %in% names(x))) return(NA_real_)
  as.numeric(x[[col]][1])
}
pickrow <- function(df, pat) {
  if (is.null(df) || !nrow(df)) return(NULL)
  head(df[grepl(pat, df$method), , drop = FALSE], 1)
}

add_eaf <- function(d, type) {                     # no allele frequencies in either GWAS
  col <- paste0("eaf.", type)
  if (!col %in% names(d)) d[[col]] <- NA_real_
  d
}
format_eaa <- function(d, type) {                  # EAA GWAS columns: SNP A1 A2 Effect SE P chr pos
  d %>%
    rename(effect_allele = A1, other_allele = A2, beta = Effect, se = SE, pval = P) %>%
    format_data(type = type, snp_col = "SNP", beta_col = "beta", se_col = "se",
                effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                pval_col = "pval", chr_col = "chr", pos_col = "pos") %>%
    add_eaf(type)
}
format_metab <- function(d, type) {                # metabolite GWAS: effect allele = ref
  as.data.frame(d) %>%
    format_data(type = type, snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
                effect_allele_col = "ref", other_allele_col = "alt",
                pval_col = "pval", chr_col = "chrom", pos_col = "pos") %>%
    add_eaf(type)
}

# Harmonise (action = 2), drop removed/ambiguous SNPs, add sample sizes and F, keep F >= F_MIN
harmonise_filter <- function(exp_dat, out_dat, n_exp, n_out) {
  dat <- harmonise_data(exp_dat, out_dat, action = 2)
  if (is.list(dat$remove))    dat$remove    <- unlist(dat$remove)
  if (is.list(dat$ambiguous)) dat$ambiguous <- unlist(dat$ambiguous)
  dat <- filter(dat, remove == FALSE, is.na(ambiguous) | ambiguous == FALSE)
  if (!nrow(dat)) return(NULL)
  dat$samplesize.exposure <- n_exp
  dat$samplesize.outcome  <- n_out
  dat$F_stat <- dat$beta.exposure^2 / dat$se.exposure^2
  attr(dat, "F_mean") <- mean(dat$F_stat, na.rm = TRUE)
  dat_f <- filter(dat, F_stat >= F_MIN)
  attr(dat_f, "F_mean") <- attr(dat, "F_mean")     # mean F before filtering, as reported
  if (nrow(dat_f)) dat_f else NULL
}

# MR methods + primary estimate (Wald ratio for 1 SNP, otherwise IVW random effects) + sensitivity tests
run_mr <- function(dat, methods) {
  ns <- length(unique(dat$SNP))
  if (ns == 1) methods <- "mr_wald_ratio"           # single instrument: Wald ratio
  mr_res <- tryCatch(mr(dat, method_list = methods), error = function(e) NULL)
  if (is.null(mr_res)) return(NULL)
  wald <- pickrow(mr_res, "^Wald ratio$")
  if (ns == 1 && !is.null(wald)) {
    main <- wald; method_main <- "Wald ratio"
  } else {
    main <- pickrow(mr_res, "Inverse variance weighted"); method_main <- "Inverse variance weighted (random effects)"
  }
  het   <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  presso <- NULL
  if (has_mrpresso && ns >= 3) {
    presso <- tryCatch(MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = as.data.frame(dat),
      NbDistribution = 1000, SignifThreshold = 0.05), error = function(e) NULL)
    if (is.null(presso$`MR-PRESSO results`)) presso <- NULL
  }
  het_ivw <- if (!is.null(het)) het[het$method == "Inverse variance weighted", , drop = FALSE] else NULL
  list(
    mr_res = mr_res, het = het, pleio = pleio, presso = presso, dat = dat, ns = ns,
    method_main = method_main,
    beta = safe1(main, "b"), se = safe1(main, "se"), pval = safe1(main, "pval"),
    Q_pval = safe1(het_ivw, "Q_pval"),
    egger_intercept_pval = safe1(pleio, "pval"),
    presso_pval = if (is.null(presso)) NA_real_ else
      tryCatch(presso$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA_real_)
  )
}

save_plots <- function(r, plot_dir, prefix = "") {
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  f <- function(name) file.path(plot_dir, paste0(prefix, name))
  try(ggsave(f("scatter_plot.png"), mr_scatter_plot(r$mr_res, r$dat)[[1]], width = 10, height = 8, dpi = 300))
  try(ggsave(f("forest_plot.png"), mr_forest_plot(mr_singlesnp(r$dat))[[1]], width = 10, height = 12, dpi = 300))
  if (r$ns >= 3) {
    try(ggsave(f("leaveoneout_plot.png"), mr_leaveoneout_plot(mr_leaveoneout(r$dat))[[1]], width = 10, height = 12, dpi = 300))
    try(ggsave(f("funnel_plot.png"), mr_funnel_plot(mr_singlesnp(r$dat))[[1]], width = 8, height = 8, dpi = 300))
  }
}

# LD clumping against the local reference (.bim rewritten with rsIDs as SNP IDs)
tmp_ref <- file.path(tempdir(), "ref_rsid")
fread(paste0(ref_bfile, ".bim"), header = FALSE) %>%
  mutate(V2 = sub(".*;", "", V2)) %>%
  fwrite(paste0(tmp_ref, ".bim"), sep = "\t", col.names = FALSE)
file.symlink(normalizePath(paste0(ref_bfile, c(".bed", ".fam"))), paste0(tmp_ref, c(".bed", ".fam")))

clump_local <- function(exp_dat) {
  cl <- ld_clump(data.frame(rsid = exp_dat$SNP, pval = exp_dat$pval.exposure),
                 clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, clump_p = CLUMP_P,
                 bfile = tmp_ref, plink_bin = plink_bin)
  out <- filter(exp_dat, SNP %in% cl$rsid)
  attr(out, "n_clumped") <- nrow(cl)
  out
}

methods <- c("mr_ivw_mre", "mr_egger_regression", "mr_weighted_median")
if (has_raps) methods <- c(methods, "mr_raps")

# =====================================================================
# Forward MR: ferulic 4-sulfate -> EAA
# =====================================================================
metab_gwas <- fread(exposure_gwas)

# Instruments
exp_sig <- filter(metab_gwas, pval < P_EXP)
exp_dat <- format_metab(exp_sig, "exposure") %>% mutate(exposure = "Ferulic 4 sulfate")
n_sig_before_clump <- nrow(exp_sig)

exp_clumped <- clump_local(exp_dat)
n_after_clump <- attr(exp_clumped, "n_clumped")
message("Instruments: ", n_sig_before_clump, " at P < ", P_EXP, "; ", n_after_clump, " after clumping")

methods_fwd <- if (has_grip) c(methods, "mr_grip") else methods
fwd <- list()
for (outcome_name in names(outcomes)) {
  out_dat <- format_eaa(readRDS(outcomes[[outcome_name]]), "outcome") %>% mutate(outcome = outcome_name)
  dat <- harmonise_filter(exp_clumped, out_dat, N_EXP, N_EAA)
  if (is.null(dat)) { message("No SNPs left for ", outcome_name); next }
  r <- run_mr(dat, methods_fwd)
  if (is.null(r)) next
  message(sprintf("%s: %d SNPs, beta = %.4f, P = %.3g", outcome_name, r$ns, r$beta, r$pval))

  res <- tibble::tibble(
    p_threshold = P_EXP, clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2,
    n_sig_before_clump = n_sig_before_clump, n_after_clump = n_after_clump,
    exposure = "Ferulic 4 sulfate", outcome = outcome_name, n_SNPs_used = r$ns,
    F_mean_exposure_strength = attr(dat, "F_mean"),
    Q_pval_cochrans_Q_IVW = r$Q_pval,
    Egger_intercept_pval_pleiotropy = r$egger_intercept_pval,
    MR_PRESSO_global_pval_outlier_test = r$presso_pval,
    MR_method_primary = r$method_main,
    MR_beta_estimate = r$beta, MR_standard_error = r$se, MR_pval_main = r$pval,
    MR_RAPS_beta = safe1(pickrow(r$mr_res, "RAPS"), "b"), MR_RAPS_pval = safe1(pickrow(r$mr_res, "RAPS"), "pval"),
    MR_GRIP_beta = safe1(pickrow(r$mr_res, "GRIP"), "b"), MR_GRIP_pval = safe1(pickrow(r$mr_res, "GRIP"), "pval")
  )
  fwd[[outcome_name]] <- c(r, list(res = res))

  od <- file.path(output_dir, paste0("MR_output_", outcome_name))
  dir.create(od, showWarnings = FALSE, recursive = TRUE)
  fwrite(res, file.path(od, "MR_sensitivity_results.csv"), sep = ";", dec = ",", quote = FALSE, na = "NA")
  write.csv(r$mr_res, file.path(od, "MR_all_methods.csv"), row.names = FALSE)
  if (!is.null(r$het))    write.csv(r$het, file.path(od, "MR_heterogeneity_detailed.csv"), row.names = FALSE)
  if (!is.null(r$pleio))  write.csv(r$pleio, file.path(od, "MR_pleiotropy_detailed.csv"), row.names = FALSE)
  if (!is.null(r$presso)) saveRDS(r$presso, file.path(od, "MR_PRESSO_full_results.rds"))
  write.csv(dat, file.path(od, "harmonised_data_with_F.csv"), row.names = FALSE)
  save_plots(r, file.path(output_dir, paste0("MR_plots_", outcome_name)))
}

combined_results <- bind_rows(lapply(fwd, `[[`, "res"))
write.csv(combined_results, file.path(combined_dir, "MR_all_outcomes_combined.csv"), row.names = FALSE)
fwrite(combined_results, file.path(combined_dir, "MR_all_outcomes_combined_EUR.csv"),
       sep = ";", dec = ",", quote = FALSE, na = "NA")

# ---- Combined forest plot: per-SNP Wald ratios + IVW summary (85 x 62 mm, 5 pt) ----
outcome_order  <- c("PhenoAge_AA", "Hannum_AA", "GrimAge_AA", "IEAA")
outcome_labels <- c(PhenoAge_AA = "PhenoAge AA", Hannum_AA = "Hannum AA", GrimAge_AA = "GrimAge AA", IEAA = "IEAA")

snp_data <- bind_rows(lapply(names(fwd), function(o) {
  fwd[[o]]$dat %>%
    mutate(beta = beta.outcome / beta.exposure,
           se_mr = abs(beta) * sqrt((se.outcome / beta.outcome)^2 + (se.exposure / beta.exposure)^2),
           CI_lower = beta - 1.96 * se_mr, CI_upper = beta + 1.96 * se_mr,
           label = SNP, outcome = o, type = "SNP", pval = NA_real_) %>%
    select(outcome, label, beta, CI_lower, CI_upper, pval, type)
}))
summary_data <- combined_results %>%
  transmute(outcome, label = paste0("IVW (n=", n_SNPs_used, ")"), beta = MR_beta_estimate,
            CI_lower = beta - 1.96 * MR_standard_error, CI_upper = beta + 1.96 * MR_standard_error,
            pval = MR_pval_main, type = "Summary")

forest_data <- bind_rows(snp_data, summary_data) %>%
  mutate(outcome_group = factor(outcome, levels = outcome_order)) %>%
  arrange(outcome_group, type) %>%
  group_by(outcome_group) %>% mutate(y_in_group = row_number()) %>% ungroup() %>%
  mutate(y_position = y_in_group * 0.6 + (as.numeric(outcome_group) - 1) * 4)

ivw_labels <- filter(forest_data, type == "Summary") %>%
  mutate(text_label = paste0("Beta=", round(beta, 3), "\n", "p=", format(round(pval, 3), nsmall = 3)))

x_max <- max(forest_data$CI_upper, na.rm = TRUE)
x_min <- min(forest_data$CI_lower, na.rm = TRUE)
x_range <- x_max - x_min
x_text_pos <- x_max + x_range * 0.08
x_limit <- x_max + x_range * 0.35
rect_data <- data.frame(xmin = x_text_pos - x_range * 0.02, xmax = x_limit, ymin = -Inf, ymax = Inf,
                        outcome_group = unique(forest_data$outcome_group))

p_forest <- ggplot(forest_data, aes(x = beta, y = y_position)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.25) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = type), height = 0.2, linewidth = 0.35) +
  geom_point(aes(color = type, size = type, shape = type)) +
  geom_rect(data = rect_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", color = NA, inherit.aes = FALSE) +
  geom_label(data = ivw_labels, aes(x = x_text_pos, y = y_position, label = text_label),
             hjust = 0, size = 1.35, color = "black", fill = "white",
             label.padding = unit(0.12, "lines"), label.size = NA, label.r = unit(0, "lines")) +
  scale_color_manual(values = c(SNP = "#2166AC", Summary = "#B2182B")) +
  scale_size_manual(values = c(SNP = 0.8, Summary = 1.4)) +
  scale_shape_manual(values = c(SNP = 16, Summary = 18)) +
  scale_x_continuous(limits = c(x_min, x_limit), expand = c(0.02, 0)) +
  scale_y_continuous(breaks = forest_data$y_position, labels = forest_data$label) +
  coord_cartesian(xlim = c(x_min, x_limit), clip = "off") +
  facet_grid(outcome_group ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(outcome_labels)) +
  labs(x = "MR effect size", y = NULL) +
  theme_minimal(base_size = 5) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray92", linewidth = 0.2),
    axis.text = element_text(size = 5, colour = "black"),
    axis.title = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(margin = margin(t = 1)),
    strip.text.y = element_text(size = 5, colour = "black", margin = margin(0, 0.5, 0, 0.5)),
    strip.background = element_rect(fill = "grey96", color = NA),
    strip.placement = "outside",
    legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = 5, colour = "black"),
    legend.key.size = unit(0.25, "cm"), legend.spacing.x = unit(0.15, "cm"),
    panel.spacing.y = unit(0.4, "lines"), plot.margin = margin(4, 12, 2, 0),
    axis.ticks = element_line(linewidth = 0.2, colour = "black"),
    axis.ticks.length = unit(0.05, "cm")
  )
ggsave(file.path(combined_dir, "combined_forest_plot_detailed_LD_nature.pdf"), p_forest,
       width = 85, height = 62, units = "mm", dpi = 300, bg = "white")

# =====================================================================
# Reverse MR: EAA -> ferulic 4-sulfate
# =====================================================================
metab_out <- format_metab(metab_gwas, "outcome") %>% mutate(outcome = "Ferulic 4 sulfate")
reverse_dir <- file.path(output_dir, "reverse_MR")

rev <- list()
for (clock in names(outcomes)) {
  clock_sig <- filter(readRDS(outcomes[[clock]]), P < P_EAA)
  if (!nrow(clock_sig)) next
  clock_exp <- format_eaa(clock_sig, "exposure") %>%
    mutate(exposure = clock) %>%
    clump_local()
  n_rev_clumped <- attr(clock_exp, "n_clumped")
  if (!nrow(clock_exp)) next
  dat <- harmonise_filter(clock_exp, metab_out, N_EAA, N_EXP)
  if (is.null(dat)) { message("No SNPs left for ", clock, " (reverse)"); next }
  r <- run_mr(dat, methods)
  if (is.null(r)) next
  message(sprintf("Reverse %s: %d SNPs, beta = %.4f, P = %.3g", clock, r$ns, r$beta, r$pval))

  res <- tibble::tibble(
    direction = "Reverse", exposure = clock, outcome = "Ferulic 4 sulfate",
    n_sig_before_clump = nrow(clock_sig), n_after_clump = n_rev_clumped, n_SNPs_used = r$ns,
    F_mean = attr(dat, "F_mean"), Q_pval_heterogeneity = r$Q_pval,
    Egger_intercept_pval = r$egger_intercept_pval, MR_PRESSO_pval = r$presso_pval,
    MR_method = r$method_main, MR_beta = r$beta, MR_se = r$se, MR_pval = r$pval,
    MR_Egger_beta = safe1(pickrow(r$mr_res, "MR Egger"), "b"), MR_Egger_pval = safe1(pickrow(r$mr_res, "MR Egger"), "pval"),
    MR_WM_beta = safe1(pickrow(r$mr_res, "Weighted median"), "b"), MR_WM_pval = safe1(pickrow(r$mr_res, "Weighted median"), "pval"),
    MR_RAPS_beta = safe1(pickrow(r$mr_res, "RAPS"), "b"), MR_RAPS_pval = safe1(pickrow(r$mr_res, "RAPS"), "pval")
  )
  rev[[clock]] <- res

  od <- file.path(reverse_dir, paste0("reverse_MR_", clock))
  dir.create(od, showWarnings = FALSE, recursive = TRUE)
  write.csv(res, file.path(od, "reverse_MR_results.csv"), row.names = FALSE)
  write.csv(r$mr_res, file.path(od, "reverse_MR_all_methods.csv"), row.names = FALSE)
  write.csv(dat, file.path(od, "reverse_harmonised_data.csv"), row.names = FALSE)
  try(ggsave(file.path(od, "reverse_scatter_plot.png"), mr_scatter_plot(r$mr_res, dat)[[1]], width = 10, height = 8, dpi = 300))
  try(ggsave(file.path(od, "reverse_forest_plot.png"), mr_forest_plot(mr_singlesnp(dat))[[1]], width = 10, height = 10, dpi = 300))
}
if (length(rev)) write.csv(bind_rows(rev), file.path(reverse_dir, "reverse_MR_all_outcomes.csv"), row.names = FALSE)
