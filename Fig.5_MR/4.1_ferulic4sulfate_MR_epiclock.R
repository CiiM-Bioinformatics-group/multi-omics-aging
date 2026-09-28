# 5.1 Bidirectional two-sample MR: ferulic 4-sulfate <-> epigenetic age acceleration (EAA)


library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ieugwasr)
library(ggplot2)

P_EXP    <- 1e-6    # metabolite instruments
N_EXP    <- 2392    # metabolite GWAS sample size
P_EAA    <- 5e-8    # EAA instruments (reverse MR)
N_EAA    <- 34710   # EAA GWAS sample size
CLUMP_KB <- 10000
CLUMP_R2 <- 0.5
CLUMP_P  <- 1
F_MIN    <- 10

ref_bfile  <- "data/MR/ref_1kg/EUR"                  # 1000G phase 3 EUR (rsID-based)
plink_bin  <- Sys.which("plink")                     # or the full path to a PLINK 1.9 binary
metab_file <- "data/MR/ferulic_4_sulfate/phenocode-X100005389_Plasma_NHW.tsv.gz"   # https://ontime.wustl.edu/pheno/X100005389_Plasma_NHW
outcomes <- c(                                       # EAA GWAS (McCartney et al. 2021)
  PhenoAge_AA = "data/MR/AA_qtl_filter.RDS",
  Hannum_AA   = "data/MR/Hannum_qtl_filter.RDS",
  IEAA        = "data/MR/IEAA_qtl_filter.RDS",
  GrimAge_AA  = "data/MR/GrimAge_qtl_filter.RDS"
)
output_dir <- "results/MR"
stopifnot(file.exists(paste0(ref_bfile, c(".bed", ".bim", ".fam"))))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

has_mrpresso <- requireNamespace("MRPRESSO", quietly = TRUE)
has_raps     <- requireNamespace("mr.raps", quietly = TRUE)
has_grip     <- "mr_grip" %in% mr_method_list()$obj
methods_rev  <- c("mr_ivw_mre", "mr_egger_regression", "mr_weighted_median", if (has_raps) "mr_raps")
methods_fwd  <- c(methods_rev, if (has_grip) "mr_grip")

# ---- Helpers ----
safe1 <- function(x, col) {
  if (is.null(x) || !nrow(x) || !(col %in% names(x))) return(NA_real_)
  as.numeric(x[[col]][1])
}
pickrow <- function(df, pat) {
  if (is.null(df) || !nrow(df)) return(NULL)
  head(df[grepl(pat, df$method), , drop = FALSE], 1)
}
add_eaf <- function(d, type) {
  col <- paste0("eaf.", type)
  if (!col %in% names(d)) d[[col]] <- NA_real_
  d
}
format_eaa <- function(d, type) {        # EAA GWAS: SNP A1 A2 Effect SE P chr pos
  d %>%
    rename(effect_allele = A1, other_allele = A2, beta = Effect, se = SE, pval = P) %>%
    format_data(type = type, snp_col = "SNP", beta_col = "beta", se_col = "se",
                effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                pval_col = "pval", chr_col = "chr", pos_col = "pos") %>%
    add_eaf(type)
}
format_metab <- function(d, type) {      # metabolite GWAS: effect allele = ref
  as.data.frame(d) %>%
    format_data(type = type, snp_col = "rsids", beta_col = "beta", se_col = "sebeta",
                effect_allele_col = "ref", other_allele_col = "alt",
                pval_col = "pval", chr_col = "chrom", pos_col = "pos") %>%
    add_eaf(type)
}
clump_1kg <- function(exp_dat) {
  cl <- ld_clump(data.frame(rsid = exp_dat$SNP, pval = exp_dat$pval.exposure),
                 clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, clump_p = CLUMP_P,
                 bfile = ref_bfile, plink_bin = plink_bin)
  filter(exp_dat, SNP %in% cl$rsid)
}

# Harmonise -> F filter -> MR (Wald ratio if 1 SNP, else IVW random effects) -> sensitivity
run_mr <- function(exp_clumped, out_dat, n_exp, n_out, methods) {
  dat <- harmonise_data(exp_clumped, out_dat, action = 2)
  if (is.list(dat$remove))    dat$remove    <- unlist(dat$remove)
  if (is.list(dat$ambiguous)) dat$ambiguous <- unlist(dat$ambiguous)
  dat <- filter(dat, remove == FALSE, is.na(ambiguous) | ambiguous == FALSE)
  if (!nrow(dat)) return(NULL)
  dat$samplesize.exposure <- n_exp
  dat$samplesize.outcome  <- n_out
  dat$F_stat <- dat$beta.exposure^2 / dat$se.exposure^2
  F_mean <- mean(dat$F_stat, na.rm = TRUE)
  dat <- filter(dat, F_stat >= F_MIN)
  ns <- length(unique(dat$SNP))
  if (ns == 0) return(NULL)
  if (ns == 1) methods <- "mr_wald_ratio"

  mr_res <- mr(dat, method_list = methods)
  main <- if (ns == 1) pickrow(mr_res, "^Wald ratio$") else pickrow(mr_res, "Inverse variance weighted")
  het   <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  presso <- NULL
  if (has_mrpresso && ns >= 3) {
    set.seed(123)
    presso <- tryCatch(MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = as.data.frame(dat),
      NbDistribution = 1000, SignifThreshold = 0.05), error = function(e) NULL)
  }
  presso_p <- tryCatch(as.numeric(presso$`MR-PRESSO results`$`Global Test`$Pvalue), error = function(e) NA_real_)
  if (!length(presso_p)) presso_p <- NA_real_
  het_ivw <- if (!is.null(het)) het[het$method == "Inverse variance weighted", , drop = FALSE] else NULL
  summary <- tibble::tibble(
    exposure = dat$exposure[1], outcome = dat$outcome[1],
    n_SNPs_used = ns, F_mean = F_mean,
    Q_pval_IVW = safe1(het_ivw, "Q_pval"),
    Egger_intercept_pval = safe1(pleio, "pval"),
    MR_PRESSO_global_pval = presso_p,
    MR_method = if (ns == 1) "Wald ratio" else "Inverse variance weighted (random effects)",
    MR_beta = safe1(main, "b"), MR_se = safe1(main, "se"), MR_pval = safe1(main, "pval"),
    Egger_beta = safe1(pickrow(mr_res, "MR Egger"), "b"), Egger_pval = safe1(pickrow(mr_res, "MR Egger"), "pval"),
    WM_beta = safe1(pickrow(mr_res, "Weighted median"), "b"), WM_pval = safe1(pickrow(mr_res, "Weighted median"), "pval"),
    RAPS_beta = safe1(pickrow(mr_res, "RAPS"), "b"), RAPS_pval = safe1(pickrow(mr_res, "RAPS"), "pval"),
    GRIP_beta = safe1(pickrow(mr_res, "GRIP"), "b"), GRIP_pval = safe1(pickrow(mr_res, "GRIP"), "pval")
  )
  list(summary = summary, mr_res = mr_res, het = het, pleio = pleio, presso = presso, dat = dat, ns = ns)
}

save_mr <- function(r, od) {
  dir.create(od, showWarnings = FALSE, recursive = TRUE)
  write.csv(r$summary, file.path(od, "MR_summary.csv"), row.names = FALSE)
  write.csv(r$mr_res,  file.path(od, "MR_all_methods.csv"), row.names = FALSE)
  write.csv(r$dat,     file.path(od, "harmonised_data_with_F.csv"), row.names = FALSE)
  if (!is.null(r$het))    write.csv(r$het,   file.path(od, "MR_heterogeneity_detailed.csv"), row.names = FALSE)
  if (!is.null(r$pleio))  write.csv(r$pleio, file.path(od, "MR_pleiotropy_detailed.csv"), row.names = FALSE)
  if (!is.null(r$presso)) saveRDS(r$presso,  file.path(od, "MR_PRESSO_full_results.rds"))
  ss <- tryCatch(mr_singlesnp(r$dat), error = function(e) NULL)
  try(ggsave(file.path(od, "scatter_plot.png"), mr_scatter_plot(r$mr_res, r$dat)[[1]], width = 10, height = 8, dpi = 300))
  if (!is.null(ss)) try(ggsave(file.path(od, "forest_plot.png"), mr_forest_plot(ss)[[1]], width = 10, height = 12, dpi = 300))
  if (r$ns >= 3) {
    try(ggsave(file.path(od, "leaveoneout_plot.png"), mr_leaveoneout_plot(mr_leaveoneout(r$dat))[[1]], width = 10, height = 12, dpi = 300))
    if (!is.null(ss)) try(ggsave(file.path(od, "funnel_plot.png"), mr_funnel_plot(ss)[[1]], width = 8, height = 8, dpi = 300))
  }
}

# ---- Load GWAS ----
metab_gwas <- fread(metab_file)
metab_out  <- format_metab(metab_gwas, "outcome") %>% mutate(outcome = "Ferulic 4 sulfate")
eaa <- lapply(outcomes, readRDS)

# ---- Forward MR: ferulic 4-sulfate -> EAA ----
exp_sig <- filter(metab_gwas, pval < P_EXP)
exp_dat <- format_metab(exp_sig, "exposure") %>% mutate(exposure = "Ferulic 4 sulfate")
exp_clumped <- clump_1kg(exp_dat)
message("Forward instruments: ", nrow(exp_sig), " at P < ", P_EXP, "; ", nrow(exp_clumped), " after clumping")

forward <- list()
for (o in names(outcomes)) {
  out_dat <- format_eaa(eaa[[o]], "outcome") %>% mutate(outcome = o)
  r <- run_mr(exp_clumped, out_dat, N_EXP, N_EAA, methods_fwd)
  if (is.null(r)) { message(o, ": no SNPs left"); next }
  r$summary <- mutate(r$summary, direction = "Forward", n_sig_before_clump = nrow(exp_sig),
                      n_after_clump = nrow(exp_clumped), .before = 1)
  save_mr(r, file.path(output_dir, paste0("MR_output_", o)))
  forward[[o]] <- r$summary
}
forward_df <- bind_rows(forward)

# ---- Reverse MR: EAA -> ferulic 4-sulfate ----
reverse <- list()
for (o in names(outcomes)) {
  clock_sig <- filter(eaa[[o]], P < P_EAA)
  clock_exp <- format_eaa(clock_sig, "exposure") %>% mutate(exposure = o)
  clock_clumped <- clump_1kg(clock_exp)
  message(o, " reverse instruments: ", nrow(clock_sig), " at P < ", P_EAA, "; ", nrow(clock_clumped), " after clumping")
  if (!nrow(clock_clumped)) next
  r <- run_mr(clock_clumped, metab_out, N_EAA, N_EXP, methods_rev)
  if (is.null(r)) { message(o, ": no SNPs left"); next }
  r$summary <- mutate(r$summary, direction = "Reverse", n_sig_before_clump = nrow(clock_sig),
                      n_after_clump = nrow(clock_clumped), .before = 1)
  save_mr(r, file.path(output_dir, "reverse_MR", paste0("reverse_MR_", o)))
  reverse[[o]] <- r$summary
}
reverse_df <- bind_rows(reverse)

# ---- Combined table ----
bidirectional <- bind_rows(forward_df, reverse_df)
combined_dir <- file.path(output_dir, "MR_combined_results")
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(bidirectional, file.path(combined_dir, "MR_bidirectional.csv"), row.names = FALSE)
bidirectional %>%
  select(direction, exposure, outcome, n_after_clump, n_SNPs_used, MR_method, MR_beta, MR_se, MR_pval,
         Q_pval_IVW, Egger_intercept_pval, MR_PRESSO_global_pval) %>%
  as.data.frame() %>% print()
