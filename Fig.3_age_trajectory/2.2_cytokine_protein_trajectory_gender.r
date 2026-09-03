library(dplyr)
library(tidyr)
library(ggplot2)
library(segmented)

load("/vol/projects/yzhang/500FG_aging/output/14_AA/immune_aging_trajectory_calculated_gender.Rdata")

merge_data <- read.csv("/vol/projects/yzhang/500FG_aging/output/14_AA/500fg_heatmap_merge_data_gender.csv",row.names=1)
head(merge_data)
dim(merge_data)
sum(is.na(merge_data))

cn <- colnames(merge_data)
cytokine <- merge_data[, match("IFNy_C.conidiaHK_WB_48h", cn):match("IL22_Bacteroides_PBMC_7days", cn)]
protein  <- merge_data[, (match("IL22_Bacteroides_PBMC_7days", cn) + 1):ncol(merge_data)]
head(cytokine)
dim(cytokine)
head(protein)
dim(protein)

###line should make by itself
Age <- merge_data$Age
# Z-score per column, then row mean (handles different scales across markers)
cyt_score <- rowMeans(scale(cytokine), na.rm = TRUE)
prot_score <- rowMeans(scale(protein), na.rm = TRUE)
df1 <- data.frame(Age, cyt_score, prot_score)
df1 <- df1[complete.cases(df1), ]

# LOESS fit (span: larger = smoother, default 0.75)
fit_cyt <- loess(cyt_score ~ Age, data = df1, span = 0.75)
fit_prot <- loess(prot_score ~ Age, data = df1, span = 0.75)
age_seq <- seq(min(df1$Age), max(df1$Age), length.out = 200)
pred_cyt <- predict(fit_cyt, newdata = data.frame(Age = age_seq))
pred_prot <- predict(fit_prot, newdata = data.frame(Age = age_seq))
# Normalize fitted values to 0-1 for schematic y-axis
pred_cyt_norm <- (pred_cyt - min(pred_cyt)) / diff(range(pred_cyt))
pred_prot_norm <- (pred_prot - min(pred_prot)) / diff(range(pred_prot))

tau_cyt  <- segmented(lm(cyt_score ~ Age, data = df1), seg.Z = ~ Age)$psi[1, 2]
tau_prot <- segmented(lm(prot_score ~ Age, data = df1), seg.Z = ~ Age)$psi[1, 2]

df_plot <- rbind(
  data.frame(Age = age_seq, value = pred_cyt_norm, type = "Immunosenescence (cytokine decline)"),
  data.frame(Age = age_seq, value = pred_prot_norm, type = "Inflammaging (proteomic shift)")
)

p1 <- ggplot(df_plot, aes(x = Age, y = value, color = type)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(tau_cyt, tau_prot), linetype = 2, color = "gray40") +
 annotate("text", x = tau_cyt+1, y = 0.2, label = paste0("Cytokine trajectory transition (", round(tau_cyt, 0), " y)"), size = 1.8, color = "black",angle = 90,  hjust = -0.1, vjust =0.5) +
  annotate("text", x = tau_prot+1, y = 0.2, label = paste0("Proteomic trajectory transition (", round(tau_prot, 0), " y)"), size = 1.8, color = "black", angle = 90, hjust = -0.1, vjust = 0.5) +
  scale_color_manual(values = c("Inflammaging (proteomic shift)" = "#56B4E9", "Immunosenescence (cytokine decline)" = "#E69F00")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Age (years)", y = "Relative trajectory (schematic)", color = NULL) +
  theme_minimal(base_size = 5) +
  theme(
   legend.position = "bottom",
    text = element_text(size = 5),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
       plot.margin = margin(3, 3, 3, 3, "mm") 
  ) +
  guides(color = guide_legend(ncol = 2))  # vertical legend (stacked)
ggsave("/vol/projects/yzhang/500FG_aging/output/14_AA/immune_aging_trajectory_calculated_gender.pdf", p1, width = 80, height = 70, units = "mm")
p1


save.image("/vol/projects/yzhang/500FG_aging/output/14_AA/immune_aging_trajectory_calculated_gender.Rdata")

### layer-specific DE-SWAN transition point
Age <- merge_data$Age

# Ensure numeric matrix for both layers
cytokine <- as.data.frame(lapply(cytokine, as.numeric))
protein  <- as.data.frame(lapply(protein,  as.numeric))

# Layer-specific DE-SWAN proportions per sample
# If values are already proportions (or 0/1 DE flags), rowMeans is appropriate.
# If your values are not proportions, set use_binarize = TRUE and tune de_threshold.
use_binarize <- FALSE
de_threshold <- 0

cyt_prop <- if (use_binarize) {
  rowMeans(abs(as.matrix(cytokine)) > de_threshold, na.rm = TRUE)
} else {
  rowMeans(as.matrix(cytokine), na.rm = TRUE)
}

prot_prop <- if (use_binarize) {
  rowMeans(abs(as.matrix(protein)) > de_threshold, na.rm = TRUE)
} else {
  rowMeans(as.matrix(protein), na.rm = TRUE)
}

df1 <- data.frame(Age, cyt_prop, prot_prop)
df1 <- df1[complete.cases(df1), ]

# Helper to extract segmented breakpoint robustly
get_tau <- function(y, x) {
  fit_lm <- lm(y ~ x)
  fit_seg <- segmented(fit_lm, seg.Z = ~x, psi = median(x, na.rm = TRUE))
  if ("Est." %in% colnames(fit_seg$psi)) {
    as.numeric(fit_seg$psi[1, "Est."])
  } else {
    as.numeric(fit_seg$psi[1, 2])
  }
}

# LOESS trajectories
age_seq <- seq(min(df1$Age), max(df1$Age), length.out = 200)
fit_cyt  <- loess(cyt_prop  ~ Age, data = df1, span = 0.75)
fit_prot <- loess(prot_prop ~ Age, data = df1, span = 0.75)

pred_cyt  <- predict(fit_cyt,  newdata = data.frame(Age = age_seq))
pred_prot <- predict(fit_prot, newdata = data.frame(Age = age_seq))

# Normalize to 0-1 for schematic y-axis
norm01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  d <- diff(r)
  if (!is.finite(d) || d == 0) return(rep(0.5, length(x)))
  (x - r[1]) / d
}
pred_cyt_norm  <- norm01(pred_cyt)
pred_prot_norm <- norm01(pred_prot)

# Transition timepoints (layer-specific)
tau_cyt  <- get_tau(df1$cyt_prop,  df1$Age)
tau_prot <- get_tau(df1$prot_prop, df1$Age)

cat("Cytokine layer transition timepoint:", round(tau_cyt, 2), "years\n")
cat("Protein layer transition timepoint :", round(tau_prot, 2), "years\n")

df_plot <- rbind(
  data.frame(Age = age_seq, value = pred_cyt_norm,  type = "Cytokine layer"),
  data.frame(Age = age_seq, value = pred_prot_norm, type = "Protein layer")
)

p1 <- ggplot(df_plot, aes(x = Age, y = value, color = type)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(tau_cyt, tau_prot), linetype = 2, color = "gray40") +
  annotate("text", x = tau_cyt + 1, y = 0.20,
           label = paste0("Cytokine transition (", round(tau_cyt, 0), " y)"),
           size = 1.8, color = "black", angle = 90, hjust = -0.1, vjust = 0.5) +
  annotate("text", x = tau_prot + 1, y = 0.80,
           label = paste0("Protein transition (", round(tau_prot, 0), " y)"),
           size = 1.8, color = "black", angle = 90, hjust = -0.1, vjust = 0.5) +
  scale_color_manual(values = c("Cytokine layer" = "#E69F00", "Protein layer" = "#56B4E9")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Age (years)", y = "Relative trajectory (schematic)", color = NULL) +
  theme_minimal(base_size = 5) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 5),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    plot.margin = margin(3, 3, 3, 3, "mm")
  ) +
  guides(color = guide_legend(ncol = 2))

ggsave(
  "/vol/projects/yzhang/500FG_aging/output/14_AA/immune_aging_trajectory_DE-SWAN.pdf",
  p1, width = 80, height = 70, units = "mm"
)

p1

Age <- merge_data$Age

# 1) Prepare numeric matrices
cytokine <- as.data.frame(lapply(cytokine, as.numeric))
protein  <- as.data.frame(lapply(protein,  as.numeric))

# 2) Z-score each feature to match heatmap-style pattern extraction
cyt_z  <- scale(as.matrix(cytokine))
prot_z <- scale(as.matrix(protein))

# 3) Age binning (5-year bins)
bin_width <- 5
age_min <- floor(min(Age, na.rm = TRUE) / bin_width) * bin_width
age_max <- ceiling(max(Age, na.rm = TRUE) / bin_width) * bin_width
breaks <- seq(age_min, age_max, by = bin_width)

df_age <- data.frame(Age = Age) |>
  mutate(bin = cut(Age, breaks = breaks, include.lowest = TRUE, right = FALSE))

# 4) Build bin-level trajectory per layer
get_bin_traj <- function(zmat, df_age) {
  split_idx <- split(seq_len(nrow(df_age)), df_age$bin)
  split_idx <- split_idx[sapply(split_idx, length) >= 5]  # minimum bin size

  centers <- sapply(names(split_idx), function(b) {
    lr <- as.numeric(gsub("\\[|\\)|\\]|\\(", "", unlist(strsplit(b, ","))))
    mean(lr)
  })

  n_bin <- sapply(split_idx, length)

  # Mean expression per feature in each age bin
  # Then summarize to one layer score per bin (mean across features)
  layer_score <- sapply(split_idx, function(ix) {
    feat_means <- colMeans(zmat[ix, , drop = FALSE], na.rm = TRUE)
    mean(feat_means, na.rm = TRUE)
  })

  data.frame(Age = as.numeric(centers), score = as.numeric(layer_score), n = as.numeric(n_bin))
}

traj_cyt  <- get_bin_traj(cyt_z, df_age)
traj_prot <- get_bin_traj(prot_z, df_age)

# 5) Breakpoint on bin-level trajectory (weighted by bin sample size)
get_tau <- function(traj, psi_init) {
  lm0 <- lm(score ~ Age, data = traj, weights = n)
  seg0 <- segmented(lm0, seg.Z = ~Age, psi = psi_init)
  if ("Est." %in% colnames(seg0$psi)) as.numeric(seg0$psi[1, "Est."]) else as.numeric(seg0$psi[1, 2])
}

# Use your visual prior as initial values
tau_cyt  <- get_tau(traj_cyt,  psi_init = 45)
tau_prot <- get_tau(traj_prot, psi_init = 60)

cat("Cytokine transition:", round(tau_cyt, 2), "years\n")
cat("Protein transition :", round(tau_prot, 2), "years\n")

df_plot <- rbind(
  data.frame(Age = age_seq, value = pred_cyt_norm,  type = "Cytokine responses(Immunosenescence)"),
  data.frame(Age = age_seq, value = pred_prot_norm, type = "Protein(Inflammaging)")
)

p1 <- ggplot(df_plot, aes(x = Age, y = value, color = type)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(tau_cyt, tau_prot), linetype = 2, color = "gray40") +
  annotate("text", x = tau_cyt + 1, y = 0.20,
           label = paste0("Cytokine transition (", round(tau_cyt, 1), " y)"),
           size = 1.8, color = "black", angle = 90, hjust = -0.1, vjust = 0.5) +
  annotate("text", x = tau_prot + 1, y = 0.20,
           label = paste0("Protein transition (", round(tau_prot, 1), " y)"),
           size = 1.8, color = "black", angle = 90, hjust = -0.1, vjust = 0.5) +
  scale_color_manual(values = c("Cytokine responses(Immunosenescence)" = "#E69F00", "Protein(Inflammaging)" = "#56B4E9")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Age (years)", y = "Relative trajectory (schematic)", color = NULL) +
  theme_minimal(base_size = 5) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 5),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    plot.margin = margin(3, 3, 3, 3, "mm")
  ) +
  guides(color = guide_legend(ncol = 2))

ggsave(
  "/vol/projects/yzhang/500FG_aging/output/14_AA/immune_aging_trajectory_Age_segmented.pdf",
  p1, width = 80, height = 70, units = "mm"
)

p1

print(tau_cyt)
print(tau_prot)

Age <- merge_data$Age
# 1) Prepare numeric matrices
cytokine <- as.data.frame(lapply(cytokine, as.numeric))
protein  <- as.data.frame(lapply(protein,  as.numeric))
# 2) Sample-level proportions for plotting the two smooth lines
cyt_prop  <- rowMeans(as.matrix(cytokine), na.rm = TRUE)
prot_prop <- rowMeans(as.matrix(protein),  na.rm = TRUE)
df1 <- data.frame(Age = Age, cyt_prop = cyt_prop, prot_prop = prot_prop)
df1 <- df1[complete.cases(df1), ]
# 3) Build age-binned trajectory for breakpoint estimation
get_bin_traj <- function(age, y, bin_width = 5, min_n = 5) {
  age_min <- floor(min(age, na.rm = TRUE) / bin_width) * bin_width
  age_max <- ceiling(max(age, na.rm = TRUE) / bin_width) * bin_width
  breaks <- seq(age_min, age_max, by = bin_width)
  d <- data.frame(Age = age, score = y)
  d <- d[complete.cases(d), ]
  d$bin <- cut(d$Age, breaks = breaks, include.lowest = TRUE, right = FALSE)
  sp <- split(d, d$bin)
  sp <- sp[sapply(sp, nrow) >= min_n]
  if (length(sp) < 4) return(data.frame(Age = numeric(), score = numeric(), n = numeric()))
  centers <- sapply(names(sp), function(b) {
    lr <- as.numeric(gsub("\\[|\\)|\\]|\\(", "", unlist(strsplit(b, ","))))
    mean(lr)
  })
  data.frame(
    Age = as.numeric(centers),
    score = sapply(sp, function(x) mean(x$score, na.rm = TRUE)),
    n = sapply(sp, nrow)
  )
}
traj_cyt  <- get_bin_traj(df1$Age, df1$cyt_prop,  bin_width = 5, min_n = 5)
traj_prot <- get_bin_traj(df1$Age, df1$prot_prop, bin_width = 5, min_n = 5)
# 4) Automatic initialization for segmented breakpoint
get_tau_auto <- function(traj) {
  traj <- traj[is.finite(traj$Age) & is.finite(traj$score) & is.finite(traj$n), ]
  if (nrow(traj) < 6) return(NA_real_)
  lm0 <- lm(score ~ Age, data = traj, weights = n)
  psi_grid <- unique(as.numeric(quantile(
    traj$Age, probs = seq(0.2, 0.8, by = 0.1), na.rm = TRUE
  )))
  fits <- lapply(psi_grid, function(p0) {
    tryCatch(segmented(lm0, seg.Z = ~Age, psi = p0), error = function(e) NULL)
  })
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) return(NA_real_)
  best_fit <- fits[[which.min(sapply(fits, AIC))]]
  if ("Est." %in% colnames(best_fit$psi)) {
    as.numeric(best_fit$psi[1, "Est."])
  } else {
    as.numeric(best_fit$psi[1, 2])
  }
}
tau_cyt  <- get_tau_auto(traj_cyt)
tau_prot <- get_tau_auto(traj_prot)
cat("Cytokine layer transition timepoint (binned):", round(tau_cyt, 2), "years\n")
cat("Protein layer transition timepoint (binned) :", round(tau_prot, 2), "years\n")
# 5) Keep original style: two LOESS lines from sample-level data
age_seq <- seq(min(df1$Age), max(df1$Age), length.out = 200)
fit_cyt  <- loess(cyt_prop ~ Age,  data = df1, span = 0.75)
fit_prot <- loess(prot_prop ~ Age, data = df1, span = 0.75)
pred_cyt  <- predict(fit_cyt,  newdata = data.frame(Age = age_seq))
pred_prot <- predict(fit_prot, newdata = data.frame(Age = age_seq))
norm01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  d <- diff(r)
  if (!is.finite(d) || d == 0) return(rep(0.5, length(x)))
  (x - r[1]) / d
}
pred_cyt_norm  <- norm01(pred_cyt)
pred_prot_norm <- norm01(pred_prot)
df_plot <- rbind(
  data.frame(Age = age_seq, value = pred_cyt_norm,  type = "Cytokine layer"),
  data.frame(Age = age_seq, value = pred_prot_norm, type = "Protein layer")
)
p1 <- ggplot(df_plot, aes(x = Age, y = value, color = type)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(tau_cyt, tau_prot), linetype = 2, color = "gray40") +
  annotate("text", x = tau_cyt + 1, y = 0.20,
           label = paste0("Cytokine transition (", sprintf("%.1f", tau_cyt), " y)"),
           size = 1.8, color = "black", angle = 90, hjust = 0, vjust = 0.5) +
  annotate("text", x = tau_prot + 1, y = 0.70,
           label = paste0("Protein transition (", sprintf("%.1f", tau_prot), " y)"),
           size = 1.8, color = "black", angle = 90, hjust = 0, vjust = 0.5) +
  scale_color_manual(values = c("Cytokine layer" = "#E69F00", "Protein layer" = "#56B4E9")) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.08))) +
  coord_cartesian(clip = "off") +
  labs(x = "Age (years)", y = "Relative trajectory (schematic)", color = NULL) +
  theme_minimal(base_size = 5) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 5),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    plot.margin = margin(3, 8, 3, 3, "mm")
  ) +
  guides(color = guide_legend(ncol = 2))
ggsave(
  "/vol/projects/yzhang/500FG_aging/output/14_AA/immune_aging_trajectory_Age_binned_segmented_auto.pdf",
  p1, width = 80, height = 70, units = "mm"
)
p1

