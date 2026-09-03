library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(networkD3)
library(tibble)
library(ggalluvial)
library(RColorBrewer)
library(scales)
library(viridis)

platelet <- read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_log2_platelet_count.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)
head(platelet)
dim(platelet)
any(is.na(platelet))
sum(is.na(platelet))

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv")
rownames(basicPhenos) = basicPhenos$ID_500fg

## Filter a feature matrix and phenotype table to shared sample IDs and align rows in the same order
prep_filter_match <- function(feature_mat, basicPhenos) {
  idx <- rownames(feature_mat)[rownames(feature_mat) %in% rownames(basicPhenos)]  # keep feature_mat order
  feature_filter <- feature_mat[idx, , drop = FALSE]
  pheno_match <- basicPhenos[rownames(feature_filter), , drop = FALSE]
  list(feature_filter = feature_filter, basicPhenos_filter_match = pheno_match)
}
res <- prep_filter_match(platelet, basicPhenos)
platelet_filter <- res$feature_filter
basicPhenos_filter_match <- res$basicPhenos_filter_match
##one missing gender value
basicPhenos_filter_match$Gender <- trimws(basicPhenos_filter_match$Gender)
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA

head(platelet_filter)
dim(platelet_filter)
head(basicPhenos_filter_match)
dim(basicPhenos_filter_match)

####model 2: age + platelet using FDR age not scale + gender
all_results <- list()

for (i in seq_len(ncol(platelet_filter))) {

  # Build a minimal modeling table
  data <- data.frame(
    feature = platelet_filter[, i],
    age = basicPhenos_filter_match$Age,
    Gender = factor(basicPhenos_filter_match$Gender)  # convert chr to factor
  )

  # Linear regression
  mod <- lm(feature ~ age + Gender, data = data)

  # Extract coefficient and p-value for age (not Gender)
  cf <- summary(mod)$coefficients
  res <- cf["age", c("Estimate", "Pr(>|t|)")]

  # Store results
  result <- data.frame(
    feature  = colnames(platelet_filter)[i],
    estimate = res["Estimate"],
    p        = res["Pr(>|t|)"]
  )

  # Signed -log10(p)
  result$value <- -log10(result$p) * result$estimate

  all_results[[i]] <- result
}

final_result <- do.call(rbind, all_results)
final_result$padj <- p.adjust(final_result$p, method = "BH")
final_result$sig <- ifelse(final_result$padj < 0.05, "sig", "no")

head(final_result)
dim(final_result)
table(final_result$sig)
write.csv(final_result, file = "/vol/projects/yzhang/500FG_aging/output/09_platelet/platelet_age_gender_FDR.csv")
