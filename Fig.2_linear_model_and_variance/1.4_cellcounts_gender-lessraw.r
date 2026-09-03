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
library(readr)

cellCounts = read.csv("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_cellCounts_cellPerc.csv",row.names=1)
head(cellCounts)
cellCounts_select = read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_inverse_rank_normalized_cellcounts.txt", header = TRUE, stringsAsFactors = FALSE) ##less
head(cellCounts_select)
dim(cellCounts_select)
cellCounts <- cellCounts[,colnames(cellCounts) %in% colnames(cellCounts_select)]
head(cellCounts)
dim(cellCounts)
any(is.na(cellCounts))

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv")
rownames(basicPhenos) = basicPhenos$ID_500fg

## Filter a feature matrix and phenotype table to shared sample IDs and align rows in the same order
prep_filter_match <- function(feature_mat, basicPhenos) {
  idx <- rownames(feature_mat)[rownames(feature_mat) %in% rownames(basicPhenos)]  # keep feature_mat order
  feature_filter <- feature_mat[idx, , drop = FALSE]
  pheno_match <- basicPhenos[rownames(feature_filter), , drop = FALSE]
  list(feature_filter = feature_filter, basicPhenos_filter_match = pheno_match)
}
res <- prep_filter_match(cellCounts, basicPhenos)
cellCounts_filter <- res$feature_filter
basicPhenos_filter_match <- res$basicPhenos_filter_match
##one missing gender value
basicPhenos_filter_match$Gender <- trimws(basicPhenos_filter_match$Gender)
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA

head(cellCounts_filter)
dim(cellCounts_filter)
head(basicPhenos_filter_match)
dim(basicPhenos_filter_match)

####model 2: age + hormone using FDR age not scale + gender
all_results <- list()

for (i in seq_len(ncol(cellCounts_filter))) {

  # Build a minimal modeling table
  data <- data.frame(
    feature =log(cellCounts_filter[, i]+ 1),
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
    feature  = colnames(cellCounts_filter)[i],
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

##replace the name
cellcounts_info <- read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_cellcounts_info.txt")
head(cellcounts_info)
match_indices <- match(final_result$feature, cellcounts_info$code)
final_result$feature <- cellcounts_info$finalName[match_indices]
head(final_result)
dim(final_result)

write.csv(final_result, file = "/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_age_rawless_gender_FDR.csv")
print(table(final_result$sig))

save.image("/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_age_lessraw_FDR_namereplace_gender.Rdata")

load("/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_age_lessraw_FDR_namereplace_gender.Rdata")

final_result_sig <- final_result %>%
  filter(sig == "sig") %>%
  select(feature, estimate, p, padj)
head(final_result_sig)
dim(final_result_sig)
write_excel_csv(final_result_sig, "/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_age_rawless_gender_FDR_sig.csv")
