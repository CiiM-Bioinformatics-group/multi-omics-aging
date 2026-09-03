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

# metabolite <- readRDS("/vol/projects/yzhang/500FG_aging/input/metabolite/metabolite_trans.rds")
# head(metabolite)
# dim(metabolite)

metabolite <- read_tsv("/vol/projects/CIIM/cohorts_old/500FG/metabolic/500FG raw metabolite_fromJianbo.tsv")
head(metabolite)
dim(metabolite)
metabolite_num <- metabolite[,c(22:479)]
metabolite_num <- as.data.frame(metabolite_num)
rownames(metabolite_num) <- metabolite$ metabolite_identification
head(metabolite_num)
dim(metabolite_num)
any(is.na(metabolite_num))
# 转置 cytokine_filter 数据框
metabolite_trans <- as.data.frame(t(metabolite_num))
colnames(metabolite_trans) <- rownames(metabolite_num)
rownames(metabolite_trans) <- colnames(metabolite_num)
dim(metabolite_trans)
head(metabolite_trans)

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv",row.names=1)
rownames(basicPhenos) = basicPhenos$ID_500fg

## Filter a feature matrix and phenotype table to shared sample IDs and align rows in the same order
prep_filter_match <- function(feature_mat, basicPhenos) {
  idx <- intersect(rownames(feature_mat), rownames(basicPhenos))
  feature_filter <- feature_mat[idx, , drop = FALSE]
  pheno_match <- basicPhenos[rownames(feature_filter), , drop = FALSE]
  list(feature_filter = feature_filter, basicPhenos_filter_match = pheno_match)
}
# usage
res <- prep_filter_match(metabolite_trans, basicPhenos)
metabolite_filter <- res$feature_filter
basicPhenos_filter_match <- res$basicPhenos_filter_match
##one missing gender value
basicPhenos_filter_match$Gender <- trimws(basicPhenos_filter_match$Gender)
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA


head(metabolite_filter)
dim(metabolite_filter)
head(basicPhenos_filter_match)
dim(basicPhenos_filter_match)

####model 2: age + protein using FDR age not scale + gender
all_results <- list()

for (i in seq_len(ncol(metabolite_filter))) {

  # Build a minimal modeling table
  data <- data.frame(
    feature = log(metabolite_filter[, i] + 1),
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
    feature  = colnames(metabolite_filter)[i],
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
write.csv(final_result, file = "/vol/projects/yzhang/500FG_aging/output/03_metabolite/metabolite_age_linear_gender_FDR.csv")

save.image("/vol/projects/yzhang/500FG_aging/output/03_metabolite/metabolite_gender.Rdata")

load("/vol/projects/yzhang/500FG_aging/output/03_metabolite/metabolite_gender.Rdata")

final_result_sig <- final_result %>%
  filter(sig == "sig") %>%
  select(feature, estimate, p, padj)
head(final_result_sig)
dim(final_result_sig)
write_excel_csv(final_result_sig,"/vol/projects/yzhang/500FG_aging/output/03_metabolite/metabolite_age_linear_gender_FDR_sig.csv")
