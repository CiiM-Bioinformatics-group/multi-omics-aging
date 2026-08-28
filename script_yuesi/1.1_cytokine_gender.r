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

cytokine <- read.csv("/vol/projects/yzhang/500FG_aging/input/cytokine/filtered_cytokine_filled_trans.csv")
rownames(cytokine) <- cytokine[,1]
cytokine <- cytokine[,-1]
head(cytokine)
dim(cytokine) #489samples 91 cytokines
sum(is.na(cytokine))

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
res <- prep_filter_match(cytokine, basicPhenos)
cytokine_filter <- res$feature_filter
basicPhenos_filter_match <- res$basicPhenos_filter_match
##one missing gender value
basicPhenos_filter_match$Gender <- trimws(basicPhenos_filter_match$Gender)
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA

head(cytokine_filter)
dim(cytokine_filter)
head(basicPhenos_filter_match)
dim(basicPhenos_filter_match)

# model 1: age + cytokine (age not scaled) + Gender
all_results <- list()

for (i in seq_len(ncol(cytokine_filter))) {

  # Build a minimal modeling table
  data <- data.frame(
    feature = log(cytokine_filter[, i] + 1),
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
    feature  = colnames(cytokine_filter)[i],
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

write.csv(final_result, file = "/vol/projects/yzhang/500FG_aging/output/01_cytokine/cytokine_age_liner_gender_FDR.csv")

###compare with previous version
protein_age <- read.csv("/vol/projects/yzhang/500FG_aging/output/01_cytokine/cytokine_age_liner_FDR.csv",row.names=1)
head(protein_age)
dim(protein_age)
protein_age_sig <- protein_age[protein_age$sig=="sig",]
head(protein_age_sig)
dim(protein_age_sig)

save.image("/vol/projects/yzhang/500FG_aging/output/01_cytokine/cytokine_gender.RData")
