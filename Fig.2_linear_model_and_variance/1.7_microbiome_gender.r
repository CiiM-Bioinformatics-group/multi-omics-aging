library(readxl)
library(readr)
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
library(car)

load("/vol/projects/yzhang/500FG_aging/output/10_microbiome/microbiome_gender.RData")

microbiome <- read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_microbiome_pathways.txt", header = TRUE, sep = "", stringsAsFactors = FALSE,check.names = FALSE)
head(microbiome)
dim(microbiome)
any(is.na(microbiome))
sum(is.na(microbiome))

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv")
rownames(basicPhenos) = basicPhenos$ID_500fg

## Filter a feature matrix and phenotype table to shared sample IDs and align rows in the same order
prep_filter_match <- function(feature_mat, basicPhenos) {
  idx <- rownames(feature_mat)[rownames(feature_mat) %in% rownames(basicPhenos)]  # keep feature_mat order
  feature_filter <- feature_mat[idx, , drop = FALSE]
  pheno_match <- basicPhenos[rownames(feature_filter), , drop = FALSE]
  list(feature_filter = feature_filter, basicPhenos_filter_match = pheno_match)
}
res <- prep_filter_match(microbiome, basicPhenos)
microbiome_filter <- res$feature_filter
basicPhenos_filter_match <- res$basicPhenos_filter_match
##one missing gender value
basicPhenos_filter_match$Gender <- trimws(basicPhenos_filter_match$Gender)
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA

head(microbiome_filter)
dim(microbiome_filter)
head(basicPhenos_filter_match)
dim(basicPhenos_filter_match)

bad_cols <- which(apply(microbiome_filter, 2, function(x) any(x <= -1, na.rm = TRUE)))
length(bad_cols)
colnames(microbiome_filter)[bad_cols][1:10]
j <- bad_cols[1]
x <- microbiome_filter[, j]
summary(x)
range(x, na.rm = TRUE)
which(x < -1)              # 这些行会导致 NaN
x[x < -1][1:20]            # 看前 20 个触发值
rownames(microbiome_filter)[which(x < -1)][1:20]  # 对应样本ID   
sum(microbiome_filter < 0, na.rm=TRUE)                       

####model 2: age + immunoglobulinusing FDR age not scale + gender
all_results <- list()

for (i in seq_len(ncol(microbiome_filter))) {

    ###rank inverse normal transform since orginial data has minus value    
   #  x <- microbiome_filter[, i]
   # feature_irnt <- qnorm((rank(x, na.last = "keep", ties.method = "average") - 0.5) /
   #                    sum(!is.na(x)))

  # Build a minimal modeling table
  data <- data.frame(
    feature = asinh(microbiome_filter[, i]),  #microbiome_filter[, i]；feature_irnt
    age = basicPhenos_filter_match$Age,
    Gender = factor(basicPhenos_filter_match$Gender)  # convert chr to factor
  )
  
  # Remove incomplete rows (handles missing Gender or feature NA)
  data <- data[complete.cases(data), ]
  # Skip if too few samples or no variation
  if (nrow(data) < 10 || sd(data$feature) == 0) {
    all_results[[i]] <- data.frame(
      feature  = colnames(microbiome_filter)[i],
      estimate = NA_real_,
      p        = NA_real_
    )
    next
  }
  # Linear regression
  mod <- lm(feature ~ age + Gender, data = data)

  # Extract coefficient and p-value for age (not Gender)
  cf <- summary(mod)$coefficients
  res <- cf["age", c("Estimate", "Pr(>|t|)")]

  # Store results
  result <- data.frame(
    feature  = colnames(microbiome_filter)[i],
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
print(table(final_result$sig))
write.csv(final_result, file = "/vol/projects/yzhang/500FG_aging/output/10_microbiome/microbiome_age_gender_FDR.csv")

final_result <- final_result %>% 
  mutate(
    sig = case_when(
      !is.na(p) & p < 0.05 & estimate > 0 ~ "upregulated",
      !is.na(p) & p < 0.05 & estimate < 0 ~ "downregulated",
      TRUE ~ "not significant"
    ),
    log_p = -log10(pmax(p, .Machine$double.eps))
  )
top3_upregulated <- final_result %>%
  filter(sig == "upregulated") %>%
  arrange(desc(log_p)) %>%
  head(3)
top3_downregulated <- final_result %>%
  filter(sig == "downregulated") %>%
  arrange(desc(log_p)) %>%
  head(3)
top4_combined <- bind_rows(top3_upregulated, top3_downregulated)
ggplot(final_result, aes(x = estimate, y = log_p)) + 
  geom_point(aes(color = sig), size = 3, alpha = 0.6) +
  scale_color_manual(values = c(
    "upregulated" = "#FFCC99", 
    "downregulated" = "#CCEBC5", 
    "not significant" = "grey"
  )) + 
  theme_classic() + 
  xlab("Estimate (Effect Size)") +
  ylab("-log10(p)") +
  ggtitle("Microbiome ~ Age Volcano Plot") +
  annotate(
    "text",
    x = max(final_result$estimate, na.rm = TRUE),
    y = max(final_result$log_p, na.rm = TRUE),
    label = paste0(
      "Upregulated: ", sum(final_result$sig == "upregulated"), "\n",
      "Downregulated: ", sum(final_result$sig == "downregulated")
    ),
    hjust = 1, vjust = 1, size = 4
  ) +
  geom_text_repel(
    data = top4_combined,
    aes(label = feature),   # <- 这里改成 feature
    size = 2.8,
    color = "black",
    box.padding = 0.5,
    point.padding = 0.5,
    max.overlaps = Inf
  )

class(final_result$p)
length(unique(final_result$p))
head(final_result$p, 20)
length(unique(final_result$padj))

save.image("/vol/projects/yzhang/500FG_aging/output/10_microbiome/microbiome_gender.RData")

load("/vol/projects/yzhang/500FG_aging/output/10_microbiome/microbiome_gender.RData")

final_result_sig <- final_result %>%
  filter(sig == "sig") %>%
  select(feature, estimate, p, padj)
head(final_result_sig)
dim(final_result_sig)
write_excel_csv(final_result_sig,"/vol/projects/yzhang/500FG_aging/output/10_microbiome/microbiome_age_gender_FDR_sig.csv")
