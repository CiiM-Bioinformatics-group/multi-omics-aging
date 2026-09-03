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

hormone <- read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_log2_hormone_levels.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)
head(hormone)
dim(hormone)
any(is.na(hormone))
sum(is.na(hormone))

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv")
rownames(basicPhenos) = basicPhenos$ID_500fg

## Filter a feature matrix and phenotype table to shared sample IDs and align rows in the same order
prep_filter_match <- function(feature_mat, basicPhenos) {
  idx <- intersect(rownames(feature_mat), rownames(basicPhenos))
  feature_filter <- feature_mat[idx, , drop = FALSE]
  pheno_match <- basicPhenos[rownames(feature_filter), , drop = FALSE]
  list(feature_filter = feature_filter, basicPhenos_filter_match = pheno_match)
}
# usage
res <- prep_filter_match(hormone, basicPhenos)
hormone_filter <- res$feature_filter
basicPhenos_filter_match <- res$basicPhenos_filter_match
##one missing gender value
basicPhenos_filter_match$Gender <- trimws(basicPhenos_filter_match$Gender)
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA

head(hormone_filter)
dim(hormone_filter)
head(basicPhenos_filter_match)
dim(basicPhenos_filter_match)


####model 2: age + hormone using FDR age not scale + gender
all_results <- list()

for (i in seq_len(ncol(hormone_filter))) {

  # Build a minimal modeling table
  data <- data.frame(
    feature = hormone_filter[, i],
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
    feature  = colnames(hormone_filter)[i],
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
write.csv(final_result, file = "/vol/projects/yzhang/500FG_aging/output/06_hormone/hormone_age_liner_gender_FDR.csv")

# 添加列：显著性分类（上调/下调/非显著）
final_result <- final_result %>% 
  mutate(
    sig = ifelse(p < 0.05 & estimate > 0, "upregulated",
                 ifelse(p < 0.05 & estimate < 0, "downregulated", "not significant")),
    log_p = -log10(pmax(p, .Machine$double.eps))  # 避免 p 为零导致计算错误
  )

# 分别筛选上调和下调中 log_p 最大的前 4 个
top3_upregulated <- final_result %>%
  filter(sig == "upregulated") %>%
  arrange(desc(log_p)) %>%
   head(3)


top3_downregulated <- final_result %>%
  filter(sig == "downregulated") %>%
  arrange(desc(log_p)) %>%
  head(3)

# 合并上调和下调的前 3 个显著点
top4_combined <- bind_rows(top3_upregulated, top3_downregulated)

# 设置保存火山图的分辨率和尺寸
ppi <- 300
png("/vol/projects/yzhang/500FG_aging/output/06_hormone/hormone_age_top3.png", 
    width = 5 * ppi, height = 4 * ppi, res = ppi)

# 绘制火山图
ggplot(final_result, aes(x = estimate, y = log_p)) + 
  geom_point(aes(color = sig), size = 3, alpha = 0.6) +  # 显著性颜色
  scale_color_manual(values = c(
    "upregulated" = "#FFCC99", 
    "downregulated" = "#B3DE69", 
    "not significant" = "grey"
  )) + 
  theme_classic() + 
  xlab("Estimate (Effect Size)") +
  ylab("-log10(p)") +
  ggtitle("Hormone ~ Age  Volcano Plot") +
  # 显示显著性统计信息
  annotate(
    "text", 
    x = max(final_result$estimate, na.rm = TRUE) -0.1, 
    y = max(final_result$log_p, na.rm = TRUE) - 4.5, 
    label = paste0(
     "Upregulated: ", sum(final_result$sig == "upregulated"), 
     " (", round(mean(final_result$sig == "upregulated") * 100, 2), "%)", 
     "\nDownregulated: ", sum(final_result$sig == "downregulated"), 
     " (", round(mean(final_result$sig == "downregulated") * 100, 2), "%)"
    ), 
    size = 4, 
    hjust = 0, 
    color = "black"
  ) +
  # 添加 top 4 显著点的名字，使用 geom_text_repel 避免重叠
  geom_text_repel(
    data = top4_combined,
    aes(x = estimate, y = log_p, label = hormone),
    size = 2, 
    color = "black",
    box.padding = 0.5,  # 控制文本与点的距离
    point.padding = 0.5,  # 控制点与文本的距离
    max.overlaps = Inf  # 避免重叠
  )

dev.off()

save.image("/vol/projects/yzhang/500FG_aging/output/06_hormone/hormone_gender.Rdata")

load("/vol/projects/yzhang/500FG_aging/output/06_hormone/hormone_gender.Rdata")

final_result_sig <- final_result %>%
  filter(sig == "sig") %>%
  select(feature, estimate, p, padj)
head(final_result_sig)
dim(final_result_sig)
write_excel_csv(final_result_sig,"/vol/projects/yzhang/500FG_aging/output/06_hormone/hormone_age_liner_gender_FDR_sig.csv")
