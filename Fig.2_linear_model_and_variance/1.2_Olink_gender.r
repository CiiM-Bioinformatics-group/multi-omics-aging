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

load("/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_gender.RData")

olink <- readRDS("/vol/projects/CIIM/cohorts_old/500FG/olinkData/NPX_FG500.RDS")
head(olink)
dim(olink)
any(is.na(olink))

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv",row.names=1)
rownames(basicPhenos) = basicPhenos$ID_500fg
# Get shared sample IDs between olink and phenotype table
idx <- intersect(rownames(olink), basicPhenos$ID_500fg)
length(idx)  # 458
# Subset olink to shared IDs (keeps olink row order)
olink_filter <- olink[rownames(olink) %in% idx, , drop = FALSE]
head(olink_filter)
dim(olink_filter)
# Subset phenotypes and reorder to exactly match olink_filter rownames
basicPhenos_filter_match <- basicPhenos %>%
  filter(ID_500fg %in% rownames(olink_filter)) %>%
  slice(match(rownames(olink_filter), ID_500fg))
# Sanity check: IDs should align 1-to-1
sum(rownames(olink_filter) == basicPhenos_filter_match$ID_500fg)
head(basicPhenos_filter_match)
dim(basicPhenos_filter_match)  # 458 x (pheno columns)

##one missing gender value
basicPhenos_filter_match$Gender <- trimws(basicPhenos_filter_match$Gender)
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA

table(basicPhenos$Gender)
idx <- which(trimws(as.character(basicPhenos$Gender)) == "")
basicPhenos[idx, ] 

# model 1: age + protein (age not scaled) + Gender
all_results <- list()

for (i in seq_len(ncol(olink_filter))) {

  # Build a minimal modeling table
  data <- data.frame(
    feature = log(olink_filter[, i] + 1),
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
    feature  = colnames(olink_filter)[i],
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

write.csv(final_result,
          file = "/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_age_liner_gender_FDR.csv")

# 添加列：显著性分类（上调/下调/非显著）
final_result <- final_result %>% 
  mutate(
    sig = ifelse(p < 0.05 & estimate > 0, "upregulated",
                 ifelse(p < 0.05 & estimate < 0, "downregulated", "not significant")),
    log_p = -log10(p) 
  )

# 分别筛选上调和下调中 log_p 最大的前 4 个
top3_upregulated <- final_result %>%
  filter(sig == "upregulated") %>%
  arrange(desc(log_p)) %>%
  slice(1:3)

top3_downregulated <- final_result %>%
  filter(sig == "downregulated") %>%
  arrange(desc(log_p)) %>%
  slice(1:3)

# 合并上调和下调的前 3 个显著点
top4_combined <- bind_rows(top3_upregulated, top3_downregulated)

# 设置保存火山图的分辨率和尺寸
ppi <- 300
png("/vol/projects/yzhang/500FG_aging/output/02_Olink/volcano_age_olink_liner_top3(2).png", 
    width = 5 * ppi, height = 4 * ppi, res = ppi)

# 绘制火山图
ggplot(final_result, aes(x = estimate, y = log_p)) + 
  geom_point(aes(color = sig), size = 3, alpha = 0.6) +  # 显著性颜色
  scale_color_manual(values = c(
    "upregulated" = "#FFB3E6", 
    "downregulated" = "#80B1D3", 
    "not significant" = "grey"
  )) + 
  theme_classic() + 
  xlab("Estimate (Effect Size)") +
  ylab("-log10(p)") +
  ggtitle("Proteomics ~ Age Volcano Plot") +
  # 显示显著性统计信息
  annotate(
    "text", 
    x = max(final_result$estimate, na.rm = TRUE) * -1.3, 
    y = max(final_result$log_p, na.rm = TRUE) - 3, 
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
    aes(x = estimate, y = log_p, label = olink),
    size = 3.5, 
    color = "black",
    box.padding = 0.5,  # 控制文本与点的距离
    point.padding = 0.5,  # 控制点与文本的距离
    #max.overlaps = Inf  # 避免重叠
  )

dev.off()

table(final_result$sig)
protein_sig <- final_result[final_result$sig=="sig",]
dim(protein_sig)

###compare with previous version
protein_age <- read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_age_linermodel_ynorm(FDR).csv",row.names=1)
head(protein_age)
dim(protein_age)
protein_age_sig <- protein_age[protein_age$sig=="sig",]
head(protein_age_sig)
dim(protein_age_sig)

length(intersect(protein_sig$feature, protein_age_sig$olink))
setdiff(protein_sig$feature, protein_age_sig$olink) # in 1st but not in 2nd
setdiff(protein_age_sig$olink, protein_sig$feature)  # in 2nd but not in 1st

save.image("/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_gender.RData")

load("/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_gender.RData")

final_result_sig <- final_result %>%
  filter(sig == "sig") %>%
  select(feature, estimate, p, padj)
head(final_result_sig)
dim(final_result_sig)
write_excel_csv(final_result_sig, "/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_age_liner_gender_FDR_sig.csv")
