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

Mvalue <- readRDS("/vol/projects/CIIM/methylation/500FG/1QC/Output/500FG.Mvalue.rds")
Mvalue_trans <- as.data.frame(t(Mvalue))
# colnames(Mvalue_trans) <- rownames(Mvalue)
# rownames(Mvalue_trans) <- colnames(Mvalue)

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv")
rownames(basicPhenos) = basicPhenos$ID_500fg
# prepare cytokine_filter and EAA_filter_match
idx = intersect(rownames(Mvalue_trans), basicPhenos$ID_500fg) 
Mvalue_filter = Mvalue_trans[which(rownames(Mvalue_trans) %in% idx),]
# head(Mvalue_filter)  #458samples,1596 metabolite
dim(Mvalue_filter)

basicPhenos_filter = basicPhenos %>% filter(ID_500fg %in% idx)
head(basicPhenos_filter) #458samples
dim(basicPhenos_filter)

sum(rownames(Mvalue_filter) == basicPhenos_filter$ID_500fg)

basicPhenos_filter_match = basicPhenos_filter[match(rownames(Mvalue_filter), basicPhenos_filter$ID_500fg), ]

sum(basicPhenos_filter_match$ID_500fg == rownames(Mvalue_filter))
head(basicPhenos_filter_match) 
dim(basicPhenos_filter_match) #458samples
##one missing gender value
basicPhenos_filter_match$Gender <- trimws(basicPhenos_filter_match$Gender)
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA

#check distribution
# res = NULL
# for(i in 1:854368){
#   nor_check = shapiro.test(Mvalue_filter[,i])$p.value
#   result = data.frame(Mvalue_filter = colnames(Mvalue_filter)[i], pval = nor_check) 
#   res = rbind(res, result)
# }
# table(res$pval < 0.05) #non-normal 1573 normal 23

#model 1: age + protein 
all_results <- list()

for (i in 1:854368) {
    
   # Build a minimal modeling table
   data <- data.frame(methylation = Mvalue_filter[, i], age = basicPhenos_filter_match$Age, Gender = factor(basicPhenos_filter_match$Gender))
       
    
   # Linear regression
    mod <- lm(formula = methylation ~ age + Gender, data = data)
    
   # Extract coefficient and p-value for age (not Gender)
    cf <- summary(mod)$coefficients
    res <- cf[2, c("Estimate", "Pr(>|t|)")]
    
     # Store results
    result <- data.frame(
        Mvalue = colnames(Mvalue_filter)[i],
        estimate = res["Estimate"],
        p = res["Pr(>|t|)"]
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
                      
saveRDS(final_result, file = "/vol/projects/yzhang/500FG_aging/output/05_methylation/methylation_age_liner_gender_FDR.RDS")
