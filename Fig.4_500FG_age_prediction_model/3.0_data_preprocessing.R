# 3.0 Input preparation for 3.1-3.13
#     Inverse rank normalization of cytokine, proteomics, metabolomics and microbiome data;
#     cell counts with descriptive column names.

# ---- Paths ----
data_dir <- "data"

# ---- Raw inputs (samples in rows) ----
cytokine   <- readRDS(file.path(data_dir, "cytokine_filter_transposed.rds"))
olink      <- readRDS(file.path(data_dir, "NPX_FG500.RDS"))
metabolite <- readRDS(file.path(data_dir, "metabolite_trans.rds"))
microbiome <- read.table(file.path(data_dir, "500FG_microbiome_pathways.txt"),
                         header = TRUE, sep = "", stringsAsFactors = FALSE)
cellcounts <- read.table(file.path(data_dir, "500FG_inverse_rank_normalized_cellcounts.txt"),
                         header = TRUE, stringsAsFactors = FALSE)
cellcounts_info <- read.table(file.path(data_dir, "500FG_cellcounts_info.txt"))

# ---- Inverse rank normalization (ties broken at random) ----
rank_normalize <- function(value) {
  qnorm((rank(value, na.last = "keep", ties.method = "random") - 0.5) / sum(!is.na(value)))
}
normalize_table <- function(x) {
  set.seed(123)
  out <- as.data.frame(apply(x, 2, rank_normalize))
  rownames(out) <- rownames(x)
  out
}

write.csv(normalize_table(metabolite), file.path(data_dir, "metabolite_rank_normalized.csv"), row.names = TRUE)
write.csv(normalize_table(cytokine),   file.path(data_dir, "cytokine_rank_normalized.csv"), row.names = TRUE)
write.csv(normalize_table(olink),      file.path(data_dir, "olink_rank_normalized.csv"), row.names = TRUE)

microbiome_rn <- normalize_table(microbiome)
microbiome_rn <- microbiome_rn[order(rownames(microbiome_rn)), ]
write.csv(microbiome_rn, file.path(data_dir, "microbiome_rank_normalized.csv"), row.names = TRUE)

# ---- Cell counts: replace codes (IT1, IT2, ...) with cell population names ----
colnames(cellcounts) <- cellcounts_info$finalName[match(colnames(cellcounts), cellcounts_info$code)]
write.csv(cellcounts, file.path(data_dir, "cellcounts_name_replace.csv"), row.names = TRUE)
