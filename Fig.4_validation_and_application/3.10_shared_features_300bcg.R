# 3.10 Features measured in both 500FG and 300BCG (used in 3.11 and 3.12)
#      Non-methylation: proteomics, cytokines, metabolites, cell counts with matching names;
#      methylation: CpGs with matching probe IDs (first 10 characters).
#      Run from the repository root: Rscript validation/<this script>

library(dplyr)
library(readxl)

# ---- Paths ----
data_dir <- "data"
bcg_dir  <- file.path(data_dir, "300bcg")

# =============================================================================
# 1. Non-methylation features
# =============================================================================
# ---- 300BCG (samples in rows) ----
cytokine_300bcg    <- read.csv(file.path(bcg_dir, "filtered_cytokine.csv"), row.names = 1)
olink_300bcg       <- readRDS(file.path(bcg_dir, "olink_baseline_filter.RDS"))
metabolite_300bcg  <- read.csv(file.path(bcg_dir, "metabolite_BCG300_fliter_data.csv"), row.names = 1)
cellCounts_300bcg  <- read.csv(file.path(bcg_dir, "300BCG_IRT_all_replaced_name.csv"), row.names = 1, check.names = FALSE)
basicPhenos_300bcg <- read_excel(file.path(bcg_dir, "201103 300BCG Age and sex data.xlsx"))

# Cytokines: baseline (corr.1) only, names harmonized with 500FG (e.g. IL6_S.aureus_PBMC_24h)
colnames(cytokine_300bcg) <- colnames(cytokine_300bcg) |>
  gsub("T3", "S.aureus_24h", x = _) |>
  gsub("T4", "M.tuberculosis_24h", x = _) |>
  gsub("W3", "S.aureus_7days", x = _) |>
  gsub("W4", "M.tuberculosis_7days", x = _)
cytokine_300bcg_base <- cytokine_300bcg[, grep("corr\\.1", colnames(cytokine_300bcg), value = TRUE)]
colnames(cytokine_300bcg_base) <- sapply(gsub("_corr\\.1$", "", colnames(cytokine_300bcg_base)), function(name) {
  parts <- unlist(strsplit(name, "_"))
  if (length(parts) >= 2) {
    paste0(paste(parts[1:2], collapse = "_"), "_PBMC",
           if (length(parts) > 2) paste0("_", paste(parts[3:length(parts)], collapse = "_")) else "")
  } else paste0(name, "_PBMC")
})

# Cell counts: sample 118 has no measurements
cellCounts_300bcg <- cellCounts_300bcg[rownames(cellCounts_300bcg) != "118", ]

common_300bcg <- Reduce(intersect, list(rownames(cytokine_300bcg_base), rownames(olink_300bcg),
                                        rownames(metabolite_300bcg), rownames(cellCounts_300bcg),
                                        basicPhenos_300bcg$PatientID))

# ---- 500FG (samples in rows) ----
cytokine_500fg    <- readRDS(file.path(data_dir, "cytokine_filter_transposed.rds"))
olink_500fg       <- readRDS(file.path(data_dir, "NPX_FG500.RDS"))
metabolite_500fg  <- readRDS(file.path(data_dir, "metabolite_trans.rds"))
cellcounts_500fg  <- read.csv(file.path(data_dir, "cellcounts_name_replace.csv"), row.names = 1, check.names = FALSE)
microbiome_500fg  <- read.table(file.path(data_dir, "500FG_microbiome_pathways.txt"),
                                header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
basicPhenos_500fg <- read.csv(file.path(data_dir, "Age_group_basicPhenos.csv"))

common_500fg <- sort(Reduce(intersect, list(rownames(cytokine_500fg), rownames(olink_500fg),
                                            rownames(metabolite_500fg), rownames(cellcounts_500fg),
                                            rownames(microbiome_500fg), basicPhenos_500fg$ID_500fg)))

# ---- Shared feature names per layer ----
shared <- list(
  Olink      = intersect(colnames(olink_500fg), colnames(olink_300bcg)),
  Cytokine   = intersect(colnames(cytokine_500fg), colnames(cytokine_300bcg_base)),
  Metabolite = intersect(colnames(metabolite_500fg), colnames(metabolite_300bcg)),
  CellCounts = intersect(colnames(cellcounts_500fg), colnames(cellCounts_300bcg))
)
print(lengths(shared))

# Combined table: layer prefix on feature names, layers in the order above
build_table <- function(layers, samples) {
  do.call(cbind, lapply(names(shared), function(l) {
    x <- layers[[l]][samples, shared[[l]], drop = FALSE]
    colnames(x) <- paste0(l, "_", colnames(x))
    x
  }))
}

fg_tab_common_feat <- build_table(list(Olink = olink_500fg, Cytokine = cytokine_500fg,
                                       Metabolite = metabolite_500fg, CellCounts = cellcounts_500fg), common_500fg)

bcg_tab_common_feat <- build_table(list(Olink = olink_300bcg, Cytokine = cytokine_300bcg_base,
                                        Metabolite = metabolite_300bcg, CellCounts = cellCounts_300bcg), common_300bcg)
rownames(bcg_tab_common_feat) <- sprintf("%03d", as.integer(rownames(bcg_tab_common_feat)))   # 4 -> "004"
bcg_tab_common_feat <- bcg_tab_common_feat[order(rownames(bcg_tab_common_feat)), , drop = FALSE]

# Phenotypes in the same sample order
basicPhenos_common_500fg_ordered <- basicPhenos_500fg[match(common_500fg, basicPhenos_500fg$ID_500fg), , drop = FALSE]
basicPhenos_common_300bcg <- as.data.frame(basicPhenos_300bcg[match(common_300bcg, basicPhenos_300bcg$PatientID), , drop = FALSE])
rownames(basicPhenos_common_300bcg) <- sprintf("%03d", as.integer(basicPhenos_common_300bcg$PatientID))
basicPhenos_common_300bcg <- basicPhenos_common_300bcg[order(rownames(basicPhenos_common_300bcg)), , drop = FALSE]

write.csv(fg_tab_common_feat, file.path(bcg_dir, "500fg_tab_common_feat_replaced.csv"), row.names = TRUE)
write.csv(bcg_tab_common_feat, file.path(bcg_dir, "300bcg_tab_common_feat_replaced.csv"), row.names = TRUE)
write.csv(basicPhenos_common_500fg_ordered, file.path(bcg_dir, "basicPhenos_tab_common_500fg_ordered.csv"), row.names = TRUE)
write.csv(basicPhenos_common_300bcg, file.path(bcg_dir, "basicPhenos_common_300bcg_feat_replaced.csv"), row.names = TRUE)

# =============================================================================
# 2. Shared CpGs
# =============================================================================
Mvalue_300bcg <- readRDS(file.path(bcg_dir, "mvalue_v1_trans.rds"))
Mvalue_500fg  <- readRDS(file.path(data_dir, "Mvalue_trans.rds"))

common_cols <- intersect(substr(colnames(Mvalue_500fg), 1, 10), colnames(Mvalue_300bcg))
length(common_cols)

write.csv(common_cols, file.path(bcg_dir, "common_cols.csv"))
write.csv(Mvalue_300bcg[, common_cols, drop = FALSE], file.path(bcg_dir, "Mvalue_300BCG_common.csv"))
