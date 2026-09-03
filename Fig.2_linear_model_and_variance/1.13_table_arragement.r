library(openxlsx)

cytokine <- read.csv("/vol/projects/yzhang/500FG_aging/output/01_cytokine/cytokine_age_liner_gender_FDR.csv")
protein <- read.csv( "/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_age_liner_gender_FDR.csv")
metabolite <- read.csv("/vol/projects/yzhang/500FG_aging/output/03_metabolite/metabolite_age_linear_gender_FDR.csv")
cellcount <- read.csv("/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_age_rawless_gender_FDR.csv")
immunoglobulin <- read.csv("/vol/projects/yzhang/500FG_aging/output/08_immunoglobulin/immunoglobulin_age_gender_FDR.csv")
hormone <- read.csv("/vol/projects/yzhang/500FG_aging/output/06_hormone/hormone_age_liner_gender_FDR.csv")
microbiome <- read.csv("/vol/projects/yzhang/500FG_aging/output/10_microbiome/microbiome_age_gender_FDR.csv")
platelet <- read.csv("/vol/projects/yzhang/500FG_aging/output/09_platelet/platelet_age_gender_FDR.csv")

## ============================================================
## 500FG aging: p < 0.05 结果整理为 supplementary Excel
## 依赖: openxlsx
## ============================================================

## ---------- 1. 读入数据 ----------
cytokine <- read.csv(
  "/vol/projects/yzhang/500FG_aging/output/01_cytokine/cytokine_age_liner_gender_FDR.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
protein <- read.csv(
  "/vol/projects/yzhang/500FG_aging/output/02_Olink/olink_age_liner_gender_FDR.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
metabolite <- read.csv(
  "/vol/projects/yzhang/500FG_aging/output/03_metabolite/metabolite_age_linear_gender_FDR.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
cellcount <- read.csv(
  "/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_age_rawless_gender_FDR.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
immunoglobulin <- read.csv(
  "/vol/projects/yzhang/500FG_aging/output/08_immunoglobulin/immunoglobulin_age_gender_FDR.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
hormone <- read.csv(
  "/vol/projects/yzhang/500FG_aging/output/06_hormone/hormone_age_liner_gender_FDR.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
microbiome <- read.csv(
  "/vol/projects/yzhang/500FG_aging/output/10_microbiome/microbiome_age_gender_FDR.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)
platelet <- read.csv(
  "/vol/projects/yzhang/500FG_aging/output/09_platelet/platelet_age_gender_FDR.csv",
  check.names = FALSE, stringsAsFactors = FALSE
)

methylation <- readRDS("/vol/projects/yzhang/500FG_aging/output/05_methylation/methylation_age_liner_gender_FDR.RDS")

## 兼容列名大小写 / 常见别名
find_col <- function(df, candidates) {
  nms <- names(df)
  nms_l <- tolower(gsub("[^a-z0-9]", "", nms))
  cand_l <- tolower(gsub("[^a-z0-9]", "", candidates))
  hit <- match(cand_l, nms_l)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) {
    stop(
      "找不到列：", paste(candidates, collapse = " / "),
      "\n实际列名：", paste(nms, collapse = ", ")
    )
  }
  nms[hit[1]]
}

format_sig_table <- function(df, p_cutoff = 0.05) {
  feat_col <- find_col(df, c("feature", "Feature", "features","Mvalue"))
  est_col  <- find_col(df, c("estimate", "Estimate", "coef", "beta"))
  p_col    <- find_col(df, c("p", "P", "pvalue", "p_value", "P value", "Pvalue"))
  fdr_col  <- find_col(df, c("padj", "FDR", "fdr", "qvalue", "q_value",
                             "FDR q value (BH)", "p.adjust"))

  out <- data.frame(
    Feature = as.character(df[[feat_col]]),
    Estimate = as.numeric(df[[est_col]]),
    `P value` = as.numeric(df[[p_col]]),
    `FDR q value (BH)` = as.numeric(df[[fdr_col]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  out <- out[!is.na(out[["P value"]]) & out[["P value"]] < p_cutoff, , drop = FALSE]
  out <- out[order(out[["P value"]], out[["FDR q value (BH)"]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

protein_sig      <- format_sig_table(protein)
metabolite_sig   <- format_sig_table(metabolite)
cytokine_sig     <- format_sig_table(cytokine)
cellcount_sig    <- format_sig_table(cellcount)
microbiome_sig   <- format_sig_table(microbiome)
hormone_sig      <- format_sig_table(hormone)
immune_sig       <- format_sig_table(rbind(immunoglobulin, platelet))
methylation_sig  <- if (!is.null(methylation)) format_sig_table(methylation) else NULL

## ---------- 3. 按图三顺序写 Excel（数字保持 numeric，避免变成文本） ----------
out_xlsx <- "/vol/projects/yzhang/500FG_aging/output/Supplementary_tables_age_associated.xlsx"

## Excel sheet 名最长 31 字符；下面名称均在限制内
sheet_specs <- list(
  list(
    sheet = "Sup table 1. Methylation_sig",
    title = "Supplementary table 1. CpG sites significantly associated with chronological age.",
    data  = methylation_sig
  ),
  list(
    sheet = "Sup table 2. Proteins_sig",
    title = "Supplementary table 2. Proteins significantly associated with chronological age.",
    data  = protein_sig
  ),
  list(
    sheet = "Sup table 3. Metabolites_sig",
    title = "Supplementary table 3. Metabolites significantly associated with chronological age.",
    data  = metabolite_sig
  ),
  list(
    sheet = "Sup table 4. Cytokine_sig",
    title = "Supplementary table 4. Cytokine responces significantly associated with chronological age.",
    data  = cytokine_sig
  ),
  list(
    sheet = "Sup table 5. Cellcount_sig",
    title = "Supplementary table 5. Immune cell counts significantly associated with chronological age.",
    data  = cellcount_sig
  ),
  list(
    sheet = "Sup table 6. Microbiomes_sig",
    title = "Supplementary table 6. Microbiome features significantly associated with chronological age.",
    data  = microbiome_sig
  ),
  list(
    sheet = "Sup table 7. endocrine_sig",
    title = "Supplementary table 7. Circulating endocrine traits significantly associated with chronological age.",
    data  = hormone_sig
  ),
  list(
    sheet = "Sup table 8. immune_sig",
    title = "Supplementary table 8. Circulating immune traits significantly associated with chronological age.",
    data  = immune_sig
  )
)

wb <- createWorkbook()

style_title <- createStyle(
  fontSize = 12, fontName = "Calibri", textDecoration = "bold",
  wrapText = TRUE, valign = "center", halign = "left"
)

## 表头：去掉底部横线
style_header <- createStyle(
  fontSize = 11, fontName = "Calibri", textDecoration = "bold",
  valign = "center", halign = "center"
)

style_text <- createStyle(
  fontSize = 11, fontName = "Calibri",
  halign = "left", valign = "center",
  numFmt = "@"   # Feature 强制文本，防止被 Excel 当成日期
)

## Estimate 列：3 位有效数字
style_est <- createStyle(
  fontSize = 11, fontName = "Calibri",
  halign = "right", valign = "center",
  numFmt = "0.00E+00"
)

## P value / FDR 列：2 位有效数字的科学计数法
style_p <- createStyle(
  fontSize = 11, fontName = "Calibri",
  halign = "right", valign = "center",
  numFmt = "0.0E+00"
)

write_sup_sheet <- function(wb, spec) {
  if (is.null(spec$data)) {
    message("跳过（无数据/无文件）：", spec$sheet)
    return(invisible(NULL))
  }

  addWorksheet(wb, spec$sheet, gridLines = TRUE)
  n <- nrow(spec$data)

  ## 第 1 行：表题（合并 A1:D1）
  writeData(wb, spec$sheet, spec$title, startCol = 1, startRow = 1, colNames = FALSE)
  mergeCells(wb, spec$sheet, cols = 1:4, rows = 1)
  addStyle(wb, spec$sheet, style_title, rows = 1, cols = 1, gridExpand = TRUE)
  setRowHeights(wb, spec$sheet, rows = 1, heights = 22)

  ## 第 2 行：列名；第 3 行起：数据
  writeData(
    wb, spec$sheet, spec$data,
    startCol = 1, startRow = 2,
    colNames = TRUE, rowNames = FALSE,
    keepNA = FALSE
  )

  addStyle(wb, spec$sheet, style_header, rows = 2, cols = 1:4, gridExpand = TRUE)

  if (n > 0) {
    data_rows <- 3:(n + 2)
    addStyle(wb, spec$sheet, style_text, rows = data_rows, cols = 1,   gridExpand = TRUE)
    addStyle(wb, spec$sheet, style_est,  rows = data_rows, cols = 2,   gridExpand = TRUE)
    addStyle(wb, spec$sheet, style_p,    rows = data_rows, cols = 3:4, gridExpand = TRUE)
  }

  ## 列宽：Feature 加宽，数值列适中
  setColWidths(wb, spec$sheet, cols = 1, widths = 36)
  setColWidths(wb, spec$sheet, cols = 2, widths = 16)
  setColWidths(wb, spec$sheet, cols = 3, widths = 16)
  setColWidths(wb, spec$sheet, cols = 4, widths = 20)

  freezePane(wb, spec$sheet, firstActiveRow = 3, firstActiveCol = 1)
  setRowHeights(wb, spec$sheet, rows = 2, heights = 18)
  if (n > 0) setRowHeights(wb, spec$sheet, rows = 3:(n + 2), heights = 16)

  invisible(NULL)
}

invisible(lapply(sheet_specs, function(spec) write_sup_sheet(wb, spec)))

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

message("已写出：", out_xlsx)
