# Age prediction from multi-omics data (500FG, external validation in 300BCG)

Scripts for predicting chronological age from single omics layers and from multi-omics models in the 500FG cohort, with external validation of TabPFN and other machine-learning models in 300BCG.

## Pipeline

| Script | Description | Main output |
|---|---|---|
| `prediction_utils.R` | Shared R functions | – |
| `3.0_data_preprocessing.R` | Inverse rank normalization; cell count names | `data/*_rank_normalized.csv`, `data/cellcounts_name_replace.csv` |
| **500FG: elastic net** | | |
| `3.1_single_omics_prediction.R` | Elastic net per layer (cytokine, proteomics, metabolomics, cell counts, microbiome), all samples of each layer | `results/prediction/*_results.rds` |
| `3.2_methylation_spearman_preparation.R` | Spearman correlation of each CpG with age in each training set (all methylation samples) | `results/methylation_spearman/` |
| `3.3_methylation_spearman_prediction.R` | Elastic net on CpGs with Spearman p < 1e-4 | `results/prediction/methylation_prediction_*.rds` |
| `3.4_multiomics_methylation_fdr_preparation.R` | Common samples of all six layers, 100 fixed splits, CpG Spearman FDR; CpGs with FDR < 0.05 | `results/fdr_common/` |
| `3.5_multiomics_prediction.R` | Multi-omics elastic net (all non-methylation features + FDR-selected CpGs) | `results/prediction/multi_prediction_*.rds` |
| `3.6_plot_r2_comparison.R` | R² across layers and the multi-omics model | `results/figures/3.6_multi-omics_r2_comparison.pdf` |
| **500FG: TabPFN and other models** | | |
| `3.7_omics_fdr_preparation.R` | Spearman FDR of non-methylation features (same samples and splits as 3.4); CSV inputs for Python | `results/fdr_common/*_FDR_common_results.json` |
| `3.8_top83_feature_selection.R` | Top 83 features per layer by FDR | `results/fdr_common/FDR_top83/` |
| `3.9_tabpfn_500fg_top83.py` | TabPFN, ElasticNet, RandomForest, XGBoost, CatBoost on top-83 features | `results/tabpfn_500fg/model_results_spear_top83.csv` |
| **500FG → 300BCG external validation** | | |
| `3.10_shared_cpgs_300bcg.R` | CpGs measured in both cohorts | `data/300bcg/common_cols.csv`, `Mvalue_300BCG_common.csv` |
| `3.11_top235_shared_cpg_selection.R` | Top 235 shared CpGs by FDR per split | `results/fdr_shared_top235/` |
| `3.12_tabpfn_500fg_300bcg_validation.py` | Train in 500FG, validate in 300BCG; saves fitted models | `results/tabpfn_300bcg/` |
| `3.13_plot_validation_two_panels.R` | Model comparison (500FG internal / 300BCG external) with paired t-tests | `results/figures/validation_two_panels_nature.pdf` |

Each model is evaluated on 100 train/test splits (70/30). In 3.4–3.13 all models share the same 100 splits of the samples with all six omics layers.

## Usage

Run the scripts in order from the repository root, each in a separate session:

```bash
Rscript 3.0_data_preprocessing.R
Rscript 3.1_single_omics_prediction.R
# ... 3.2 to 3.8
python 3.9_tabpfn_500fg_top83.py
Rscript 3.10_shared_cpgs_300bcg.R
Rscript 3.11_top235_shared_cpg_selection.R
python 3.12_tabpfn_500fg_300bcg_validation.py
Rscript 3.13_plot_validation_two_panels.R
```

3.2, 3.4 and 3.7 use `parallel::mclapply` (`n_cores` at the top of each script) and need substantial memory for the genome-wide methylation data.

## Input data (`data/`)

| File | Description |
|---|---|
| `Age_group_basicPhenos.csv`, `500FG_basicPhenos.csv` | 500FG phenotypes (age) |
| `cytokine_filter_transposed.rds`, `NPX_FG500.RDS`, `metabolite_trans.rds` | Cytokine responses, Olink proteomics, metabolomics |
| `500FG_inverse_rank_normalized_cellcounts.txt`, `500FG_cellcounts_info.txt` | Immune cell counts and cell population names |
| `500FG_microbiome_pathways.txt` | Gut microbiome pathways |
| `Mvalue_trans.rds` | DNA methylation M values (samples × CpGs) |
| `300bcg/mvalue_v1_trans.rds` | 300BCG DNA methylation M values |
| `300bcg/500fg_tab_common_feat_replaced.csv`, `300bcg/300bcg_tab_common_feat_replaced.csv` | Non-methylation features measured in both cohorts |
| `300bcg/basicPhenos_tab_common_500fg_ordered.csv`, `300bcg/basicPhenos_common_300bcg_feat_replaced.csv` | Phenotypes for the shared-feature tables |

Raw data are not included in this repository.

## Dependencies

- R: `caret`, `glmnet`, `dplyr`, `tidyr`, `purrr`, `ggplot2`, `ggpubr`, `tidyverse`, `jsonlite`, `parallel`
- Python: `tabpfn`, `scikit-learn`, `xgboost`, `catboost`, `pandas`, `numpy`, `joblib`
