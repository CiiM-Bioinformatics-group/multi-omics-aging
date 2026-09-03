"""
Predict HIV age using saved models from 500FG training.

Minimal version for fixed HIV paths:
1) load HIV common features (CSV)
2) load HIV methylation (RDS)
3) methylation harmonization:
   - keep sample IDs as string
   - remove suffix after "_" in probe names
   - if duplicated probe names appear after suffix removal, take mean
4) for each saved model artifact:
   - rebuild external matrix with artifact["mvalue_cols"]
   - align to artifact["feature_columns"]
   - TabPFN: keep missing values as NA for native handling
   - other models: impute missing values by artifact["impute_mean"]
   - predict and save per-iteration CSV
5) save per-iteration performance on HIV (R2/RMSE/MAE)
"""

import glob
import os
import time
from typing import Optional

import numpy as np
import pandas as pd
from joblib import load
from sklearn.metrics import mean_squared_error


MODEL_DIR = "/vol/projects/yzhang/300BCG/output/05_foundation_model/with_methylation/saved_models"
HIV_COMMON_FEAT_CSV = "/vol/projects/yzhang/2000HIV/multi_omics_model/output/hiv_aligned_to_feat_500fg_shared.csv"
HIV_MVALUE_RDS = "/vol/projects/yzhang/2000HIV/methylation/output/Mvalue_trans.rds"
HIV_PHENO_CSV = "/vol/projects/yzhang/2000HIV/multi_omics_model/output/hiv_basicPhenos_common.csv"
# `X` is the original R rownames (for example, EMC001) exported by read.csv().
HIV_PHENO_ID_COL = "Unnamed: 0"
HIV_PHENO_AGE_COL = "age"
OUTPUT_DIR = "/vol/projects/yzhang/2000HIV/multi_omics_model/output/external_prediction_with_saved_models"
MODEL_NAME_FILTER = None  # Run all models in saved_models

# Auto-detect methylation orientation by overlap with common feature sample IDs.
# Set to True/False only if you want to force behavior.
MVALUE_SAMPLES_ARE_ROWS = None
PREVIEW_ONLY = False


def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def print_preview(name: str, df: pd.DataFrame) -> None:
    print(f"\n[{name}] shape: {df.shape}", flush=True)
    print(f"{name} index examples: {list(df.index[:5])}", flush=True)
    print(f"{name} column examples: {list(df.columns[:5])}", flush=True)
    print(f"{name} 3x5 preview:", flush=True)
    print(df.iloc[:3, :5], flush=True)


def clean_feature_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).replace("[", "").replace("]", "").replace("<", "") for c in df.columns]
    return df


def strip_suffix(name: str) -> str:
    return str(name).split("_")[0]


def load_HIV_common_feat(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.map(str)
    return clean_feature_columns(df)


def load_HIV_mvalue_rds(
    path: str,
    common_sample_ids: pd.Index,
    samples_are_rows: Optional[bool] = None,
) -> pd.DataFrame:
    try:
        import pyreadr
    except ImportError as exc:
        raise ImportError("Please install pyreadr first: pip install pyreadr") from exc

    result = pyreadr.read_r(path)
    if not result:
        raise ValueError(f"No object found in RDS: {path}")

    df = pd.DataFrame(next(iter(result.values())))
    df.index = df.index.map(str)
    df.columns = [str(c) for c in df.columns]
    print_preview("mvalue_raw", df)

    common_id_set = set(common_sample_ids.map(str))
    overlap_rows = len(set(df.index) & common_id_set)
    overlap_cols = len(set(df.columns) & common_id_set)

    if samples_are_rows is None:
        # Auto choose axis with larger overlap as sample axis.
        samples_are_rows = overlap_rows >= overlap_cols

    if samples_are_rows:
        # already samples x probes
        pass
    else:
        # probes x samples -> transpose to samples x probes
        df = df.T
        df.index = df.index.map(str)
    print_preview("mvalue_after_orientation", df)

    print(
        f"Methylation ID overlap check: rows={overlap_rows}, cols={overlap_cols}; "
        f"samples_are_rows={samples_are_rows}"
    , flush=True)

    # Harmonize probe names: remove suffix and merge duplicates by mean
    df.columns = [strip_suffix(c) for c in df.columns]
    df = df.T.groupby(level=0).mean().T
    print_preview("mvalue_after_probe_harmonization", df)

    return df


def is_tabpfn_model(model_name: str) -> bool:
    return "tabpfn" in str(model_name).lower()


def load_HIV_age(path: str) -> pd.Series:
    # Use the exported R rownames as sample IDs, rather than relying on row order.
    pheno = pd.read_csv(path)
    required_cols = {HIV_PHENO_ID_COL, HIV_PHENO_AGE_COL}
    missing_cols = required_cols - set(pheno.columns)
    if missing_cols:
        raise ValueError(
            f"Required phenotype column(s) {sorted(missing_cols)} not found: {path}. "
            f"Available columns: {list(pheno.columns)}"
        )

    raw_sample_ids = pheno[HIV_PHENO_ID_COL]
    if raw_sample_ids.isna().any() or (raw_sample_ids.astype(str).str.strip() == "").any():
        raise ValueError(f"Phenotype ID column '{HIV_PHENO_ID_COL}' contains missing IDs")
    sample_ids = raw_sample_ids.astype(str)
    if sample_ids.duplicated().any():
        duplicates = sample_ids[sample_ids.duplicated()].unique().tolist()
        raise ValueError(f"Phenotype ID column '{HIV_PHENO_ID_COL}' has duplicate IDs: {duplicates[:5]}")

    age = pd.to_numeric(pheno[HIV_PHENO_AGE_COL], errors="coerce")
    age.index = sample_ids
    return age


def select_common_samples(common_feat: pd.DataFrame, mvalue: pd.DataFrame) -> list[str]:
    shared = set(common_feat.index) & set(mvalue.index)
    ordered_shared = [sid for sid in common_feat.index if sid in shared]
    if not ordered_shared:
        raise ValueError(
            "No shared samples between common feature and methylation.\n"
            f"common examples: {list(common_feat.index[:5])}\n"
            f"mvalue examples: {list(mvalue.index[:5])}"
        )
    return ordered_shared


def build_external_matrix(
    common_feat: pd.DataFrame,
    mvalue: pd.DataFrame,
    feature_columns: list[str],
    mvalue_cols: list[str],
    impute_mean: dict,
    keep_na: bool = False,
) -> tuple[pd.DataFrame, pd.Index, int]:
    shared = select_common_samples(common_feat, mvalue)

    x_common = common_feat.loc[shared]
    x_mvalue = mvalue.loc[shared].reindex(columns=mvalue_cols)

    x_external = pd.concat([x_common, x_mvalue], axis=1)
    x_external = x_external.reindex(columns=feature_columns)
    if not keep_na:
        x_external = x_external.fillna(pd.Series(impute_mean))
    x_external = x_external.astype(np.float64)
    n_missing = int(x_external.isna().sum().sum())

    return x_external, pd.Index(shared), n_missing


def parse_model_file(path: str) -> tuple[int, str]:
    # Example: iter_001_TabPFN.joblib
    base = os.path.basename(path).replace(".joblib", "")
    parts = base.split("_", 2)
    if len(parts) != 3:
        raise ValueError(f"Unexpected model filename: {base}")
    return int(parts[1]), parts[2]


def compute_metrics(y_true: pd.Series, y_pred: pd.Series) -> tuple[float, float, float, int]:
    aligned = pd.concat([y_true.rename("y_true"), y_pred.rename("y_pred")], axis=1).dropna()
    if aligned.empty:
        return np.nan, np.nan, np.nan, 0

    y_true_arr = aligned["y_true"].to_numpy(dtype=np.float64)
    y_pred_arr = aligned["y_pred"].to_numpy(dtype=np.float64)

    mse = float(mean_squared_error(y_true_arr, y_pred_arr))
    rmse = float(np.sqrt(mse))
    if len(y_true_arr) < 2:
        r2 = np.nan
    else:
        r = np.corrcoef(y_true_arr, y_pred_arr)[0, 1]
        r2 = float(r ** 2) if np.isfinite(r) else np.nan

    return r2, rmse, mse, int(len(aligned))


def main() -> None:
    log(f"PID={os.getpid()} started.")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    log("Loading HIV data...")
    common_feat = load_HIV_common_feat(HIV_COMMON_FEAT_CSV)
    print_preview("common_feat", common_feat)
    log("Loaded common feature CSV")
    mvalue = load_HIV_mvalue_rds(
        HIV_MVALUE_RDS,
        common_sample_ids=common_feat.index,
        samples_are_rows=MVALUE_SAMPLES_ARE_ROWS,
    )
    log("Loaded and processed methylation RDS")
    y_HIV = load_HIV_age(HIV_PHENO_CSV)
    log("Loaded phenotype")
    print(f"\n[phenotype] shape: {y_HIV.shape}", flush=True)
    print(f"phenotype index examples: {list(y_HIV.index[:5])}", flush=True)
    print("phenotype Age examples:", flush=True)
    print(y_HIV.head(5), flush=True)

    log(f"HIV common feature shape: {common_feat.shape}")
    log(f"HIV methylation shape (processed): {mvalue.shape}")
    log(f"HIV phenotype valid Age count: {int(y_HIV.notna().sum())}")

    if PREVIEW_ONLY:
        log("PREVIEW_ONLY=True, stop after dataframe previews.")
        return

    model_files = sorted(glob.glob(os.path.join(MODEL_DIR, "iter_*_*.joblib")))
    if not model_files:
        raise FileNotFoundError(f"No saved model files found in {MODEL_DIR}")

    parsed = []
    for p in model_files:
        it, mn = parse_model_file(p)
        parsed.append((p, it, mn))

    if MODEL_NAME_FILTER is not None:
        model_name_set = set(MODEL_NAME_FILTER)
        parsed = [x for x in parsed if x[2] in model_name_set]

    if not parsed:
        raise ValueError(f"No model files left after MODEL_NAME_FILTER={MODEL_NAME_FILTER}")

    unique_models = sorted({x[2] for x in parsed})
    unique_iterations = sorted({x[1] for x in parsed})
    log(
        f"Model files selected: {len(parsed)} | "
        f"models={unique_models} ({len(unique_models)}) | "
        f"iterations={len(unique_iterations)} ({min(unique_iterations)}-{max(unique_iterations)})"
    )
    metrics_rows = []
    prediction_rows = []
    t0 = time.time()

    for i, (model_path, iteration, model_name) in enumerate(parsed, start=1):
        log(f"[{i}/{len(parsed)}] Start {model_name} iter {iteration:03d}")
        artifact = load(model_path)
        keep_na = is_tabpfn_model(model_name)

        x_external, sample_ids, n_missing = build_external_matrix(
            common_feat=common_feat,
            mvalue=mvalue,
            feature_columns=artifact["feature_columns"],
            mvalue_cols=artifact["mvalue_cols"],
            impute_mean=artifact["impute_mean"],
            keep_na=keep_na,
        )
        if keep_na:
            log(
                f"{model_name} iter {iteration:03d}: kept {n_missing} NA values "
                "for TabPFN native missing-value handling"
            )
        elif n_missing:
            log(
                f"{model_name} iter {iteration:03d}: {n_missing} NA remain after mean imputation"
            )

        preds = np.asarray(artifact["model"].predict(x_external.values), dtype=np.float64)
        y_pred = pd.Series(preds, index=sample_ids.astype(str), dtype=float)
        r2, rmse, mse, n_eval = compute_metrics(y_HIV, y_pred)

        for sid, pred in zip(sample_ids.astype(str), preds):
            prediction_rows.append(
                {
                    "Model": model_name,
                    "Iteration": iteration,
                    "sample_id": sid,
                    "predicted_age": float(pred),
                    "true_age": float(y_HIV.get(sid, np.nan)),
                }
            )
        elapsed_min = (time.time() - t0) / 60.0
        log(
            f"[{i}/{len(parsed)}] Done {model_name} iter {iteration:03d} | "
            f"N={n_eval}, R2={r2:.4f}, RMSE={rmse:.4f}, MSE={mse:.4f}, "
            f"elapsed={elapsed_min:.1f}m, PID={os.getpid()}"
        )

        metrics_rows.append(
            {
                "Model": model_name,
                "Iteration": iteration,
                "N_eval": n_eval,
                "R2": r2,
                "RMSE": rmse,
                "MSE": mse,
            }
        )

    metrics_df = pd.DataFrame(metrics_rows).sort_values(["Model", "Iteration"])
    pred_df = pd.DataFrame(prediction_rows).sort_values(["Model", "Iteration", "sample_id"])

    # 1) all predicted ages (all models, all iterations, all samples)
    pred_df.to_csv(
        os.path.join(OUTPUT_DIR, "HIV_all_models_all_iterations_predicted_age.csv"),
        index=False,
    )

    # 2) all performance rows (all models, all iterations)
    metrics_df.to_csv(
        os.path.join(OUTPUT_DIR, "HIV_all_models_all_iterations_performance.csv"),
        index=False,
    )

    log(f"Done. Output saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
