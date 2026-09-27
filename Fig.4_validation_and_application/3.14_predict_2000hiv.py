# 3.14 Predict age in 2000HIV with the models saved in 3.12 (5 models x 100 iterations)
#      TabPFN keeps missing values (native handling); other models use the training means.
#      Run from the repository root: python validation/3.14_predict_2000hiv.py

import glob
import os
import numpy as np
import pandas as pd
import pyreadr
from joblib import load
from sklearn.metrics import mean_squared_error

# ---- Paths ----
model_dir = "results/tabpfn_300bcg/saved_models"
hiv_dir = "data/2000hiv"
common_feat_csv = f"{hiv_dir}/hiv_aligned_to_feat_500fg_shared.csv"   # features shared with 500FG/300BCG
mvalue_rds = f"{hiv_dir}/Mvalue_trans.rds"
pheno_csv = f"{hiv_dir}/hiv_basicPhenos_common.csv"
output_dir = "results/validation_2000hiv"
os.makedirs(output_dir, exist_ok=True)

# ---- Load data ----
common_feat = pd.read_csv(common_feat_csv, index_col=0)
common_feat.index = common_feat.index.map(str)
common_feat.columns = [str(c).replace("[", "").replace("]", "").replace("<", "") for c in common_feat.columns]

mvalue = pd.DataFrame(next(iter(pyreadr.read_r(mvalue_rds).values())))
mvalue.index = mvalue.index.map(str)
mvalue.columns = [str(c) for c in mvalue.columns]
ids = set(common_feat.index)
if len(set(mvalue.index) & ids) < len(set(mvalue.columns) & ids):   # probes x samples -> samples x probes
    mvalue = mvalue.T
    mvalue.index = mvalue.index.map(str)
mvalue.columns = [c.split("_")[0] for c in mvalue.columns]            # probe ID without suffix
mvalue = mvalue.T.groupby(level=0).mean().T                          # average duplicated probes

pheno = pd.read_csv(pheno_csv)
age = pd.to_numeric(pheno["age"], errors="coerce")
age.index = pheno["Unnamed: 0"].astype(str)                          # R row names = sample IDs

samples = [s for s in common_feat.index if s in set(mvalue.index)]
print("2000HIV samples with features and methylation:", len(samples))

def metrics(y_true, y_pred):
    d = pd.concat([y_true.rename("t"), y_pred.rename("p")], axis=1).dropna()
    if d.empty:
        return np.nan, np.nan, np.nan, 0
    mse = float(mean_squared_error(d["t"], d["p"]))
    r2 = float(np.corrcoef(d["t"], d["p"])[0, 1] ** 2) if len(d) > 1 else np.nan
    return r2, np.sqrt(mse), mse, len(d)

# ---- Predict with every saved model ----
metric_rows, pred_rows = [], []
for path in sorted(glob.glob(f"{model_dir}/iter_*_*.joblib")):
    _, iteration, model_name = os.path.basename(path).replace(".joblib", "").split("_", 2)
    iteration = int(iteration)
    artifact = load(path)

    X = pd.concat([common_feat.loc[samples], mvalue.loc[samples].reindex(columns=artifact["mvalue_cols"])], axis=1)
    X = X.reindex(columns=artifact["feature_columns"])
    if "tabpfn" not in model_name.lower():
        X = X.fillna(pd.Series(artifact["impute_mean"]))
    X = X.astype(np.float64)

    pred = np.asarray(artifact["model"].predict(X.values), dtype=np.float64)
    y_pred = pd.Series(pred, index=samples, dtype=float)
    r2, rmse, mse, n = metrics(age, y_pred)
    print(f"{model_name} iter {iteration:03d}: N={n}, R2={r2:.4f}, RMSE={rmse:.4f}")

    metric_rows.append({"Model": model_name, "Iteration": iteration, "N_eval": n, "R2": r2, "RMSE": rmse, "MSE": mse})
    pred_rows += [{"Model": model_name, "Iteration": iteration, "sample_id": s,
                   "predicted_age": float(p), "true_age": float(age.get(s, np.nan))} for s, p in zip(samples, pred)]

# ---- Save ----
pd.DataFrame(pred_rows).sort_values(["Model", "Iteration", "sample_id"]).to_csv(
    f"{output_dir}/HIV_all_models_all_iterations_predicted_age.csv", index=False)
pd.DataFrame(metric_rows).sort_values(["Model", "Iteration"]).to_csv(
    f"{output_dir}/HIV_all_models_all_iterations_performance.csv", index=False)
