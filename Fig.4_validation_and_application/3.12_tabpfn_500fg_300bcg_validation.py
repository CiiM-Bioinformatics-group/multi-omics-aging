# 3.12 Train on 500FG, validate externally in 300BCG
#      Features shared by both cohorts + top 235 shared CpGs per split (3.11); same 100 splits as 3.4.
#      Saves per-iteration models (for application to other cohorts), metrics and mean predicted ages.
#      Run from the repository root: python validation/3.12_tabpfn_500fg_300bcg_validation.py

import os
import json
import numpy as np
import pandas as pd
from collections import defaultdict
from joblib import dump
from sklearn.base import clone
from sklearn.metrics import mean_squared_error
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from catboost import CatBoostRegressor
from tabpfn import TabPFNRegressor

# ---- Paths ----
data_dir = "data"
bcg_dir = f"{data_dir}/300bcg"
top235_dir = "results/fdr_shared_top235/Methylation_FDR_shared_top235"
output_dir = "results/tabpfn_300bcg"
model_dir = f"{output_dir}/saved_models"
os.makedirs(model_dir, exist_ok=True)

# ---- Load data ----
common_feat_500fg = pd.read_csv(f"{bcg_dir}/500fg_tab_common_feat_replaced.csv", index_col=0)
common_feat_300bcg = pd.read_csv(f"{bcg_dir}/300bcg_tab_common_feat_replaced.csv", index_col=0)
basicPhenos_300bcg = pd.read_csv(f"{bcg_dir}/basicPhenos_common_300bcg_feat_replaced.csv", index_col="PatientID")
basicPhenos_500fg = pd.read_csv(f"{bcg_dir}/basicPhenos_tab_common_500fg_ordered.csv", index_col="ID_500fg")
mvalue_500fg = pd.read_csv(f"{data_dir}/Mvalue_trans.csv", index_col=0)
mvalue_300bcg = pd.read_csv(f"{bcg_dir}/Mvalue_300BCG_common.csv", index_col=0)

# 500FG probe names: average duplicated probes, keep the part before "_"
mvalue_500fg = mvalue_500fg.T.groupby(level=0).mean().T
mvalue_500fg.columns = [c.split("_")[0] for c in mvalue_500fg.columns]

# ---- Samples with both shared features and methylation (order of the feature tables) ----
common_feat_300bcg.index = common_feat_300bcg.index.map(str)
samples_500fg = [s for s in common_feat_500fg.index if s in set(mvalue_500fg.index)]
common_feat_500fg = common_feat_500fg.loc[samples_500fg]
mvalue_500fg = mvalue_500fg.loc[samples_500fg]

mvalue_300bcg.index = mvalue_300bcg.index.map(lambda s: str(int(s[6:9])))   # sample ID -> patient number
samples_300bcg = [s for s in common_feat_300bcg.index if s in set(mvalue_300bcg.index)]
common_feat_300bcg = common_feat_300bcg.loc[samples_300bcg]
mvalue_300bcg = mvalue_300bcg.loc[samples_300bcg]

clean = lambda cols: [str(c).replace("[", "").replace("]", "").replace("<", "") for c in cols]
common_feat_500fg.columns = clean(common_feat_500fg.columns)
common_feat_300bcg.columns = clean(common_feat_300bcg.columns)
print("500FG:", common_feat_500fg.shape, "| 300BCG:", common_feat_300bcg.shape)

# ---- Models ----
models = [
    ("TabPFN", TabPFNRegressor(random_state=42)),
    ("ElasticNet", Pipeline([
        ("scaler", StandardScaler()),
        ("elasticnet", ElasticNetCV(l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 1.0], alphas=None,
                                    cv=5, max_iter=100000, random_state=42))])),
    ("RandomForest", RandomForestRegressor(random_state=42)),
    ("XGBoost", XGBRegressor(random_state=42)),
    ("CatBoost", CatBoostRegressor(random_state=42, verbose=0)),
]

def evaluate(model, sets):
    """Fit on the training set; R2 (squared Pearson r), RMSE, MSE and predictions per set."""
    X_train, y_train = sets["train"]
    model.fit(X_train, y_train)
    out = {}
    for name, (X, yy) in sets.items():
        pred = np.asarray(model.predict(X), dtype=np.float64)
        out[f"{name}_R2"] = np.corrcoef(yy, pred)[0, 1] ** 2
        out[f"{name}_MSE"] = mean_squared_error(yy, pred)
        out[f"{name}_RMSE"] = np.sqrt(out[f"{name}_MSE"])
        out[f"{name}_pred"] = pred
    return out

results = {name: [] for name, _ in models}
preds = {s: {name: defaultdict(list) for name, _ in models} for s in ["train", "test", "external"]}
y_500fg = basicPhenos_500fg.loc[common_feat_500fg.index, "Age"]
y_external = basicPhenos_300bcg.loc[[int(s) for s in samples_300bcg], "Age"].values.astype(np.float64)

# ---- 100 iterations ----
for i in range(1, 101):
    with open(f"{top235_dir}/iteration_{i}_fdr_info.json") as f:
        info = json.load(f)
    split_indices = np.array(info["split_indices"], dtype=int).ravel() - 1   # R 1-based -> 0-based

    # Methylation features: probe ID (first 10 characters), present in 500FG
    mvalue_cols = [c for c in dict.fromkeys(f[:10] for f in info["selected_features"]) if c in mvalue_500fg.columns]

    X_merged = pd.concat([common_feat_500fg, mvalue_500fg[mvalue_cols]], axis=1)
    if X_merged.shape[1] > 500:   # TabPFN feature limit
        X_merged = X_merged.iloc[:, :500]
    X_merged = X_merged.fillna(X_merged.mean())

    train_mask = np.zeros(len(X_merged), dtype=bool)
    train_mask[split_indices] = True

    X_external = pd.concat([common_feat_300bcg, mvalue_300bcg.reindex(columns=mvalue_cols)], axis=1)
    X_external = X_external.reindex(columns=X_merged.columns).fillna(X_merged.mean())

    sets = {
        "train": (X_merged.iloc[train_mask].values.astype(np.float64), y_500fg.iloc[train_mask].values.astype(np.float64)),
        "test": (X_merged.iloc[~train_mask].values.astype(np.float64), y_500fg.iloc[~train_mask].values.astype(np.float64)),
        "external": (X_external.values.astype(np.float64), y_external),
    }
    sample_ids = {"train": X_merged.index[train_mask], "test": X_merged.index[~train_mask], "external": X_external.index}
    print(f"Iteration {i}: {X_merged.shape[1]} features ({len(mvalue_cols)} CpGs)")

    for name, template in models:
        model = clone(template)
        r = evaluate(model, sets)
        dump({"model": model, "feature_columns": X_merged.columns.tolist(),
              "mvalue_cols": mvalue_cols, "impute_mean": X_merged.mean().to_dict()},
             f"{model_dir}/iter_{i:03d}_{name}.joblib")
        for s in preds:
            for sid, p in zip(sample_ids[s].astype(str), r[f"{s}_pred"]):
                preds[s][name][sid].append(p)
        results[name].append({k: v for k, v in r.items() if not k.endswith("_pred")})

# ---- Save metrics ----
rows = [{"Model": name, "Iteration": i,
         "Train_R2": r["train_R2"], "Test_R2": r["test_R2"], "External_R2": r["external_R2"],
         "Train_RMSE": r["train_RMSE"], "Test_RMSE": r["test_RMSE"], "External_RMSE": r["external_RMSE"],
         "Train_MSE": r["train_MSE"], "Test_MSE": r["test_MSE"], "External_MSE": r["external_MSE"]}
        for name, res in results.items() for i, r in enumerate(res, 1)]
pd.DataFrame(rows).to_csv(f"{output_dir}/model_results_300bcg_external_with_methy_share_feattop235.csv", index=False)

# ---- Save mean predicted age per sample ----
file_tag = {"train": "500fg_train", "test": "500fg_test", "external": "300bcg"}
order = {"train": common_feat_500fg.index, "test": common_feat_500fg.index, "external": pd.Index(samples_300bcg)}
for s, by_model in preds.items():
    for name, d in by_model.items():
        rows = [{"sample_id": sid, "predicted_age": float(np.mean(d[sid]))} for sid in order[s].astype(str) if d.get(sid)]
        pd.DataFrame(rows).to_csv(f"{output_dir}/{name}_{file_tag[s]}_predicted_age_ensemble.csv", index=False)

for name, res in results.items():
    print(f"{name}: test R2 {np.mean([r['test_R2'] for r in res]):.4f}, "
          f"external R2 {np.mean([r['external_R2'] for r in res]):.4f}")
