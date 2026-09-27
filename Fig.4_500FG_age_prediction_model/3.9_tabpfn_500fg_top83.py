# 3.9 Multi-omics age prediction in 500FG with TabPFN, ElasticNet, RandomForest, XGBoost and CatBoost
#     Features: top 83 per layer from 3.8; same 100 fixed splits as 3.4.

import os
import json
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from catboost import CatBoostRegressor
from tabpfn import TabPFNRegressor

# ---- Paths ----
data_dir = "data"
top_dir = "results/fdr_common/FDR_top83"
output_dir = "results/tabpfn_500fg"
os.makedirs(output_dir, exist_ok=True)

# ---- Load data (samples in rows) ----
cytokine = pd.read_csv(f"{data_dir}/cytokine_filter_transposed.csv", index_col=0)
olink = pd.read_csv(f"{data_dir}/olink_NPX_FG500.csv", index_col=0)
metabolite = pd.read_csv(f"{data_dir}/metabolite_trans.csv", index_col=0)
mvalue = pd.read_csv(f"{data_dir}/Mvalue_trans.csv", index_col=0)
cellcounts = pd.read_csv(f"{data_dir}/500FG_inverse_rank_normalized_cellcounts.txt", sep=" ", index_col=0)
microbiome = pd.read_csv(f"{data_dir}/500FG_microbiome_pathways.txt", sep=" ", index_col=0)
basicPheno = pd.read_csv(f"{data_dir}/Age_group_basicPhenos.csv", index_col="ID_500fg")
basicPheno = basicPheno.drop(columns=[c for c in ["X", "Unnamed: 0"] if c in basicPheno.columns])

# ---- Common samples and combined table ----
layers = [cytokine, olink, metabolite, cellcounts, microbiome, mvalue]
common_samples = set(cytokine.index)
for df in layers[1:] + [basicPheno]:
    common_samples &= set(df.index)
common_samples = sorted(common_samples)

combined_table = pd.concat([df.loc[common_samples] for df in layers], axis=1)
combined_table.index = common_samples
combined_table.columns = [str(c).replace("[", "").replace("]", "").replace("<", "") for c in combined_table.columns]
y = basicPheno.loc[common_samples, "Age"].values
print("Combined table:", combined_table.shape)

# ---- Models ----
models = [
    ("ElasticNet", GridSearchCV(
        make_pipeline(StandardScaler(), ElasticNet(random_state=42, max_iter=100000, tol=1e-4)),
        param_grid={"elasticnet__alpha": [0.1, 0.5, 1.0],
                    "elasticnet__l1_ratio": [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]},
        cv=5, scoring="r2", n_jobs=7, verbose=0)),
    ("TabPFN", TabPFNRegressor(random_state=42)),
    ("RandomForest", RandomForestRegressor(random_state=42)),
    ("XGBoost", XGBRegressor(random_state=42)),
    ("CatBoost", CatBoostRegressor(random_state=42, verbose=0)),
]

def evaluate(model, X_train, X_test, y_train, y_test):
    model.fit(X_train, y_train)
    out = {}
    for name, X, yy in [("train", X_train, y_train), ("test", X_test, y_test)]:
        pred = model.predict(X)
        out[f"{name}_R2"] = np.corrcoef(yy, pred)[0, 1] ** 2
        out[f"{name}_MSE"] = mean_squared_error(yy, pred)
        out[f"{name}_RMSE"] = np.sqrt(out[f"{name}_MSE"])
    return out

layer_dirs = ["Cytokine", "Proteomics", "Metabolite", "Cellcounts", "Microbiome", "Methylation_merged"]
results = {name: [] for name, _ in models}

# ---- 100 iterations ----
for i in range(1, 101):
    infos = []
    for layer in layer_dirs:
        with open(f"{top_dir}/{layer}_spearman_preparation_top83/iteration_{i}_fdr_info.json") as f:
            infos.append(json.load(f))

    # Same split for all layers (1-based in R -> 0-based)
    splits = [np.sort(np.array(info["split_indices"], dtype=int).ravel() - 1) for info in infos]
    if any(not np.array_equal(s, splits[0]) for s in splits[1:]):
        raise ValueError(f"Split indices mismatch in iteration {i}")
    train_mask = np.zeros(len(combined_table), dtype=bool)
    train_mask[splits[0]] = True

    # Union of selected features (layer order kept), only those present in the table
    selected = list(dict.fromkeys(f for info in infos for f in info["selected_features"] if f != "Age"))
    available = [f for f in selected if f in combined_table.columns]
    print(f"Iteration {i}: {len(selected)} selected, {len(available)} available")

    X = combined_table[available]
    X_train, X_test = X.iloc[train_mask].values, X.iloc[~train_mask].values
    y_train, y_test = y[train_mask], y[~train_mask]

    for name, model in models:
        results[name].append(evaluate(model, X_train, X_test, y_train, y_test))

# ---- Save ----
rows = [{"Model": name, "Iteration": i,
         "Train_R2": r["train_R2"], "Test_R2": r["test_R2"],
         "Train_RMSE": r["train_RMSE"], "Test_RMSE": r["test_RMSE"],
         "Train_MSE": r["train_MSE"], "Test_MSE": r["test_MSE"]}
        for name, res in results.items() for i, r in enumerate(res, 1)]
pd.DataFrame(rows).to_csv(f"{output_dir}/model_results_spear_top83.csv", index=False)

for name, res in results.items():
    print(f"{name}: train R2 {np.mean([r['train_R2'] for r in res]):.4f}, "
          f"test R2 {np.mean([r['test_R2'] for r in res]):.4f} ± {np.std([r['test_R2'] for r in res]):.4f}")
