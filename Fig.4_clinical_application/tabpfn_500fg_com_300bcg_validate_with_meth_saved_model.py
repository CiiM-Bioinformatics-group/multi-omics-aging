import os
import torch
import json
# Setup Imports
import pandas as pd
import numpy as np

from sklearn.datasets import load_breast_cancer, load_diabetes, load_iris
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.metrics import (
    accuracy_score,
    mean_absolute_error,
    mean_squared_error,
    root_mean_squared_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RepeatedKFold 

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.inspection import DecisionBoundaryDisplay

from sklearn.datasets import fetch_openml
from sklearn.preprocessing import LabelEncoder
from IPython.display import display, Markdown, Latex

# Baseline Imports
from xgboost import XGBClassifier, XGBRegressor
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from catboost import CatBoostClassifier, CatBoostRegressor

from tabpfn import TabPFNClassifier, TabPFNRegressor
from tabpfn_extensions.post_hoc_ensembles.sklearn_interface import AutoTabPFNClassifier, AutoTabPFNRegressor

from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, Lasso
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from catboost import CatBoostRegressor
from sklearn.model_selection import cross_val_score
from joblib import dump, load
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from sklearn.pipeline import make_pipeline
import seaborn as sns
    
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import glob
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import ElasticNetCV 
from sklearn.base import clone
from collections import defaultdict

###500fg
# Read csv files
common_feat_500fg = pd.read_csv("/vol/projects/yzhang/300BCG/input/foundation_model/500fg_tab_common_feat_replaced.csv",index_col=0)
common_feat_300bcg = pd.read_csv("/vol/projects/yzhang/300BCG/input/foundation_model/300bcg_tab_common_feat_replaced.csv",index_col=0)
basicPhenos_300bcg = pd.read_csv("/vol/projects/yzhang/300BCG/input/foundation_model/basicPhenos_common_300bcg_feat_replaced.csv",index_col='PatientID')
basicPhenos_500fg = pd.read_csv('/vol/projects/yzhang/300BCG/input/foundation_model/basicPhenos_tab_common_500fg_ordered.csv', index_col='ID_500fg')
mvalue_500fg = pd.read_csv("/vol/projects/yzhang/500FG_aging/input/methylation/Mvalue_trans.csv",index_col=0)
mvalue_300bcg = pd.read_csv("/vol/projects/yzhang/300BCG/input/methylation/Mvalue_300BCG_common.csv",index_col=0)
common = pd.read_csv("/vol/projects/yzhang/300BCG/input/methylation/common_cols.csv",index_col='x')

#mvalue_500fg.iloc[:, :20].to_csv("/vol/projects/yzhang/500FG_aging/input/methylation/Mvalue_trans_first20cols.csv")
#mvalue_300bcg.iloc[:, :20].to_csv("/vol/projects/yzhang/300BCG/input/methylation/Mvalue_300BCG_common_first20cols.csv")

# Use first column as index
print("common_feat_300bcg index (first 10):", common_feat_300bcg.index)
print("basicPhenos_300bcg index (first 10):", basicPhenos_300bcg.index)
# Print head of each dataframe to verify
# Display data summary
def display_data_info(df, name):
    print(f"\n{name} shape: {df.shape}")
    print(f"\n{name} head (first 4 columns):")
    print(df.iloc[:, :4].head())  # show first 4 columns only
    print("-" * 50)  # separator

# Display info for each dataframe
display_data_info(common_feat_500fg, "common_feat_500fg")
display_data_info(basicPhenos_500fg, "basicPhenos_500fg")
display_data_info(common_feat_300bcg, "common_feat_300bcg")
display_data_info(basicPhenos_300bcg, "basicPhenos_300bcg")
display_data_info(mvalue_500fg, "mvalue_500fg")
display_data_info(mvalue_300bcg, "mvalue_300bcg")
display_data_info(common, "common_col")

# Harmonize methylation probe names across cohorts
def strip_suffix(col):
    return col.split('_')[0]  # keep prefix before underscore

# Preprocess 500FG methylation matrix
mvalue_500fg = mvalue_500fg.T.groupby(level=0).mean().T
mvalue_500fg.columns = [strip_suffix(col) for col in mvalue_500fg.columns]
mvalue_500fg.to_csv("/vol/projects/yzhang/500FG_aging/input/methylation/500fg_Mvalue_trans_no_suffix.csv")

# Optional: fill NaN with column means
#mvalue_500fg = mvalue_500fg.fillna(mvalue_500fg.mean())
#mvalue_300bcg = mvalue_300bcg.fillna(mvalue_300bcg.mean())
#mvalue_300bcg.to_csv("/vol/projects/yzhang/300BCG/input/methylation/Mvalue_300BCG_common_filled.csv")

###select common samples with methylation 
# 1) Find shared samples (index intersection)
common_feat_300bcg.index = common_feat_300bcg.index.map(str)
shared_samples_500fg = list(set(common_feat_500fg.index) & set(mvalue_500fg.index))
shared_samples_500fg_sorted = [idx for idx in common_feat_500fg.index if idx in shared_samples_500fg]

# 2) Subset and sort by template order
common_feat_500fg_common = common_feat_500fg.loc[shared_samples_500fg_sorted]
mvalue_500fg_common = mvalue_500fg.loc[shared_samples_500fg_sorted]

def extract_number(idx):
    # chars 7-9 (Python index 6:9), strip leading zeros
    return str(int(idx[6:9]))

# Normalize index type
mvalue_300bcg.index = mvalue_300bcg.index.map(extract_number)

# Sanity-check indices
print(mvalue_300bcg.index[:5])
print("common_feat_300bcg index:", list(common_feat_300bcg.index[:10]))
print("mvalue_300bcg index:", list(mvalue_300bcg.index[:10]))

# External cohort (300BCG): shared samples
shared_samples_300bcg = list(set(common_feat_300bcg.index) & set(mvalue_300bcg.index))
shared_samples_300bcg_sorted = [idx for idx in common_feat_300bcg.index if idx in shared_samples_300bcg]
print("shared_samples_300bcg (first 20):", shared_samples_300bcg[:20])
print("Intersection size:", len(shared_samples_300bcg))

common_feat_300bcg_common = common_feat_300bcg.loc[shared_samples_300bcg_sorted]
mvalue_300bcg_common = mvalue_300bcg.loc[shared_samples_300bcg_sorted]


print("mvalue_500fg_common shape:", mvalue_500fg_common.shape)
print("mvalue_500fg_common head:")
print(mvalue_500fg_common.head())

print("mvalue_300bcg_common shape:", mvalue_300bcg_common.shape)
print("mvalue_300bcg_common head:")
print(mvalue_300bcg_common.head())

common_feat_500fg_common.columns = [str(col).replace('[', '').replace(']', '').replace('<', '') for col in common_feat_500fg_common.columns]
print("\nAfter cleaning feature names:")
print("Column examples:", common_feat_500fg_common.columns[:5].tolist())

common_feat_300bcg_common.columns = [str(col).replace('[', '').replace(']', '').replace('<', '') for col in common_feat_300bcg_common.columns]
print("\nAfter cleaning feature names:")
print("Column examples:", common_feat_300bcg_common.columns[:5].tolist())



# Define models
models = [
    ('TabPFN', TabPFNRegressor(random_state=42)),
    ('ElasticNet', Pipeline([
        ('scaler', StandardScaler()),
        ('elasticnet', ElasticNetCV(
           l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 1.0],
            alphas=None,  # auto-select alpha
            cv=5,
            max_iter=100000,
            random_state=42,
            verbose=1
       ))
    ])),    
    ('RandomForest', RandomForestRegressor(random_state=42)),
    ('XGBoost', XGBRegressor(random_state=42)),
    ('CatBoost', CatBoostRegressor(random_state=42, verbose=0))
]

# Initial X/y (sanity check)
X = common_feat_500fg_common.values  # convert to numpy array
y = basicPhenos_500fg['Age'].values  # convert to numpy array

# Check X/y shapes
print("X shape:", X.shape)
print("y shape:", y.shape)

#X_external = common_feat_300bcg_common.values
#y_external = basicPhenos_300bcg['Age'].values

results = {name: [] for name, _ in models}
successful_iterations = 0
failed_iterations = 0

#  Initialize arrays to store results
n_iterations = 100
train_r2_values = {name: [] for name, _ in models}
val_r2_values = {name: [] for name, _ in models}
train_preds_by_model = {name: defaultdict(list) for name, _ in models}
test_preds_by_model = {name: defaultdict(list) for name, _ in models}
external_preds_by_model = {name: defaultdict(list) for name, _ in models}
external_index_for_save = common_feat_300bcg_common.index.intersection(mvalue_300bcg_common.index)

output_dir = "/vol/projects/yzhang/300BCG/output/05_foundation_model/with_methylation"
model_dir = os.path.join(output_dir, "saved_models")
os.makedirs(model_dir, exist_ok=True)

def read_json_file(path):
    with open(path, 'r') as f:
        return json.load(f)
base_path ="/vol/projects/yzhang/300BCG/input/methylation/FDR_shared_top235"        
# Loop through each iteration
for i in range(1, n_iterations + 1):
    try:      
        # Split data into training and validation sets
        #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=i)
        # 1) Load iteration JSON
        mvalue_info = read_json_file(f"{base_path}/Methylation_FDR_shared_top235/iteration_{i}_fdr_info.json")
        split_indices = np.array(mvalue_info['split_indices'], dtype=int) - 1  # 0-based

        # 2) Methylation features from JSON
        #mvalue_features = mvalue_info['selected_features'] if 'selected_features' in mvalue_info else []
        mvalue_features = [f[:10] for f in mvalue_info.get('selected_features', [])]
        unique_mvalue_features = list(dict.fromkeys(mvalue_features))
        mvalue_cols = [f for f in unique_mvalue_features if f in mvalue_500fg_common.columns]
        if len(mvalue_cols) == 0:
            raise ValueError(f"Iteration {i}: no methylation features found in 500FG")

        X_merged = pd.concat([common_feat_500fg_common, mvalue_500fg_common[mvalue_cols]], axis=1)
        if X_merged.shape[1] > 500:
            print(f"Warning: {X_merged.shape[1]} features detected. For TabPFN compatibility, using top 500 features.")
            X_merged = X_merged.iloc[:, :500]
        # Impute missing values after feature selection
        X_merged = X_merged.fillna(X_merged.mean())
        
        # 6. y
        y = basicPhenos_500fg.loc[common_feat_500fg_common.index, 'Age']

        # Build train mask from split_indices
        n_samples = X_merged.shape[0]
        train_mask = np.zeros(n_samples, dtype=bool)
        train_mask[split_indices] = True

        X_train = X_merged.iloc[train_mask]
        y_train = y.iloc[train_mask]
        
        # Test set: complement of train mask in 500FG
        test_mask = ~train_mask  # complement of training mask
        X_test = X_merged.iloc[test_mask]
        y_test = y.iloc[test_mask]

        # Build 300BCG external feature matrix
        external_index = common_feat_300bcg_common.index.intersection(mvalue_300bcg_common.index)
        X_external_common = common_feat_300bcg_common.loc[external_index]
        X_external_mvalue = mvalue_300bcg_common.loc[external_index].reindex(columns=mvalue_cols)
        X_external = pd.concat([X_external_common, X_external_mvalue], axis=1)
        
        # Align external columns to training columns
        X_external = X_external.reindex(columns=X_merged.columns)
        X_external = X_external.fillna(X_merged.mean())
        
        external_index_int = [int(idx) for idx in external_index]
        y_external = basicPhenos_300bcg.loc[external_index_int, 'Age']

        print(f"Iteration {i}:")
        print(f"Training set shape: {X_train.shape}")
        print(f"Test set shape: {X_test.shape}")
        print(f"External validation set shape: {X_external.shape}")
        print(f"Number of methylation features: {len(mvalue_cols)}")
        
        # Verify column order consistency
        print(f"Training columns (first 5): {list(X_train.columns[:5])}")
        print(f"Test columns (first 5): {list(X_test.columns[:5])}")
        print(f"External columns (first 5): {list(X_external.columns[:5])}")
        print(f"Columns match: {list(X_train.columns) == list(X_test.columns) == list(X_external.columns)}")    
        
        def multi_omics_model_train(X_train, X_test, y_train, y_test, model, X_external, y_external):
            """Train and evaluate a model, including external validation"""
            # Use numpy arrays to avoid pandas column-name issues
            X_train_array = X_train.values if hasattr(X_train, 'values') else X_train
            X_test_array = X_test.values if hasattr(X_test, 'values') else X_test
            X_external_array = X_external.values if hasattr(X_external, 'values') else X_external
            
            # Ensure y is numpy array
            y_train_array = y_train.values if hasattr(y_train, 'values') else y_train
            y_test_array = y_test.values if hasattr(y_test, 'values') else y_test
            y_external_array = y_external.values if hasattr(y_external, 'values') else y_external
            
            # Cast to float64
            X_train_array = X_train_array.astype(np.float64)
            X_test_array = X_test_array.astype(np.float64)
            X_external_array = X_external_array.astype(np.float64)
            
            y_train_array = y_train_array.astype(np.float64)
            y_test_array = y_test_array.astype(np.float64)
            y_external_array = y_external_array.astype(np.float64)
            
            if isinstance(model, GridSearchCV):
                model.fit(X_train_array, y_train_array)
                train_r2 = model.best_score_
            else:
                train_r2 = cross_val_score(model, X_train_array, y_train_array, cv=5, scoring='r2', n_jobs=7).mean()
                model.fit(X_train_array, y_train_array)
            
            # Training predictions
            train_pred = model.predict(X_train_array)
            r = np.corrcoef(y_train_array, train_pred)[0, 1]
            train_r2 = r ** 2
            train_rmse = np.sqrt(mean_squared_error(y_train_array, train_pred))
            train_mse = mean_squared_error(y_train_array, train_pred)
            
            # Internal test predictions
            test_pred = model.predict(X_test_array)
            r = np.corrcoef(y_test_array, test_pred)[0, 1]
            test_r2 = r ** 2
            test_rmse = np.sqrt(mean_squared_error(y_test_array, test_pred))
            test_mse = mean_squared_error(y_test_array, test_pred)

            # External validation
            external_pred = model.predict(X_external_array)
            r = np.corrcoef(y_external_array, external_pred)[0, 1]
            external_r2 = r ** 2
            external_rmse = np.sqrt(mean_squared_error(y_external_array, external_pred))
            external_mse = mean_squared_error(y_external_array, external_pred)
            
            return {
                'train_R2': train_r2,
                'test_R2': test_r2,
                'train_RMSE': train_rmse,
                'test_RMSE': test_rmse,
                'train_MSE': train_mse,
                'test_MSE': test_mse,
                'external_R2': external_r2,
                'external_RMSE': external_rmse,
                'external_MSE': external_mse,
                'model': model,
                'train_pred': np.asarray(train_pred, dtype=np.float64),
                'test_pred': np.asarray(test_pred, dtype=np.float64),
                'external_pred': np.asarray(external_pred, dtype=np.float64),
            }
        
        # Train and evaluate each model
        for name, model_template in models:
            try:
                model = clone(model_template)
                result = multi_omics_model_train(X_train, X_test, y_train, y_test, model, X_external, y_external)

                safe_name = name.replace(' ', '_')
                dump(
                    {
                        'model': result['model'],
                        'feature_columns': X_merged.columns.tolist(),
                        'mvalue_cols': mvalue_cols,
                        'impute_mean': X_merged.mean().to_dict(),
                    },
                    os.path.join(model_dir, f'iter_{i:03d}_{safe_name}.joblib'),
                )

                for sid, pred in zip(X_train.index.astype(str), result['train_pred']):
                    train_preds_by_model[name][sid].append(pred)
                for sid, pred in zip(X_test.index.astype(str), result['test_pred']):
                    test_preds_by_model[name][sid].append(pred)
                for sid, pred in zip(external_index.astype(str), result['external_pred']):
                    external_preds_by_model[name][sid].append(pred)

                result_out = {k: v for k, v in result.items() if k not in ('model', 'train_pred', 'test_pred', 'external_pred')}
                results[name].append(result_out)
                print(f"{name} performance:")
                print(f"  Training R²: {result['train_R2']:.4f}")
                print(f"  Test R²: {result['test_R2']:.4f}")
                print(f"  External Validation R²: {result['external_R2']:.4f}")
  
            except Exception as e:
                print(f"Error in model {name}: {str(e)}")
                continue

        successful_iterations += 1

    except Exception as e:
        print(f"Error in iteration {i}: {str(e)}")
        failed_iterations += 1
        continue

# Summary
print("\nModeling completed:")
print(f"Total iterations: {n_iterations}")
print(f"Successful iterations: {successful_iterations}")
print(f"Failed iterations: {failed_iterations}")

# Save outputs
print("\nSaving results...")


# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

def save_pred_ensemble(pred_dict, sample_ids, out_path):
    rows = []
    for sid in sample_ids.astype(str):
        preds = pred_dict.get(sid, [])
        if preds:
            rows.append({'sample_id': sid, 'predicted_age': float(np.mean(preds))})
    pd.DataFrame(rows).to_csv(out_path, index=False)

for model_name, pred_dict in train_preds_by_model.items():
    if not pred_dict:
        continue
    safe = model_name.replace(' ', '_')
    save_pred_ensemble(
        pred_dict,
        common_feat_500fg_common.index,
        os.path.join(output_dir, f'{safe}_500fg_train_predicted_age_ensemble.csv'),
    )

for model_name, pred_dict in test_preds_by_model.items():
    if not pred_dict:
        continue
    safe = model_name.replace(' ', '_')
    save_pred_ensemble(
        pred_dict,
        common_feat_500fg_common.index,
        os.path.join(output_dir, f'{safe}_500fg_test_predicted_age_ensemble.csv'),
    )

for model_name, pred_dict in external_preds_by_model.items():
    if not pred_dict:
        continue
    safe = model_name.replace(' ', '_')
    save_pred_ensemble(
        pred_dict,
        external_index_for_save,
        os.path.join(output_dir, f'{safe}_300bcg_predicted_age_ensemble.csv'),
    )

# Save metrics to CSV
results_df = pd.DataFrame()
for name, model_results in results.items():
    for i, result in enumerate(model_results, 1):
        row = {
            'Model': name,
            'Iteration': i,
            'Train_R2': result['train_R2'],
            'Test_R2': result['test_R2'],
            'External_R2': result['external_R2'],
            'Train_RMSE': result['train_RMSE'],
            'Test_RMSE': result['test_RMSE'],
            'External_RMSE': result['external_RMSE'],
            'Train_MSE': result['train_MSE'],
            'Test_MSE': result['test_MSE'],
            'External_MSE': result['external_MSE']
        }
        results_df = pd.concat([results_df, pd.DataFrame([row])], ignore_index=True)

# Save detailed metrics table
#results_df.to_csv(os.path.join(output_dir, 'model_results_300bcg_external_with_methy_share_feattop235.csv'), index=False)

# Create DataFrame for plotting
plot_data = []
for name in results.keys():
    for result in results[name]:
        plot_data.extend([
            {'Model': name, 'Set': 'Training', 'R2': result['train_R2']},
            {'Model': name, 'Set': 'Test', 'R2': result['test_R2']},
            {'Model': name, 'Set': 'External', 'R2': result['external_R2']}
        ])

plot_df = pd.DataFrame(plot_data)

# Save performance distribution plot
plt.figure(figsize=(15, 8))
sns.boxplot(data=plot_df, x='Model', y='R2', hue='Set', 
            palette=['lightblue', 'lightgreen', 'salmon'], 
            showfliers=False)
sns.stripplot(data=plot_df, x='Model', y='R2', hue='Set',
              palette=['lightblue', 'lightgreen', 'salmon'],
              dodge=True,
              alpha=0.6,
              size=3)
plt.title('R² Distribution Across Training, Validation, and External Validation Sets', fontsize=14, pad=20)
plt.xlabel('Model', fontsize=12)
plt.ylabel('R²', fontsize=12)
plt.xticks(rotation=45)
plt.legend(title='Set', bbox_to_anchor=(1.15, 1), loc='upper left')
plt.tight_layout(rect=[0, 0, 0.9, 1])
#plt.savefig(os.path.join(output_dir, 'tabpfn_model_300bcg_external_with_methy_share_feattop235.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"Results have been saved to {output_dir}")
print("Files saved:")
print(f"- model_results.csv: Detailed results in CSV format")
print(f"- saved_models/: Fitted models per iteration (joblib)")
print(f"- *_500fg_train_predicted_age_ensemble.csv")
print(f"- *_500fg_test_predicted_age_ensemble.csv")
print(f"- *_300bcg_predicted_age_ensemble.csv")
print(f"- model_comparison.png: Visualization of model performance")


# Create the plot with larger figure size
plt.figure(figsize=(15, 8))

# Create boxplot with jitter
sns.boxplot(data=plot_df, x='Model', y='R2', hue='Set', 
            palette=['lightblue', 'lightgreen', 'salmon'], 
            showfliers=False)

# Add jittered points
sns.stripplot(data=plot_df, x='Model', y='R2', hue='Set',
              palette=['lightblue', 'lightgreen', 'salmon'],
              dodge=True,
              alpha=0.6,
              size=3)

# Customize the plot
plt.title('R² Distribution Across Training, Validation, and External Validation Sets', fontsize=14, pad=20)
plt.xlabel('Model', fontsize=12)
plt.ylabel('R²', fontsize=12)
plt.xticks(rotation=45)

# Adjust legend
plt.legend(title='Set', bbox_to_anchor=(1.15, 1), loc='upper left')

# Adjust layout with more space for legend
plt.tight_layout(rect=[0, 0, 0.9, 1])

# Save the plot with high resolution
#plt.savefig('/vol/projects/yzhang/300BCG/output/tab_model/tabpfn_model_500fg_300bcg_with_methy_share_feattop235.png', dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

# Print summary statistics
print("\nSummary Statistics:")
print("-" * 50)
for name in results.keys():  # iterate model names directly
    train_r2 = [r['train_R2'] for r in results[name]]
    test_r2 = [r['test_R2'] for r in results[name]]
    external_r2 = [r['external_R2'] for r in results[name]]
    print(f"\n{name}:")
    print(f"\n{name}:")
    print("Training Set:")
    print(f"Mean R²: {np.mean(train_r2):.4f}")
    print(f"Std R²: {np.std(train_r2):.4f}")
    print("Validation Set:")
    print(f"Mean R²: {np.mean(test_r2):.4f}")
    print(f"Std R²: {np.std(test_r2):.4f}")
    print("External Validation Set:")
    print(f"Mean R²: {np.mean(external_r2):.4f}")
    print("-" * 50)