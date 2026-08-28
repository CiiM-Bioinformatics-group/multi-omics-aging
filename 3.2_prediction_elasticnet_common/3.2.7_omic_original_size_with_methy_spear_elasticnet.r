library(pROC)
library(caret)
library(dplyr)
library(glmnet)
library(readr)
library(jsonlite)

Mvalue <-  readRDS("/vol/projects/yzhang/500FG_aging/input/methylation/Mvalue_trans.rds")
head(Mvalue[,1:4])
dim(Mvalue)
#genotype <- readRDS("/vol/projects/yzhang/500FG_aging/input/genotype/genotype_trans.rds")
cytokine <-  read.csv("/vol/projects/yzhang/500FG_aging/input/cytokine/cytokine_rank_normalized.csv",row.names=1)
olink  <-  read.csv("/vol/projects/yzhang/500FG_aging/input/olink/olink_rank_normalized.csv",row.names=1)
metabolite  <- read.csv("/vol/projects/yzhang/500FG_aging/input/metabolite/metabolite_rank_normalized.csv",row.names=1) 
#cellcounts  <- load_and_process("/vol/projects/yzhang/500FG_aging/input/variance/cumulative_variance_input/cellCounts_variance_input.csv")  
cellcounts = read.csv("/vol/projects/yzhang/500FG_aging/input/cellCounts/cellcounts_name_replace.csv", row.names=1) ##less ##less
#microbiome <- read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecux flar_phenotype/500FG_microbiome_pathways.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)
microbiome <- read.csv("/vol/projects/yzhang/500FG_aging/input/microbiome/microbiome_rank_normalized.csv",row.names=1)
basicPhenos <- read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv",row.names=1)
head(basicPhenos)
dim(basicPhenos)

# Step 1: Identify the common samples across all datasets and basicPhenos$ID_500fg
common_samples <- Reduce(intersect, list(
  rownames(cytokine), 
  rownames(olink), 
  rownames(metabolite), 
  rownames(cellcounts), 
  rownames(microbiome),
 rownames(Mvalue), 
    basicPhenos$ID_500fg
))

# Step 2: Subset all datasets to keep only common samples
common_samples <- sort(common_samples)  # Ensure consistent order

cytokine_common <- cytokine[common_samples, , drop = FALSE]
olink_common <- olink[common_samples, , drop = FALSE]
metabolite_common <- metabolite[common_samples, , drop = FALSE]
cellcounts_common <- cellcounts[common_samples, , drop = FALSE]
microbiome_common <- microbiome[common_samples, , drop = FALSE]
Mvalue_common <- Mvalue[common_samples, , drop = FALSE]
basicPhenos_common <- basicPhenos[basicPhenos$ID_500fg %in% common_samples, , drop = FALSE]

# Add Gender as a predictor (minimal change):
# - drop the (one) missing-Gender sample to keep all layers aligned
basicPhenos_common$Gender <- trimws(basicPhenos_common$Gender)
basicPhenos_common$Gender[basicPhenos_common$Gender == ""] <- NA
non_missing_gender_ids <- basicPhenos_common$ID_500fg[!is.na(basicPhenos_common$Gender)]
common_samples <- intersect(common_samples, non_missing_gender_ids)
common_samples <- sort(common_samples)

# Re-subset again after Gender filtering to keep everything consistent
cytokine_common <- cytokine[common_samples, , drop = FALSE]
olink_common <- olink[common_samples, , drop = FALSE]
metabolite_common <- metabolite[common_samples, , drop = FALSE]
cellcounts_common <- cellcounts[common_samples, , drop = FALSE]
microbiome_common <- microbiome[common_samples, , drop = FALSE]
Mvalue_common <- Mvalue[common_samples, , drop = FALSE]
basicPhenos_common <- basicPhenos[basicPhenos$ID_500fg %in% common_samples, , drop = FALSE]
basicPhenos_common$Gender <- trimws(basicPhenos_common$Gender)
basicPhenos_common$Gender[basicPhenos_common$Gender == ""] <- NA
basicPhenos_common$Gender <- factor(basicPhenos_common$Gender)


#saveRDS(common_samples, "/vol/projects/yzhang/500FG_aging/input/prediction/alllayer_common_samples.rds")     
combined_table <- cbind(
  cytokine_common, olink_common, metabolite_common, 
   cellcounts_common, microbiome_common, Mvalue_common
)
# Sort basicPhenos_common by common_samples
basicPhenos_common <- basicPhenos_common[order(match(basicPhenos_common$ID_500fg, common_samples)), ]
# Step 6: Create a separate table for common samples
common_samples_table <- data.frame(Sample_ID = common_samples)

# Output the combined table and the table of common samples
rownames(combined_table) <-common_samples
dim(combined_table)  # Check dimensions of the combined table #291samples 
head(combined_table[,1:4])  # Preview the combined table
head(basicPhenos_common[,1:4])


#combined_age <- cbind(combined_table, Age = basicPhenos_common$Age)
# combined_age <- cbind(combined_table, Gender = basicPhenos_common$Gender, Age = basicPhenos_common$Age)
combined_age <- cbind(combined_table, Age = basicPhenos_common$Age)
head(combined_age[,1:4])
dim(combined_age)
cat("na:", sum(is.na(combined_age)), "\n")
# Calculate the variance of each variable
#variances <- apply(combined_age, 2, var)
# Set a threshold to remove low-variance variables
#threshold <- 0.01
#filtered_data <- combined_age[, variances > threshold]
#dim(filtered_data)
#combined_age <- filtered_data

##100 times
# Function to build multi-omics model using stored results
multi_omics_model_train <- function(training_data, validation_data, iteration) {
  # Model training using glmnet
  #netGrid <- expand.grid(.alpha = 1, .lambda = seq(0, 1, by = 0.01))  
  netGrid <- expand.grid(    # elastic net
  .alpha  = seq(0.1, 0.9, by = 0.1),          
  .lambda = 10^seq(-4, 1, length.out = 100)  
)
  netctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
  
  netFit <- train(Age ~ ., data = training_data, method = "glmnet", metric = "RMSE",
                  tuneGrid = netGrid, trControl = netctrl, preProcess = "scale")
  
  # Make predictions
  train_predictions <- predict(netFit, newdata = training_data)
  val_predictions <- predict(netFit, newdata = validation_data)
  
  
   pred_df <- rbind(
    data.frame(
      iteration = iteration,
      set = "train",
      sample_id = rownames(training_data),
      Age_true = training_data$Age,
      Age_pred = as.numeric(train_predictions)
    ),
    data.frame(
      iteration = iteration,
      set = "validation",
      sample_id = rownames(validation_data),
      Age_true = validation_data$Age,
      Age_pred = as.numeric(val_predictions)
    )
  )
  # Calculate performance metrics
  train_rmse <- sqrt(mean((training_data$Age - train_predictions)^2))
  train_mse <- mean((training_data$Age - train_predictions)^2)
  train_r_squared <- cor(training_data$Age, train_predictions)^2
  
  val_rmse <- sqrt(mean((validation_data$Age - val_predictions)^2))
  val_mse <- mean((validation_data$Age - val_predictions)^2)
  val_r_squared <- cor(validation_data$Age, val_predictions)^2
  
   # Extract coefficients of the best model
  model_coefficients <- as.matrix(coef(netFit$finalModel, s = netFit$bestTune$.lambda))
    
   # Get selected features (non-zero coefficients)
  selected_features <- rownames(model_coefficients)[model_coefficients != 0]
  selected_features <- selected_features[selected_features != "(Intercept)"]
    
  # Return metrics
  list(
    train_RMSE = train_rmse,
    train_MSE = train_mse,
    train_R2 = train_r_squared,
    val_RMSE = val_rmse,
    val_MSE = val_mse,
    val_R2 = val_r_squared, 
    best_model = netFit,
    coefficients = model_coefficients,
    selected_features = selected_features,
    predictions = pred_df
  )
}
cat("Starting multi-omics modeling with 100 iterations...\n")

# Track successful and failed iterations
successful_iterations <- 0
failed_iterations <- 0

# Run multi-omics modeling using stored results
results <- lapply(1:100, function(i) {
  # Load split indices and selected features from each omics model
    tryCatch({
    Mvalue_info <- readRDS(paste0("/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction/spearman_preparation_FDR/FDR_common005/methylation_combined_spearman_preparation_rds/iteration_", i, "_fdr_info.rds"))
    
    # Assume Mvalue_info contains split_indices1 and split_indices2
   split_indices1 <- as.integer(unlist(Mvalue_info$split_indices1))
   split_indices2 <- as.integer(unlist(Mvalue_info$split_indices2))

# Check whether they are identical
if (identical(sort(split_indices1), sort(split_indices2))) {
  cat("split_indices1 and split_indices2 are identical.\n")
} else {
  cat("split_indices1 and split_indices2 are NOT identical!\n")
  cat("split_indices1 (head):", head(split_indices1), "\n")
  cat("split_indices2 (head):", head(split_indices2), "\n")
}

# Use split_indices1
split_indices <- split_indices1  
        
    # Use selected_features only for Mvalue; keep all original features for other layers
    cytokine_features    <- colnames(cytokine)      
    olink_features       <- colnames(olink)
    metabolite_features  <- colnames(metabolite)
    cellcounts_features  <- colnames(cellcounts)
    microbiome_features  <- colnames(microbiome)
    Mvalue_features      <- Mvalue_info$selected_features

  # if Mvalue_features more than 20000，choose FDR top 15000 
if (length(Mvalue_info$selected_features) > 20000) {
  names(Mvalue_info$fdr) <- Mvalue_info$selected_features
  top_features <- names(sort(Mvalue_info$fdr, decreasing = FALSE))[1:10000]
  Mvalue_features <- top_features
  cat("Mvalue_features > 20000, using FDR top 10000 features.\n")
} else {
  Mvalue_features <- Mvalue_info$selected_features
}
        
  # Combine all selected features
  combined_features <- unique(c(cytokine_features, olink_features, metabolite_features,
                             cellcounts_features, microbiome_features,
                              Mvalue_features))
  
  # Ensure Age is not included in combined_features
  if ("Age" %in% combined_features) {
    combined_features <- setdiff(combined_features, "Age")
  }

  # Validate that features exist in combined_age
    available_features <- intersect(combined_features, colnames(combined_age))
    if (length(available_features) == 0) {
      cat("No valid features found in iteration", i, "\n")
      failed_iterations <<- failed_iterations + 1
      return(NULL)
    }
    
    cat("Iteration", i, ":\n")
    cat("Total features:", length(combined_features), "\n")
    cat("Available features:", length(available_features), "\n")
    
        
     ##count features by layer   
    features_count <- c(
    cytokine    = length(intersect(cytokine_features, available_features)),
    olink       = length(intersect(olink_features, available_features)),
    metabolite  = length(intersect(metabolite_features, available_features)),
    cellcounts  = length(intersect(cellcounts_features, available_features)),
    microbiome  = length(intersect(microbiome_features, available_features)),
    Mvalue      = length(intersect(Mvalue_features, available_features))
    )
     # Use available features
    training <- combined_age[split_indices, ]
    validation <- combined_age[-split_indices, ]
    
     # Use available_features (not combined_features)
    training_selected <- training[, c(available_features, "Age"), drop = FALSE]
    validation_selected <- validation[, c(available_features, "Age"), drop = FALSE]
    # # available_features_plus <- unique(c(available_features, "Gender"))
    # training_selected <- training[, c(available_features_plus, "Age"), drop = FALSE]
    # validation_selected <- validation[, c(available_features_plus, "Age"), drop = FALSE]
  
        
    # Train model
    result <- multi_omics_model_train(training_selected, validation_selected, iteration = i)
    result$features_count <- features_count   
    successful_iterations <<- successful_iterations + 1
    return(result)
    
  }, error = function(e) {
    cat("Error in iteration", i, ":", e$message, "\n")
    failed_iterations <<- failed_iterations + 1
    return(NULL)
  })
  rm(list = setdiff(ls(), c("results", "i", "successful_iterations", "failed_iterations"))) # 保留主变量
  gc()
})

# Summary
cat("\nModeling completed:\n")
cat("Total iterations:", 100, "\n")
cat("Successful iterations:", successful_iterations, "\n")
cat("Failed iterations:", failed_iterations, "\n")
cat("Mean:", mean(rowSums(do.call(rbind, lapply(results, function(x) x$features_count)))), 
    "Max:", max(rowSums(do.call(rbind, lapply(results, function(x) x$features_count)))), 
    "Min:", min(rowSums(do.call(rbind, lapply(results, function(x) x$features_count)))), "\n")
                                              
# Remove NULL results
results <- results[!sapply(results, is.null)]
saveRDS(results, "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/multi_prediction_results_original_methy_separ_elasticnet.rds")  

# Save all predicted ages from 100 iterations in one file (same format as single-omics scripts)
all_predictions <- do.call(rbind, lapply(results, function(x) x$predictions))
write.csv(
  all_predictions,
  file = "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/multi_prediction_age_original_methy_separ_elasticnet.csv",
  row.names = FALSE,
  quote = TRUE
)
                                         
# Extract metrics for training and validation sets
train_rmse_values <- sapply(results, function(x) x$train_RMSE)
train_mse_values <- sapply(results, function(x) x$train_MSE)
train_r2_values <- sapply(results, function(x) x$train_R2)

val_rmse_values <- sapply(results, function(x) x$val_RMSE)
val_mse_values <- sapply(results, function(x) x$val_MSE)
val_r2_values <- sapply(results, function(x) x$val_R2)
                        
# Extract number of selected features for each iteration
num_selected_features <- sapply(results, function(x) length(x$selected_features))
                                
# Calculate mean and standard deviation for each metric
metrics_summary <- data.frame(
  Metric = c("Train RMSE", "Validation RMSE", 
             "Train MSE", "Validation MSE", 
             "Train R²", "Validation R²",
             "Number of Selected Features"),
  Mean = c(mean(train_rmse_values), mean(val_rmse_values),
           mean(train_mse_values), mean(val_mse_values),
           mean(train_r2_values), mean(val_r2_values),
           mean(num_selected_features, na.rm = TRUE)), 
  SD = c(sd(train_rmse_values), sd(val_rmse_values),
         sd(train_mse_values), sd(val_mse_values),
         sd(train_r2_values), sd(val_r2_values),
         sd(num_selected_features, na.rm = TRUE)) 
)
                              
# Print summary
print(metrics_summary)
                                                        
# Combine R² metrics into a single data frame for plotting
plot_data_r2 <- data.frame(
  Set = rep(c("Training", "Validation"), each = 100),
  R2 = c(train_r2_values, val_r2_values)
)

# Plot only R²
ppi <- 300
png("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/multi_cross_validation_original_methy_separ_elasticnet.png",  width = 5 * ppi, height = 4 * ppi, res = ppi)

g<-ggplot(plot_data_r2, aes(x = Set, y = R2, fill = Set)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(position = position_jitterdodge(), size = 1, alpha = 0.6, color = "black") +  # Scatter points
  labs(title = "R² Distribution Across Training and Validation in Alldata Sets",
       x = "Alldata",
       y = "R²") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

print(g)
dev.off()
print(g)
                                               
# Step 2: Subset all datasets to keep only common samples
common_samples <- sort(common_samples)  # Ensure consistent order

cytokine_common <- cytokine[common_samples, , drop = FALSE]
olink_common <- olink[common_samples, , drop = FALSE]
metabolite_common <- metabolite[common_samples, , drop = FALSE]
cellcounts_common <- cellcounts[common_samples, , drop = FALSE]
microbiome_common <- microbiome[common_samples, , drop = FALSE]
Mvalue_common <- Mvalue[common_samples, , drop = FALSE]
basicPhenos_common <- basicPhenos[basicPhenos$ID_500fg %in% common_samples, , drop = FALSE]

# Step 1: Define different Layer column names (assuming each dataset is structured correctly)
layer_columns <- list(
  Cytokine = colnames(cytokine_common),
  Protein = colnames(olink_common),
  Metabolite = colnames(metabolite_common),
  CellCounts = colnames(cellcounts_common),
  Microbiome = colnames(microbiome_common),
  Methylation = colnames(Mvalue_common)
)

# Step 1: Extract feature rankings from each model
feature_ranks_list <- list()

for (i in seq_along(results)) {
  model_result <- results[[i]]
  model <- model_result$best_model
  
  # Extract coefficients
  coefs <- coef(model$finalModel, s = model$bestTune$lambda)
  coefs_df <- as.data.frame(as.matrix(coefs))
  coefs_df$Feature <- rownames(coefs_df)
  colnames(coefs_df)[1] <- "Coefficient"
  
  # Remove zero coefficients and intercept
  coefs_df <- coefs_df %>%
    filter(Coefficient != 0 & Feature != "(Intercept)") %>%
    mutate(AbsValue = abs(Coefficient)) %>%
    arrange(desc(AbsValue)) %>%
    mutate(Rank = row_number()) # Compute ranking
  
  # Compute normalized rank
  N <- nrow(coefs_df)
  coefs_df <- coefs_df %>%
    mutate(Normalized_Rank = Rank / N)
  
  # Store in list
  feature_ranks_list[[i]] <- coefs_df %>%
    select(Feature, Normalized_Rank)
}

# Step 2: Compute Rank Product
feature_ranks <- list()

for (i in seq_along(feature_ranks_list)) {
  ranks <- feature_ranks_list[[i]]
  for (j in seq_along(ranks$Feature)) {
    feature <- ranks$Feature[j]
    if (feature %in% names(feature_ranks)) {
      feature_ranks[[feature]] <- c(feature_ranks[[feature]], ranks$Normalized_Rank[j])
    } else {
      feature_ranks[[feature]] <- ranks$Normalized_Rank[j]
    }
  }
}

# Compute Rank Product (geometric mean of ranks)
rank_product <- sapply(feature_ranks, function(ranks) {
  prod(ranks)^(1 / length(ranks))
})

rank_product_df <- data.frame(
  Feature = names(rank_product),
  RankProduct = rank_product
) %>%
  arrange(RankProduct)

# Step 3: Permutation Test for Significance
n_permutations <- 5000
random_rp_values <- matrix(NA, nrow = n_permutations, ncol = nrow(rank_product_df))

for (i in 1:n_permutations) {
  # For each permutation
  random_ranks <- lapply(feature_ranks_list, function(df) {
    # Permute the expression values
    permuted_ranks <- sample(df$Normalized_Rank)
    return(permuted_ranks)
  })
  
  # Calculate RP values for permuted data
  for (j in 1:nrow(rank_product_df)) {
    ranks <- sapply(random_ranks, function(r) r[j])
    random_rp_values[i, j] <- prod(ranks)^(1/length(ranks))
  }
}

# Compute p-values with continuity correction
p_values <- sapply(1:nrow(rank_product_df), function(j) {
  observed_rp <- rank_product_df$RankProduct[j]
  # Add continuity correction
  (sum(random_rp_values[, j] <= observed_rp) + 0.5) / (n_permutations + 1)
})

rank_product_df$PValue <- p_values

# Compute False Discovery Rate (FDR)
rank_product_df <- rank_product_df %>%
  arrange(PValue) %>%
 mutate(
    FDR = p.adjust(PValue, method = "fdr"),
    # Add -log10 transformation for visualization
    neg_log10_pvalue = -log10(PValue)
  )

# Step 4: Select significant features (FDR < 0.05)
significant_features <- rank_product_df %>%
  filter(FDR < 0.05)

# Output results
print(significant_features)

# Save results
saveRDS(rank_product_df, "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/rank_product_df_elasticnet.rds")


layer_counts <- rank_product_df %>%
  mutate(
    # Identify layer based on feature name and layer_columns
    Layer = case_when(
      Feature %in% layer_columns$Cytokine ~ "Cytokine",
      Feature %in% layer_columns$Olink ~ "Olink",
      Feature %in% layer_columns$Metabolite ~ "Metabolite",
      Feature %in% layer_columns$CellCounts ~ "CellCounts",
      Feature %in% layer_columns$Microbiome ~ "Microbiome",
      Feature %in% layer_columns$Mvalue ~ "Mvalue",
      TRUE ~ "Unknown"  # For any features that don't match
    )
  ) %>%
  group_by(Layer) %>%
  summarise(
    Count = n()
  ) %>%
  arrange(desc(Count))

# Print debugging information
print("\nAll features in rank_product_df:")
print(rank_product_df$Feature)

print("\nAll metabolite features in layer_columns:")
print(layer_columns$Metabolite)

# Check for unmatched features
unmatched_features <- rank_product_df %>%
  filter(!Feature %in% unlist(layer_columns)) %>%
  pull(Feature)

print("\nUnmatched features:")
print(unmatched_features)

# Print layer statistics
print("\nFeatures by layer:")
print(layer_counts)

plot_data <- layer_counts %>%
  mutate(Layer = case_when(
    Layer == "Mvalue" ~ "Methylation",
    Layer == "Olink" ~ "Protein",
    TRUE ~ Layer
  ))
# Create a bar plot of layer distribution
p_bar <- ggplot(plot_data, aes(x = reorder(Layer, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Distribution of Features by Layer",
    x = "Layer",
    y = "Number of Features"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )

# Display the bar plot
print(p_bar)

# Save the pie chart
# ggsave("/vol/projects/yzhang/500FG_aging/output/12_prediction/spearman/original_size_with_methy_spear/multi_rank_top10_original_size_with_methy_spear_gender.png", p_bar, width = 8, height = 8, dpi = 300)
                                
# Step 2: Extract the top 100 important features from each model
top_features_list <- lapply(results, function(model_result) {
  model <- model_result$best_model
  coefs <- coef(model$finalModel, s = model$bestTune$lambda)

  # Convert coefficients to a dataframe
  coefs_df <- as.data.frame(as.matrix(coefs))
  coefs_df$Feature <- rownames(coefs_df)
  colnames(coefs_df)[1] <- "Coefficient"

  # Remove intercept and zero coefficients
  coefs_df <- coefs_df %>%
    filter(Coefficient != 0 & Feature != "(Intercept)")

  # Sort by absolute coefficient values
  coefs_df <- coefs_df %>%
    mutate(AbsValue = abs(Coefficient)) %>%
    arrange(desc(AbsValue))

  # Select top 100 features
  top_100 <- head(coefs_df$Feature, 100)
  return(top_100)
})

# Step 3: Count the frequency of all selected features
all_top_features <- unlist(top_features_list)
feature_frequency <- sort(table(all_top_features), decreasing = TRUE)

# Convert to DataFrame
feature_freq_df <- as.data.frame(feature_frequency)
colnames(feature_freq_df) <- c("Feature", "Frequency")

# Step 4: Assign Layer Information to Each Feature
feature_freq_df$Layer <- sapply(feature_freq_df$Feature, function(f) {
  matched_layers <- names(layer_columns)[sapply(layer_columns, function(layer) f %in% layer)]
  if (length(matched_layers) > 0) {
    return(matched_layers[1])  # Assign the first matched layer
  } else {      
    print(paste("Unknown feature:", f)) 
    return("Unknown")  # Assign "Unknown" if no layer matches
  }
})
                                         
feature_freq_df <- feature_freq_df %>% filter(Layer != "Unknown")
write.csv(feature_freq_df,
          file = "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/feature_freq_df_elasticnet.csv",
          row.names = FALSE)
print(unique(feature_freq_df$Layer))
                                         
# Step 5: Count feature contribution per Layer
layer_distribution <- feature_freq_df %>%
  group_by(Layer) %>%
  summarise(
    Total_Features = n(),
    Total_Frequency = sum(Frequency)/100
  ) %>%
  arrange(desc(Total_Features))

 png("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/rank_product_original_size_with_methy_spear_elasticnet.png", width = 600, height = 600)
# Step 6: Visualize Layer Contribution
ggplot(layer_distribution, aes(x = reorder(Layer, -Total_Features), y = Total_Features)) +
  geom_bar(stat = "identity", fill = "pink") +
  labs(title = "Feature Contribution by Layer",
       x = "Layer",
       y = "Number of Features") +
  theme_minimal() +
  coord_flip()

# Step 7: Print the layer distribution summary
print(layer_distribution)
dev.off() 
                                         
# Extract top 10 features from each model's best coefficients
top_features_list <- lapply(results, function(model_result) {
  # Extract best model from each run
  model <- model_result$best_model
  coefs <- coef(model$finalModel, s = model$bestTune$lambda)
  
  coefs_df <- as.data.frame(as.matrix(coefs))
  coefs_df$Feature <- rownames(coefs_df)
  colnames(coefs_df)[1] <- "Coefficient"

  # Remove intercept and zero coefficients
  coefs_df <- coefs_df[coefs_df$Coefficient != 0 & coefs_df$Feature != "(Intercept)", ]

  # Sort by absolute coefficient values and select top 10 features
  coefs_df$AbsValue <- abs(coefs_df$Coefficient)
  top_10 <- coefs_df[order(-coefs_df$AbsValue), ][1:15, "Feature"]

  return(top_10)
})

# Count the frequency of each feature appearing in the top 10 across 100 runs
all_top_features <- unlist(top_features_list)
feature_frequency <- sort(table(all_top_features), decreasing = TRUE)

# Print the 10 most frequently occurring features
print(head(feature_frequency, 15))

                                         png("/vol/projects/yzhang/500FG_aging/output/12_prediction/spearman/elasticnet/multi_frequency_top15_original_size_with_methy_spear_elasticnet.png", width = 800, height = 600)
# Visualize the 10 most frequent features
barplot(head(feature_frequency, 10), las = 2, col = "lightblue",
        main = "Top 15 Most Frequent Features",
        xlab = "Features", ylab = "Frequency")
 dev.off()                                        