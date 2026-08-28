library(pROC)
library(caret)
library(dplyr)
library(glmnet)
library(readr)
library(ggplot2)
library(parallel)


Mvalue_trans <- readRDS("/vol/projects/yzhang/500FG_aging/input/methylation/Mvalue_trans.rds")
# spearman1 <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction/Mvalue_spearman_preparation/all_iterations_common_results1.rds") 
# spearman2 <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction/Mvalue_spearman_preparation/all_iterations_common_results2.rds")
#Mvalue_trans <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/methy_5000+clock_predicition/Mvalue_clock.rds") #use the filter mvalue and epi_clock cpgs sites 
dim(Mvalue_trans)
head(Mvalue_trans[,1:4])

# Calculate the variance of each variable
#variances <- apply(Mvalue_trans, 2, var)
# Set a threshold to remove low-variance variables
#threshold <- 0.1
#filtered_data <- Mvalue_trans[, variances > threshold]
#dim(filtered_data)
#Mvalue <- filtered_data

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv")

###use 192 common samples
common_samples <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/alllayer_common_samples_nogenptype.rds") 
Mvalue <- Mvalue_trans[common_samples, , drop = FALSE]
dim(Mvalue)
head(Mvalue[,1:4])

# prepare cytokine_filter and EAA_filter_match
idx = intersect(rownames(Mvalue), basicPhenos$ID_500fg) 
length(idx) #305samples
Mvalue_filter = Mvalue[which(rownames(Mvalue) %in% idx),]
head(Mvalue_filter[,1:4])  # 305 854369
dim(Mvalue_filter)

basicPhenos_filter = basicPhenos %>% filter(ID_500fg %in% idx)
head(basicPhenos_filter[,1:4]) #458samples
dim(basicPhenos_filter)

sum(rownames(Mvalue_filter) == basicPhenos_filter$ID_500fg)

basicPhenos_filter_match = basicPhenos_filter[match(rownames(Mvalue_filter), basicPhenos_filter$ID_500fg), ]

sum(basicPhenos_filter_match$ID_500fg == rownames(Mvalue_filter))
head(basicPhenos_filter_match[,1:4]) 
dim(basicPhenos_filter_match) #458samples

# Mvalue_age <- cbind(Mvalue_filter, Age = basicPhenos_filter_match$Age)
 # Add Gender as a predictor (same logic as omic_original_size_with_methy_spear_gender):
# - trim whitespace; treat "" as missing
# - drop samples with missing Gender so M-value and phenotypes stay aligned
basicPhenos_filter_match$Gender <- trimws(as.character(basicPhenos_filter_match$Gender))
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA
keep_gender <- !is.na(basicPhenos_filter_match$Gender)
Mvalue_filter <- Mvalue_filter[keep_gender, , drop = FALSE]
basicPhenos_filter_match <- basicPhenos_filter_match[keep_gender, , drop = FALSE]
stopifnot(all(rownames(Mvalue_filter) == basicPhenos_filter_match$ID_500fg))

#Mvalue_age <- cbind(Mvalue_filter, Gender = basicPhenos_filter_match$Gender, Age = basicPhenos_filter_match$Age)
Mvalue_age <- cbind(Mvalue_filter,Age = basicPhenos_filter_match$Age)
head(Mvalue_age[,1:4])
dim(Mvalue_age)

###save each model iteration to use in the mutli-omics data
# Create a folder to save data for each experiment
dir.create("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/search/elasticnet/Mvalue_spearman_pvalue_common_select_prediction_results_elasticnet", showWarnings = FALSE)

# Custom function to train and evaluate, saving split indices and selected features
train_and_evaluate <- function(data, seed, iteration, pvalue_threshold) {
    
# # Gender cleanup (one missing value as blank string); align with omic gender pipeline
#   if (!"Gender" %in% names(data)) {
#     stop("Gender column not found in input data.")
#   }
#   data$Gender <- trimws(as.character(data$Gender))
#   data$Gender[data$Gender == ""] <- NA

  Mvalue_info <- readRDS(paste0(
    "/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction/spearman_preparation_FDR/FDR_common005/methylation_combined_spearman_preparation_rds/iteration_",
    iteration, "_fdr_info.rds"
  ))

  split_indices1 <- as.integer(unlist(Mvalue_info$split_indices1))
  split_indices2 <- as.integer(unlist(Mvalue_info$split_indices2))

  if (!identical(sort(split_indices1), sort(split_indices2))) {
    warning("Split indices mismatch. Using original split indices from spearman results.")
  }
  split_indices <- split_indices1

  training <- data[split_indices, ]
  validation <- data[-split_indices, ]

  selected_features <- Mvalue_info$selected_features
  if (length(selected_features) > 20000) {
    names(Mvalue_info$fdr) <- selected_features
    selected_features <- names(sort(Mvalue_info$fdr, decreasing = FALSE))[1:10000]
    cat("Mvalue_features > 20000, using FDR top 15000 features.\n")
  }
  selected_features <- unique(selected_features)
    
  # Ensure we have selected features
  if (length(selected_features) == 0) {
    stop(paste("No features selected after filtering by p-value threshold", pvalue_threshold, "in iteration", iteration))
  }
  
  # Save split indices and selected features
  saveRDS(list(
    split_indices = split_indices,
    selected_features = selected_features,
    fdr = Mvalue_info$fdr
  ),
          file = paste0("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/search/elasticnet/Mvalue_spearman_pvalue_common_select_prediction_results_elasticnet/iteration_", iteration, "_info.rds"))
  
  # Subset data to selected features
  training_selected <- training[, c(selected_features, "Age"), drop = FALSE]
  validation_selected <- validation[, c(selected_features, "Age"), drop = FALSE]
  
  # Drop rows with missing Gender after split, then encode Gender as numeric for glmnet
  # training_selected <- training_selected[!is.na(training_selected$Gender), , drop = FALSE]
  # validation_selected <- validation_selected[!is.na(validation_selected$Gender), , drop = FALSE]
  # gender_levels <- sort(unique(na.omit(data$Gender)))
  # training_selected$Gender <- as.numeric(factor(training_selected$Gender, levels = gender_levels))
  # validation_selected$Gender <- as.numeric(factor(validation_selected$Gender, levels = gender_levels))
  
  rm(training)
  rm(validation)
  gc()

  # Model training using glmnet
  # netGrid <- expand.grid(.alpha = 1, .lambda = seq(0, 1, by = 0.01))  
  # Model training using glmnet (Elastic Net), consistent with multi-omics gender pipeline
  netGrid <- expand.grid(.alpha = seq(0.1, 0.9, by = 0.1), .lambda = seq(0, 1, by = 0.01))  
  netctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
  
  netFit <- train(Age ~ ., data = training_selected, method = "glmnet", metric = "RMSE",
                  tuneGrid = netGrid, trControl = netctrl, preProcess = "scale")
  
  # Make predictions
  train_predictions <- predict(netFit, newdata = training_selected)
  val_predictions <- predict(netFit, newdata = validation_selected)
  
  rm(netFit)
  gc()
      
      # Calculate performance metrics for the training set
train_rmse <- sqrt(mean((training_selected$Age - train_predictions)^2))  # RMSE
train_mse <- mean((training_selected$Age - train_predictions)^2)         # MSE
train_r_squared <- cor(training_selected$Age, train_predictions)^2       # R²

# Calculate performance metrics for the validation set
val_rmse <- sqrt(mean((validation_selected$Age - val_predictions)^2))    # RMSE
val_mse <- mean((validation_selected$Age - val_predictions)^2)           # MSE
val_r_squared <- cor(validation_selected$Age, val_predictions)^2         # R²
  
  rm(training_selected, validation_selected, train_predictions, val_predictions)
  gc()
      
      # Return results
return(list(
  train_RMSE = train_rmse,
  train_MSE = train_mse,               # Added MSE
  train_R2 = train_r_squared,
  val_RMSE = val_rmse,
  val_MSE = val_mse,                   # Added MSE
  val_R2 = val_r_squared,
  selected_features = selected_features
))
}
                         
# Run 100 experiments, saving split indices and selected features for each iteration
set.seed(1)

batch_size <- 10  # 每个 batch 运行 10 次
num_batches <- 10  # 总 batch 数
results <- list() 

#for (batch in 1:num_batches) {
#  batch_results <- mclapply(
 #   ((batch - 1) * batch_size + 1):(batch * batch_size), 
 #   function(i) {
 #      res <- train_and_evaluate(Mvalue_age, seed = 123 + i, iteration = i, pvalue_threshold = 1e-04)
 #     return(res)
 #   }, 
 #   mc.cores = 2
 # )
for (batch in 1:num_batches) {
  batch_results <- lapply(
    ((batch - 1) * batch_size + 1):(batch * batch_size), 
    function(i) {
      tryCatch({
        res <- train_and_evaluate(Mvalue_age, seed = 123 + i, iteration = i, pvalue_threshold = 0.05)
        return(res)
      }, error = function(e) {
        cat("Error in iteration", i, ":", e$message, "\n")
        return(NULL)
      })
    }
  )
# Add batch results to results
  results <- c(results, batch_results)
  
  rm(batch_results)
  gc()
}

# Save all results in a single file
saveRDS(results, file = "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/methylation_prediction_pvalue_common_results_spearman_elasticnet.rds")   
#results <- lapply(1:100, function(i) train_and_evaluate(Mvalue_age, seed = i, iteration = i))

                                 
# results <- mclapply(1:100, function(i) {
#  res <- train_and_evaluate(Mvalue_age, seed = i, iteration = i)
#  gc()  # 释放内存
#  return(res)
#}, mc.cores = 7) 
  results <- results[!sapply(results, is.null)]

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
                              
 plot_data_r2 <- data.frame(
  Set = rep(c("Training", "Validation"), each = 100),
  R2 = c(train_r2_values, val_r2_values)
)                                            
# Plot only R²
ppi <- 300
png("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/methylation_cross_validation_pvalue_common_spearman_elasticnet.png", 
    width = 5 * ppi, height = 4 * ppi, res = ppi)

g<-ggplot(plot_data_r2, aes(x = Set, y = R2, fill = Set)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(position = position_jitterdodge(), size = 1, alpha = 0.6, color = "black") +  # Scatter points
  labs(title = "R² Distribution Across Training and Validation in Methylation Sets",
       x = "Methylation",
       y = "R²") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
print(g)
dev.off()                        

                        