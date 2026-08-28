library(pROC)
library(caret)
library(dplyr)
library(glmnet)
library(readr)
library(ggplot2)
library(parallel)

#Mvalue <- readRDS("/vol/projects/CIIM/methylation/500FG/1QC/Output/500FG.Mvalue.rds")
#Mvalue_trans <- as.data.frame(t(Mvalue))
#colnames(Mvalue_trans) <- rownames(Mvalue)
#rownames(Mvalue_trans) <- colnames(Mvalue)
#dim(Mvalue_trans)
#head(Mvalue_trans[,1:4])
#saveRDS(Mvalue_trans, "/vol/projects/yzhang/500FG_aging/input/methylation/Mvalue_trans.rds")
Mvalue_trans <- readRDS("/vol/projects/yzhang/500FG_aging/input/methylation/Mvalue_trans.rds")
spearman1 <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction/Mvalue_spearman_preparation/all_iterations_results1.rds") 
spearman2 <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction/Mvalue_spearman_preparation/all_iterations_results2.rds")
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
Mvalue <- Mvalue_trans
#common_samples <- readRDS("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/multi-omics_common_sample_without_microbiome.rds") 
#Mvalue <- Mvalue[common_samples, , drop = FALSE]
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

Mvalue_age <- cbind(Mvalue_filter, Age = basicPhenos_filter_match$Age)
head(Mvalue_age[,1:4])
dim(Mvalue_age)

###save each model iteration to use in the mutli-omics data
# Create a folder to save data for each experiment
dir.create("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/Mvalue_spearman_pvalue_select1e04_prediction_elastic_results", showWarnings = FALSE)

# Custom function to train and evaluate, saving split indices and selected features
train_and_evaluate <- function(data, seed, iteration, pvalue_threshold) {
  # Set RNG kind and seed at the start, exactly as in code 2
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  # Data splitting
  split_indices <- createDataPartition(data$Age, p = 0.7, list = FALSE)
  training <- data[split_indices, ]
  validation <- data[-split_indices, ]
  
 # Get the corresponding iteration results from both spearman results
  spearman1_iter <- spearman1[[iteration]]
  spearman2_iter <- spearman2[[iteration]]
  
  # Combine p-values from both results
  all_pvalues <- c(spearman1_iter$spearman_pvalues, spearman2_iter$spearman_pvalues)
  all_features <- c(spearman1_iter$selected_features, spearman2_iter$selected_features)
  
  # Ensure we're using the same split indices as in the spearman results
  if (!identical(split_indices, spearman1_iter$split_indices)) {
    warning("Split indices mismatch. Using original split indices from spearman results.")
    split_indices <- spearman1_iter$split_indices
    training <- data[split_indices, ]
    validation <- data[-split_indices, ]
  }
  
  # Filter features based on p-value threshold
  selected_features <- all_features[all_pvalues < pvalue_threshold]
  
  # Ensure we have selected features
  if (length(selected_features) == 0) {
    stop(paste("No features selected after filtering by p-value threshold", pvalue_threshold, "in iteration", iteration))
  }
  
  # # Save split indices and selected features
  # saveRDS(list(split_indices = split_indices, selected_features = selected_features),
  #         file = paste0("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/Mvalue_spearman_pvalue_select1e04_prediction_results/iteration_", iteration, "_info.rds"))
  
  # Subset data to selected features
  training_selected <- training[, c(selected_features, "Age"), drop = FALSE]
  validation_selected <- validation[, c(selected_features, "Age"), drop = FALSE]
  
  rm(training)
  rm(validation)
  gc()

  # Model training using glmnet
  #netGrid <- expand.grid(.alpha = 1, .lambda = seq(0, 1, by = 0.01))
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
       res <- train_and_evaluate(Mvalue_age, seed = 123 + i, iteration = i, pvalue_threshold = 1e-04)
      return(res)
    }
)
# Add batch results to results
  results <- c(results, batch_results)
  
  rm(batch_results)
  gc()
}

# Save all results in a single file
saveRDS(results, file = "/vol/projects/yzhang/500FG_aging/output/12_prediction/methylation_prediction_pvalue_select1e04_results_spearman_elasticnet.rds")   
#results <- lapply(1:100, function(i) train_and_evaluate(Mvalue_age, seed = i, iteration = i))

                                 
# results <- mclapply(1:100, function(i) {
#  res <- train_and_evaluate(Mvalue_age, seed = i, iteration = i)
#  gc()  # 释放内存
#  return(res)
#}, mc.cores = 7) 
                                 
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
png("/vol/projects/yzhang/500FG_aging/output/12_prediction/methylation_cross_validation_pvalue_select1e04_spearman_elasticnet.png", 
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

                        