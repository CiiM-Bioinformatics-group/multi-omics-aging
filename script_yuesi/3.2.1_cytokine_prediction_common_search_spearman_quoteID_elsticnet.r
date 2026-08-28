library(pROC)
library(caret)
library(dplyr)
library(glmnet)
library(parallel)
library(ggplot2)

load("/vol/projects/yzhang/500FG_aging/output/12_prediction/cytokine_prediction_results_spearman_select_gender.Rdata")

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv")
head(basicPhenos)
#basicPhenos <- basicPhenos[!is.na(basicPhenos$Age), ] None missing value
dim(basicPhenos)

#cytokine <- read.csv("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/pheno_91cytokines_4Raul.csv", sep = ",")
# cytokine <- read.csv("/vol/projects/yzhang/500FG_aging/output/01_cytokine/filtered_stimulation_cytokine_filled.csv") #use mean value filled the missing value
# rownames(cytokine) <- cytokine[,1]
# cytokine <- cytokine[,-1]
# head(cytokine)
# dim(cytokine) #489samples 91 cytokines
# sum(is.na(cytokine))

# cytokine_filter_transposed <- as.data.frame(t(cytokine))
# colnames(cytokine_filter_transposed) <- rownames(cytokine)
# rownames(cytokine_filter_transposed) <- colnames(cytokine)
# head(cytokine_filter_transposed)
# #saveRDS(cytokine_filter_transposed,file="/vol/projects/yzhang/500FG_aging/input/cytokine/cytokine_filter_transposed.rds")
# cytokine<-cytokine_filter_transposed

cytokine <-  read.csv("/vol/projects/yzhang/500FG_aging/input/cytokine/cytokine_rank_normalized.csv",row.names=1)
head(cytokine)
dim(cytokine)
sum(is.na(cytokine))

##fliter common sample in all layer
#common_samples <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/alllayer_common_samples.rds")   
#common_samples <- readRDS("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/multi-omics_common_sample_without_microbiome.rds") 
common_samples <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/alllayer_common_samples_nogenptype.rds")
cytokine <- cytokine[common_samples, , drop = FALSE]
head(cytokine)
dim(cytokine)

# prepare cytokine_filter and pheno_filter_match
idx = intersect(rownames(cytokine), basicPhenos$ID_500fg) 
length(idx) #458samples
cytokine_filter = cytokine[which(rownames(cytokine) %in% idx),]
head(cytokine_filter)  #458samples,1596 metabolite
dim(cytokine_filter)

basicPhenos_filter = basicPhenos %>% filter(ID_500fg %in% idx)
head(basicPhenos_filter) #458samples
dim(basicPhenos_filter)

sum(rownames(cytokine_filter) == basicPhenos_filter$ID_500fg)

basicPhenos_filter_match = basicPhenos_filter[match(rownames(cytokine_filter), basicPhenos_filter$ID_500fg), ]

sum(basicPhenos_filter_match$ID_500fg == rownames(cytokine))
head(basicPhenos_filter_match) 
dim(basicPhenos_filter_match) #458samples
#saveRDS(cytokine_filter, "/vol/projects/yzhang/500FG_aging/input/cytokine/alllayer_common_cytokine_filter_age_high_var.rds")
#saveRDS(basicPhenos_filter_match,"/vol/projects/yzhang/500FG_aging/input/cytokine/alllayer_common_basicPhenos_filter_match_cytokine_age_high_var.rds")

cytokine_age <- cbind(cytokine_filter, Age = basicPhenos_filter_match$Age)
head(cytokine_age)
dim(cytokine_age)

##calculate the spearman correlation resluts first, then search
# Create directory to store results
dir.create("/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction_all_feature/cytokine_spearman_preparation", showWarnings = FALSE)

# Load the all_iterations_common_results1.rds file
all_iterations_common_results1 <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction/Mvalue_spearman_preparation/all_iterations_common_results1.rds")

# Initialize list to store all iterations results
all_results <- list()

# Function to train and evaluate
train_and_evaluate <- function(data, seed, iteration) {
  #set.seed(seed)
   # Get the split indices from all_iterations_common_results1.rds
   iteration_result <- all_iterations_common_results1[[iteration]]
  split_indices <- iteration_result$split_indices[, 1]  # Extract the indices from the matrix
  # Split data
 # split_indices <- createDataPartition(data$Age, p = 0.7, list = FALSE)
    
  training <- data[split_indices, ]
  validation <- data[-split_indices, ]
  print(dim(training))  # Should have ~70% of the original rows
  # Should be 5000 (excluding Age) 
 
  # Get feature names
  feature_names <- setdiff(names(training), "Age")  
  feature_batches <- split(feature_names, ceiling(seq_along(feature_names) / 500))  
    print(length(feature_names))

  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
    
  # Compute Spearman correlation
  spearman_results <- mclapply(feature_batches, function(batch) {
    sapply(batch, function(feature) {
      if (!feature %in% names(training)) return(c(cor = NA, pvalue = NA))
      if (var(training[[feature]], na.rm = TRUE) == 0) return(c(cor = NA, pvalue = NA))
      
      test <- tryCatch(
        cor.test(training[[feature]], training$Age, method = "spearman"),
        error = function(e) return(NULL)
      )

      if (is.null(test)) return(c(cor = NA, pvalue = NA))
      return(c(cor = test$estimate, pvalue = test$p.value))
    }, simplify = FALSE)
  }, mc.cores = 8) 
      
       if (length(spearman_results) == 0 || is.null(spearman_results[[1]])) {
    stop("Error: spearman_results is empty or NULL. Check execution.")
  }
  print(length(spearman_results))
  # Extract p-values
  spearman_pvalues <- unlist(lapply(spearman_results, function(batch) {
    sapply(batch, function(x) if (!is.null(x)) x["pvalue"] else NA)
  }))
   print(summary(spearman_pvalues))  

  # Get corresponding feature names
  #feature_names <- names(spearman_pvalues)
 # names(spearman_pvalues) <- feature_names
   feature_names <- unlist(lapply(spearman_results, names))
  print(length(spearman_pvalues))  # p-value 总数
  print(head(spearman_pvalues))    # 预览前几个 p-values
  print(names(spearman_pvalues)[1:10])       
 
 
if (length(feature_names) != length(spearman_pvalues)) {
    stop("Error: Mismatch between feature names and spearman_pvalues length.")
  }
   names(spearman_pvalues) <- feature_names
  
    # Filter features with p-value < 0.05
   # selected_features <- names(spearman_pvalues[spearman_pvalues < 0.1])        
   selected_features <- names(spearman_pvalues) 
           
  # Store results in list
  results <- list(
    split_indices = split_indices,
    selected_features = selected_features,
    spearman_pvalues = spearman_pvalues[selected_features] # 只存筛选出的 p-value
  )

  rm(training, validation, feature_batches, spearman_results, spearman_pvalues, feature_names)
  gc()  

  return(results)
}

           
# Run 100 iterations
#for (i in 1:2) {
#  cat("Running iteration", i, "\n")
#  train_and_evaluate(Mvalue_age, seed = 123 + i, iteration = i)
#}
#iteration_results <- mclapply(1:100, function(i) {
iteration_results <- lapply(1:100, function(i) {
  train_and_evaluate(cytokine_age, seed = 123 + i, iteration = i)
})

# 合并结果
for (i in 1:100) {
  all_results[[i]] <- iteration_results[[i]]
}

           # Save all results in a single RDS file
if (length(all_results) > 0) {
  saveRDS(all_results, "/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction_all_feature/cytokine_spearman_preparation/cytokine_iterations_results_quoteID.rds")
} else {
  cat("Error: all_results is empty, nothing to save.\n")
}


str(all_results)

# 先清洗 Gender
basicPhenos_filter_match$Gender <- trimws(as.character(basicPhenos_filter_match$Gender))
basicPhenos_filter_match$Gender[basicPhenos_filter_match$Gender == ""] <- NA

# 用 ID_500fg 建立映射并按 cytokine_age 的 rownames 对齐
id2gender <- setNames(basicPhenos_filter_match$Gender, basicPhenos_filter_match$ID_500fg)
cytokine_age$Gender <- id2gender[rownames(cytokine_age)]

head(cytokine_age)
dim(cytokine_age)

###save each model iteration to use in the mutli-omics data
# Create a folder to save data for each experiment
dir.create("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/search/elasticnet/cytokine_spearman_pvalue_common_select_prediction_results_quoteID_elasticnet", showWarnings = FALSE)

# Load the all_iterations_common_results1.rds file
all_iterations_common_results1 <- readRDS("/vol/projects/yzhang/500FG_aging/input/prediction/spearman_prediction/Mvalue_spearman_preparation/all_iterations_common_results1.rds")

# Custom function to train and evaluate, saving split indices and selected features
train_and_evaluate <- function(data, seed, iteration, pvalue_threshold) {
  # Set RNG kind and seed at the start, exactly as in code 2
  #RNGkind("L'Ecuyer-CMRG")
  #set.seed(seed)
  # Data splitting
  #split_indices <- createDataPartition(data$Age, p = 0.7, list = FALSE)

    # Gender cleanup (one missing value as blank string)
  # if (!"Gender" %in% names(data)) {
  #   stop("Gender column not found in input data.")
  # }
  # data$Gender <- trimws(as.character(data$Gender))
  # data$Gender[data$Gender == ""] <- NA
    
    # Get the split indices from all_iterations_common_results1.rds
   iteration_result <- all_iterations_common_results1[[iteration]]
  split_indices <- iteration_result$split_indices[, 1]  # Extract the indices from the matrix
    
  training <- data[split_indices, ]
  validation <- data[-split_indices, ]
  
 # Get the corresponding iteration results from both spearman results
  iteration_result <- all_results[[iteration]]

# Get p-values and features from the iteration result
all_pvalues <- iteration_result$spearman_pvalues
all_features <- iteration_result$selected_features
  
  # Ensure we're using the same split indices as in the spearman results
  if (!identical(split_indices, iteration_result$split_indices)) {
  warning("Split indices mismatch. Using original split indices from results.")
  split_indices <- iteration_result$split_indices
  training <- data[split_indices, ]
  validation <- data[-split_indices, ]
}

# Filter features based on p-value threshold
selected_features <- all_features[all_pvalues < pvalue_threshold]
    
 # Force Gender as a fixed feature
# selected_features <- unique(c(selected_features, "Gender"))
    
  # Ensure we have selected features
  if (length(selected_features) == 0) {
    stop(paste("No features selected after filtering by p-value threshold", pvalue_threshold, "in iteration", iteration))
  }
  
  # Save split indices and selected features
  saveRDS(list(split_indices = split_indices, selected_features = selected_features),
          file = paste0("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/search/elasticnet/cytokine_spearman_pvalue_common_select_prediction_results_quoteID_elasticnet/iteration_", iteration, "_info.rds"))
  
  # Subset data to selected features
  training_selected <- training[, c(selected_features, "Age"), drop = FALSE]
  validation_selected <- validation[, c(selected_features, "Age"), drop = FALSE]
    
  # Drop rows with missing Gender after split, then encode Gender as numeric
  # training_selected <- training_selected[!is.na(training_selected$Gender), , drop = FALSE]
  # validation_selected <- validation_selected[!is.na(validation_selected$Gender), , drop = FALSE]
  # gender_levels <- sort(unique(na.omit(data$Gender)))
  # training_selected$Gender <- as.numeric(factor(training_selected$Gender, levels = gender_levels))
  # validation_selected$Gender <- as.numeric(factor(validation_selected$Gender, levels = gender_levels))
   
  rm(training)
  rm(validation)
  gc()

  # Model training using glmnet
  #netGrid <- expand.grid(.alpha = 1, .lambda = seq(0, 1, by = 0.01)) 
  # Model training using glmnet (Elastic Net)
  netGrid <- expand.grid(.alpha = seq(0.1, 0.9, by = 0.1), .lambda = seq(0, 1, by = 0.01))  
  netctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
  
  netFit <- train(Age ~ ., data = training_selected, method = "glmnet", metric = "RMSE",
                  tuneGrid = netGrid, trControl = netctrl, preProcess = "scale")
  
  # Make predictions
  train_predictions <- predict(netFit, newdata = training_selected)
  val_predictions <- predict(netFit, newdata = validation_selected)
  
  rm(netFit)
  gc()
    
 # Save predicted age for each sample
  pred_df <- rbind(
    data.frame(
      iteration = iteration,
      set = "train",
      sample_id = rownames(training_selected),
      Age_true = training_selected$Age,
      Age_pred = as.numeric(train_predictions)
    ),
    data.frame(
      iteration = iteration,
      set = "validation",
      sample_id = rownames(validation_selected),
      Age_true = validation_selected$Age,
      Age_pred = as.numeric(val_predictions)
    )
  )
     
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
   selected_features = selected_features,
  predictions = pred_df
))
}
                      

set.seed(1)

batch_size <- 10  # 每个 batch 运行 10 次
num_batches <- 10  # 总 batch 数
results <- list() 

for (batch in 1:num_batches) {
    batch_results <- lapply(
    ((batch - 1) * batch_size + 1):(batch * batch_size), 
    function(i) {
       res <- train_and_evaluate(cytokine_age, seed = 123 + i, iteration = i, pvalue_threshold = 1)
      return(res)
    }
)
# Add batch results to results
  results <- c(results, batch_results)
  
  rm(batch_results)
  gc()
}

# Save all results in a single file
saveRDS(results, file = "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/cytokine_prediction_pvalue_common_results_spearman_quoteID_elasticnet.rds")   
#results <- lapply(1:100, function(i) train_and_evaluate(Mvalue_age, seed = i, iteration = i))

# Save all predicted ages from 100 iterations in one file
all_predictions <- do.call(rbind, lapply(results, function(x) x$predictions))
write.csv(all_predictions, file = "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/cytokine_prediction_pvalue_common_select_age_predictions_quoteID_elasticnet.csv", row.names = FALSE)


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
# ppi <- 300
# png("/vol/projects/yzhang/500FG_aging/output/12_prediction/cytokine_predicition_pvalue_common_select_spearman.png", 
#     width = 5 * ppi, height = 4 * ppi, res = ppi)

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
# dev.off()    
# print(g)

save.image("/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet/cytokine_prediction_results_spearman_select_elasticnet.Rdata")

info <- readRDS("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/Mvalue_spearman_prediction_results/iteration_1_info.rds")
str(info)

info <- readRDS("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/Mvalue_spearman_prediction_results/iteration_2_info.rds")
str(info)

info <- readRDS("/vol/projects/yzhang/500FG_aging/code/prediction/prediction_vari_spearman/Mvalue_spearman_prediction_results/iteration_3_info.rds")
str(info)
