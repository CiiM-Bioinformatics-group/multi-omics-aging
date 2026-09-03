library(pROC)
library(caret)
library(dplyr)
library(glmnet)

# cellCounts = read.csv("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_cellCounts_cellPerc.csv")
# rownames(cellCounts) <- cellCounts[,1]
# cellCounts <- cellCounts[,-1]
# head(cellCounts)
# dim(cellCounts)
# any(is.na(cellCounts))

cellCounts = read.table("/vol/projects/CIIM/cohorts_old/500FG/Molecular_phenotype/500FG_inverse_rank_normalized_cellcounts.txt", header = TRUE, stringsAsFactors = FALSE)
head(cellCounts)
dim(cellCounts)
any(is.na(cellCounts))

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv")
# prepare cytokine_filter and pheno_filter_match
idx = intersect(rownames(cellCounts), basicPhenos$ID_500fg) 
length(idx) #458samples
cellCounts_filter = cellCounts[which(rownames(cellCounts) %in% idx),]
head(cellCounts_filter)  #458samples,1596 metabolite
dim(cellCounts_filter)

basicPhenos_filter = basicPhenos %>% filter(ID_500fg %in% idx)
head(basicPhenos_filter) #458samples
dim(basicPhenos_filter)

sum(rownames(cellCounts_filter) == basicPhenos_filter$ID_500fg)

basicPhenos_filter_match = basicPhenos_filter[match(rownames(cellCounts_filter), basicPhenos_filter$ID_500fg), ]

sum(basicPhenos_filter_match$ID_500fg == rownames(cellCounts))
head(basicPhenos_filter_match) 
dim(basicPhenos_filter_match) #458samples
# write.csv(cellCounts_filter, file = "/vol/projects/yzhang/500FG_aging/input/cellCounts/common_cellCounts_filter_age.csv")
# write.csv(basicPhenos_filter_match, file = "/vol/projects/yzhang/500FG_aging/input/cellCounts/common_basicPhenos_filter_match_cellCounts_age.csv")

cellCounts_age <- cbind(cellCounts_filter, Age = basicPhenos_filter_match$Age)
head(cellCounts_age)
dim(cellCounts_age)

##replaced the Inf value
# Count the number of infinite values in each column
inf_counts <- sapply(cellCounts_age, function(col) sum(is.infinite(col)))
inf_counts[inf_counts > 0] 

# Replace infinite values with the maximum finite value in each numeric column
cellCounts_age <- lapply(cellCounts_age, function(col) {
  if (is.numeric(col)) {  # Ensure the column is numeric
    col[is.infinite(col)] <- max(col[!is.infinite(col)], na.rm = TRUE)  # Replace Inf
  }
  return(col)
})

# Convert back to a data frame 
cellCounts_age <- as.data.frame(cellCounts_age)
rownames(cellCounts_age) <- rownames(cellCounts_filter)

head(cellCounts_age)

##perform 100 times
# Custom function to train and evaluate the model
train_and_evaluate <- function(data, seed, iteration) {
  set.seed(seed)
  
  # Split data into training and validation sets
  splitSample <- createDataPartition(data$Age, p = 0.7, list = FALSE)
  training <- data[splitSample, ]
  validation <- data[-splitSample, ]
   
    # Set up the grid for hyperparameter tuning
  netGrid <- expand.grid(.alpha = seq(0.1, 0.9, by = 0.1), .lambda = seq(0, 1, by = 0.01))    
  netctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5) 
  
    # Train the model
  netFit <- train(Age ~ ., data = training, method = "glmnet", metric = "RMSE",
                  tuneGrid = netGrid, trControl = netctrl, preProcess = "scale")
  
  # Make predictions on training and validation sets
  train_predictions <- predict(netFit, newdata = training)
  val_predictions <- predict(netFit, newdata = validation)
  
    # Save predicted age for each sample   
   pred_df <- rbind(
    data.frame(
      iteration = iteration,
      set = "train",
      sample_id = rownames(training),
      Age_true = training$Age,
      Age_pred = as.numeric(train_predictions),
      stringsAsFactors = FALSE
    ),
    data.frame(
      iteration = iteration,
      set = "validation",
      sample_id = rownames(validation),
      Age_true = validation$Age,
      Age_pred = as.numeric(val_predictions),
      stringsAsFactors = FALSE
    )
  )
  # Compute performance metrics for training set
  train_rmse <- sqrt(mean((training$Age - train_predictions)^2))  # RMSE
  train_mse <- mean((training$Age - train_predictions)^2)         # MSE
  train_r_squared <- cor(training$Age, train_predictions)^2       # R²
  
  # Compute performance metrics for validation set
  val_rmse <- sqrt(mean((validation$Age - val_predictions)^2))    # RMSE
  val_mse <- mean((validation$Age - val_predictions)^2)           # MSE
  val_r_squared <- cor(validation$Age, val_predictions)^2         # R²
  
  # Return metrics
  list(
    train_RMSE = train_rmse,
    train_MSE = train_mse,
    train_R2 = train_r_squared,
    val_RMSE = val_rmse,
    val_MSE = val_mse,
    val_R2 = val_r_squared,
    best_model = netFit,
    # selected_features = selected_features,
    predictions = pred_df
  )
}


# Run the training process 100 times
set.seed(1)
results <- lapply(1:100, function(i) train_and_evaluate(cellCounts_age, seed = i, iteration = i))

# Extract metrics for training and validation sets
train_rmse_values <- sapply(results, function(x) x$train_RMSE)
train_mse_values <- sapply(results, function(x) x$train_MSE)
train_r2_values <- sapply(results, function(x) x$train_R2)

val_rmse_values <- sapply(results, function(x) x$val_RMSE)
val_mse_values <- sapply(results, function(x) x$val_MSE)
val_r2_values <- sapply(results, function(x) x$val_R2)

# Calculate mean and standard deviation for each metric
metrics_summary <- data.frame(
  Metric = c("Train RMSE", "Validation RMSE", 
             "Train MSE", "Validation MSE", 
             "Train R²", "Validation R²"),
  Mean = c(mean(train_rmse_values), mean(val_rmse_values),
           mean(train_mse_values), mean(val_mse_values),
           mean(train_r2_values), mean(val_r2_values)),
  SD = c(sd(train_rmse_values), sd(val_rmse_values),
         sd(train_mse_values), sd(val_mse_values),
         sd(train_r2_values), sd(val_r2_values))
)

# Print summary
head(metrics_summary)

saveRDS(results, "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet_orign_sample_feature/cellcount_prediction_elastic_results.rds")

all_predictions <- do.call(rbind, lapply(results, function(x) x$predictions))
write.csv(all_predictions, file = "/vol/projects/yzhang/500FG_aging/output/12_prediction/elasticnet_orign_sample_feature/cellcount_pre_age_elasticnet_orign_sample_feature.csv", row.names = FALSE)

# Combine R² metrics into a single data frame for plotting
plot_data_r2 <- data.frame(
  Set = rep(c("Training", "Validation"), each = 100),
  R2 = c(train_r2_values, val_r2_values)
)

# Plot only R²
ppi <- 300
# png("/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_cross_validation.png", 
#     width = 5 * ppi, height = 4 * ppi, res = ppi)

g<-ggplot(plot_data_r2, aes(x = Set, y = R2, fill = Set)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(position = position_jitterdodge(), size = 1, alpha = 0.6, color = "black") +  # Scatter points
  labs(title = "R² Distribution Across Training and Validation in Cellcounts Sets",
       x = "Cellcounts",
       y = "R²") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

print(g)
dev.off()
print(g)

save.image("/vol/projects/yzhang/500FG_aging/output/07_cellcounts/cellcounts_elasticnet_orign_sample.Rdata")
