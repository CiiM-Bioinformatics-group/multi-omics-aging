library(pROC)
library(caret)
library(dplyr)
library(glmnet)
library(readr)
library(ggplot2)

R.home()
system("which Rscript", intern = TRUE)
system("which -a Rscript", intern = TRUE)


genotype <- read.csv("/vol/projects/yzhang/500FG_aging/input/variance/cumulative_variance_input/genotype_pca_scores.csv")
rownames(genotype) <- genotype[,1]
genotype <- genotype[,-1]
head(genotype)
dim(genotype)

basicPhenos = read.csv("/vol/projects/yzhang/500FG_aging/output/02_Olink/Age_group_basicPhenos.csv")
head(basicPhenos)
# prepare cytokine_filter and EAA_filter_match
idx = intersect(rownames(genotype), basicPhenos$ID_500fg) 
length(idx) #458samples
genotype_filter = genotype[which(rownames(genotype) %in% idx),]
head(genotype_filter)  #458samples,1596 metabolite
dim(genotype_filter)

basicPhenos_filter = basicPhenos %>% filter(ID_500fg %in% idx)
head(basicPhenos_filter) #458samples
dim(basicPhenos_filter)

sum(rownames(genotype_filter) == basicPhenos_filter$ID_500fg)

basicPhenos_filter_match = basicPhenos_filter[match(rownames(genotype_filter), basicPhenos_filter$ID_500fg), ]

sum(basicPhenos_filter_match$ID_500fg == rownames(genotype_filter))
head(basicPhenos_filter_match) 
dim(basicPhenos_filter_match) #458samples

genotype_age <- cbind(genotype_filter, Age = basicPhenos_filter_match$Age)
head(genotype_age)
dim(genotype_age)

#Age
#split data to trianing and testing set
set.seed(1)
splitSample <- createDataPartition(genotype_age$Age, p = 0.7, list = FALSE)
training <- genotype_age[splitSample,]
validation <- genotype_age[-splitSample,]

#settings grid tuning and cv
netGrid <- expand.grid(.alpha=1,.lambda=seq(0,1,by=0.01))  
netctrl <- trainControl(method="repeatedcv", number=10, repeats=5) 

#training
netFit <- train(Age~.,data = training, method = 'glmnet', metric = "RMSE",tuneGrid=netGrid,trControl = netctrl,preProcess = "scale")

#look at the model
head(arrange(netFit$results, RMSE) )

netFit$bestTune

my.glmnet.model<- netFit$finalModel
aa<- as.matrix(coef(my.glmnet.model,s=netFit$bestTune$lambda))
b<-as.data.frame(aa[which(aa[,1]!=0),])
colnames(b)<-"coef"
b

#Prediction on validation set
predictions <- predict(netFit, validation)

cor.test(predictions,validation$Age)

RMSE(predictions,validation$Age)

plot(predictions, validation$Age,
     xlab = "predictions", 
     ylab = "validation$Age",
     main = "Scatter Plot of Predictions vs Validation Age",
     pch = 1)  # pch = 1 present circle

##perform 100 times
# Custom function to train and evaluate the model
train_and_evaluate <- function(data, seed) {
  set.seed(seed)
  
  # Split data into training and validation sets
  splitSample <- createDataPartition(data$Age, p = 0.7, list = FALSE)
  training <- data[splitSample, ]
  validation <- data[-splitSample, ]
  
  # Set up the grid for hyperparameter tuning
  netGrid <- expand.grid(.alpha = 1, .lambda = seq(0, 1, by = 0.01))  
  netctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5) 
  
  # Train the model
  netFit <- train(Age ~ ., data = training, method = "glmnet", metric = "RMSE",
                  tuneGrid = netGrid, trControl = netctrl, preProcess = "scale")
  
  # Make predictions on training and validation sets
  train_predictions <- predict(netFit, newdata = training)
  val_predictions <- predict(netFit, newdata = validation)
  
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
    best_model = netFit
  )
}

# Run the training process 100 times
set.seed(1)
results <- lapply(1:100, function(i) train_and_evaluate(genotype_age, seed = i))

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

# Combine R² metrics into a single data frame for plotting
plot_data_r2 <- data.frame(
  Set = rep(c("Training", "Validation"), each = 100),
  R2 = c(train_r2_values, val_r2_values)
)

# Plot only R²
ppi <- 300
png("/vol/projects/yzhang/500FG_aging/output/04_genotype/genotype_cross_validation.png", 
    width = 5 * ppi, height = 4 * ppi, res = ppi)

g<-ggplot(plot_data_r2, aes(x = Set, y = R2, fill = Set)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(position = position_jitterdodge(), size = 1, alpha = 0.6, color = "black") +  # Scatter points
  labs(title = "R² Distribution Across Training and Validation in Genotype Sets",
       x = "Cellcounts",
       y = "R²") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

print(g)
dev.off()
print(g)

save.image("/vol/projects/yzhang/500FG_aging/output/04_genotype/genotype_prediction.Rdata")

### raw data with high variance
results<-readRDS("/vol/projects/yzhang/500FG_aging/input/methylation/prediction/genotype_prediction_results.rds")  

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

# Combine R² metrics into a single data frame for plotting
plot_data_r2 <- data.frame(
  Set = rep(c("Training", "Validation"), each = 100),
  R2 = c(train_r2_values, val_r2_values)
)

# Plot only R²
ppi <- 300
png("/vol/projects/yzhang/500FG_aging/output/04_genotype/genotype_cross_validation(raw data).png", 
    width = 5 * ppi, height = 4 * ppi, res = ppi)

g<-ggplot(plot_data_r2, aes(x = Set, y = R2, fill = Set)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(position = position_jitterdodge(), size = 1, alpha = 0.6, color = "black") +  # Scatter points
  labs(title = "R² Distribution Across Training and Validation in Genotype Sets",
       x = "Genotype",
       y = "R²") +
  theme_minimal() +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

print(g)
dev.off()
print(g)                        

save.image("/vol/projects/yzhang/500FG_aging/output/04_genotype/genotype_prediction_raw.Rdata")
