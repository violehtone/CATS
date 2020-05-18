#  LOAD DEPENDENCIES & SOURCE UTILITY FUNCTIONS
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("neuralnet")
library("caret")
library(tidyverse)
source("CATS_utility_functions_2.3.R")
library(dplyr)


# LOAD & PREPROCESS DATA
# Load data
Train <- read.delim("Train_call.txt",header=FALSE,sep = "\t")
Labels <- read.delim("Train_clinical.txt",header=TRUE,sep="\t")

#preprocess
Train <- t(Train[,-(1:4)])
Train <- as.data.frame(Train)
names(Train)[1] <- "Sample"

# Rename column names to keep chromosome + position information
feature_information <- read.delim("Train_call.txt", header=F, sep="\t")[, 1:3]
for(i in seq_along(names(Train[-1]))){
  new_col_name <- toString(paste("C", feature_information[i + 1, 1], "S", feature_information[i + 1, 2], "E", feature_information[i + 1, 3], sep = ""))
  colnames(Train)[i + 1] <- new_col_name
}

data <- merge(Labels, Train, by="Sample")
row.names(data) <- data$Sample
data$Sample <-NULL

# Convert dataset values to numeric
for (i in 2:ncol(data)) {
  data[[i]] <- as.numeric(as.character(data[[i]]))
}


# DEFINING VARIABLES
n_features <- 20      # Number of features to use
p_cutoff <- 0.0005     # P-cutoff value to filter features to use
cv_fold_inner <- 10   # Number of folds in inner loop
cv_fold_outer <- 10   # Number of folds in outer loop
names_of_models <- c("rpart",'cforest','svmLinear',"NNs") # Names of the models used
gene_features <- c("C17S35076296E35282086", "C11S98637002E101368726", "C6S152086678E152409827") # name of the important features that contain HER2, ESR1, and PGR genes
set.seed(42)
list_of_models <- create_nested_list(names_of_models, cv_fold_outer) # list to store models
list_of_pred <- create_nested_list(names_of_models, cv_fold_outer) # list to store predictions
list_of_prediction_scores <- create_nested_list(names_of_models, cv_fold_outer)  # list to store pred. scores
feature_set <- c() #list to store selected features
data <- partition(data, cv_fold_outer, "Subgroup") #partition used in the outer loop

# DOUBLE-CROSS-VALIDATION LOOP
# OUTER LOOP

for (i in 1:cv_fold_outer) {
  # Combines all folds except i-th into train set and i-th fold is used as validation
  train_set <- data[data$partition != i,]
  train_set$partition <- NULL
  validation_set <- data[data$partition == i, ]
  validation_set$partition <- NULL

  # Feature selection
  training_filtered <- fisher_feature_selection(train_set, 0, p_cutoff, gene_features)
  features <- colnames(training_filtered)
  feature_set <- c(feature_set, features)
  validation <- validation_set[ ,features]
  
  #partition in the inner loop
  training_filtered <- partition(training_filtered, cv_fold_inner, "Subgroup")
  
  
  # INNER LOOP
  # Training and inner CV - Neural network
  for (j in 1:cv_fold_inner) {
    train_9 <- training_filtered[training_filtered$partition != j, ]
    train_9$partition <- NULL
    
    ifelse(j >= 2, 
           model_fit <- NN_model(train_9, s = "Subgroup", weight = model_fit$weights),
           model_fit <- NN_model(train_9, s = "Subgroup")
    )}
  
  list_of_models[["NNs"]][[i]] <- model_fit
  
  # Training and inner CV - Other models (DT, RF, and SVM)
  training_filtered$partition <- NULL
  
  for (name in names_of_models[1:3]) {
    list_of_models[[name]][[i]] <- build_model(training_filtered, name, cv_fold_inner)
  }

  
  # PREDICTIONS
  # Neural networks
  pred <- predict(model_fit, validation)
  t <- table(validation$Subgroup, apply(pred, 1, which.max))
  pred <- apply(pred, 1, which.max)
  
  for (n in 1:3){
    pred[pred == n] <- rownames(t)[n]
  }
  
  list_of_pred[["NNs"]][[i]] <- pred
  
  # Other models (DT, RF, and SVM)
  real_class <- validation$Subgroup
  validation$Subgroup <- NULL
  for (name in names_of_models[1:3]){
    list_of_pred[[name]][[i]] <- predict(list_of_models[[name]][[i]], 
                                         newdata = validation,
                                         type = "raw")
  }
  
  
  # ASSESSING THE ACCURACY OF MODELS
  # Create a list of model accuracies'
  names <- c('obs', 'pred')
  for (name in names_of_models){
    predicted_data <- setNames(data.frame(real_class, list_of_pred[[name]][[i]]), names)
    summary <- defaultSummary(predicted_data, lev = 3)
    list_of_prediction_scores[[name]][[i]] <- summary
  }
}


# RESULTS
# Average accuracies
avg_accuracies <- data.frame(name = numeric(0),
                             mean_accuracy = numeric(0))

for (name in names_of_models){
  mean_accuracy<-mean(sapply(1:cv_fold_outer,
                             function(i){list_of_prediction_scores[[name]][[i]][1]}))
  
  # print and store the results
  print(paste("The average accuracy of ",name,"model is",mean_accuracy))
  avg_accuracies[nrow(avg_accuracies) + 1,] <- c(name, mean_accuracy)
}

# Find the top features/biomarkers
feature_table <- table(feature_set)
actual_features <- feature_table[-nrow(feature_table)] #assuming that 'Subgroup' is the last column
actual_features <- sort(actual_features, decreasing = FALSE)
actual_features <- as.data.frame(actual_features)
actual_features <- actual_features[actual_features$Freq >= 5,]


# VISUALIZING THE RESULTS
# Plot the importance of biomarkers
features0 <- ggplot(data = actual_features,
       mapping = aes(x = Freq,
                     y = feature_set,
                     fill = feature_set)) +
  geom_col() +
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10)) +
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 8))+
  labs(x = "Frequency", y = "Feature", title = "Frequency of features selected")


# Plot the avg. accuracies of each method
avg_accuracies$mean_accuracy <- as.numeric(avg_accuracies$mean_accuracy)

ggplot(data = avg_accuracies,
       mapping = aes(x = name,
                     y = mean_accuracy,
                     fill = name)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_discrete(name = "Method",
                      breaks = c("rpart", "cforest", "svmLinear", "NNs"),
                      labels = c("Decision tree", "Random Forest", "SVM", "Neural Network")) +
  expand_limits(y = 1) +
  labs(x = "method name", y = "avg. accuracy (%)", title = "Average accuracy for each method")


# Plot the accuracies vs. iterations
iteration_accuracies <- c()

# Create a list of prediction accuracies per iteration
for (method in names(list_of_prediction_scores)) {
  scores <- nested_list_to_tibble(list_of_prediction_scores[[method]])
  scores$method <- method
  if(length(iteration_accuracies) == 0) {
    iteration_accuracies <- scores
  }else{
    iteration_accuracies <- full_join(iteration_accuracies, scores)
  }
}

iteration_accuracies$iteration <- rep(c(1:10), 4)

ggplot(data = iteration_accuracies,
       mapping = aes(x = iteration,
                     y = Accuracy,
                     color = method)) +
  geom_line() +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  scale_color_discrete(name = "Method",
                      breaks = c("rpart", "cforest", "svmLinear", "NNs"),
                      labels = c("Decision tree", "Random Forest", "SVM", "Neural Network")) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10)) +
  labs(x = "Iteration",
       y = "Accuracy (%)",
       title = "Accuracies in each iteration")


# Calculate the variances between iterations in different models as a violin plot
print("The variation for DT, RF, SVM, and NN are respectively: ")
var(iteration_accuracies[1:10,]$Accuracy) #rpart variance
var(iteration_accuracies[11:20,]$Accuracy) #cforest variance
var(iteration_accuracies[21:30,]$Accuracy) #svmlinear variance
var(iteration_accuracies[31:40,]$Accuracy) #NNs variance

# Change the method names for iteration_accuracies table
iteration_accuracies$method <- ifelse(iteration_accuracies$method == 'rpart', 'DT',
                                      ifelse(iteration_accuracies$method == 'cforest', 'RF',
                                             ifelse(iteration_accuracies$method == 'svmLinear', "SVM",
                                                    ifelse(iteration_accuracies$method == 'NNs', 'NN', NA))))

# Plot the density of accuracy scores using a violin plot
ggplot(data = iteration_accuracies,
       mapping = aes(x = method,
                     y = Accuracy,
                     fill = method)) +
  geom_violin() +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  scale_fill_discrete(name = "Method",
                      breaks = c("DT", "RF", "SVM", "NN"),
                      labels = c("Decision tree", "Random Forest", "SVM", "Neural Network")) +
  labs(x = "Classification method", 
       y = "accuracy (%)", 
       title = "Density of prediction accuracies for each method")


# FEATURES PER SUBGROUP
# Split the important biomarkers based on cancer subtype
most_important_features <- c(as.character(actual_features$feature_set), "Subgroup")
data_filtered <- data[, most_important_features]

head(data_filtered[,1:length(data_filtered)])

feature_result <- lapply(apply(data_filtered[1:(length(data_filtered)-1)], 2, function(x) {
  tapply(x, data_filtered[[(length(data_filtered))]], function(x) {
    c(mean(x), sd(x))
  })
}), function(x) {
  do.call(rbind, x)
})

# Transform nested list of matrices to a single dataframe
df_features <- t(data.frame(feature_result))
# Filter std. dev. rows out (leaving only mean)
sd_rows <- seq(0, (2 * (length(data_filtered) - 1)), 2)
df_features <- df_features[-sd_rows, ]
df_features <- as.data.frame(df_features)

# Wrangle dataframe
library(data.table)
df_features <- data.frame(setDT(df_features, keep.rownames = "DNA_region"))

df_features_new <- gather(data = df_features,
                          key = subgroup,
                          value = mean_aberration,
                          -DNA_region)

# Plot the means of features
features1 <- ggplot(data = df_features_new,
       mapping = aes(x = DNA_region,
                     y = mean_aberration,
                     fill = subgroup)) +
  geom_col(position = 'dodge') +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 8)) +
  labs(x = "DNA region",
       y = "Mean aberration",
       title = "Mean aberration in different DNA regions by breast cancer subgroups")


require(gridExtra)
grid.arrange(features0, features1, ncol=2)


# Find out the best model

## 1)  Save the best performing method
varlist = c()
for (i in seq(1, 40, 10)){varlist = append(varlist, var(iteration_accuracies[as.numeric(i):as.numeric(i+9),]$Accuracy))}
names(varlist) = names_of_models

best_score = 0
best_model = ""
new_score = 0
for (model in names_of_models){
  new_score = avg_accuracies[avg_accuracies$name == model,]$mean_accuracy / as.numeric(varlist[model])
  if (as.numeric(new_score) > best_score){
    best_score = new_score
    best_model = model
  }
}


cat("The best model is", best_model, "with a score (acc. / var) of", best_score)
list_of_models$best_model
#-> rpart is the best method

## Find the best model of the 10 rpart models
best_model = 0
best_acc = 0

for (i in 1:10) {
  # Calculate the avg. accuracy
  avg_acc = mean(list_of_models$rpart[[i]]$results$Accuracy)
  # Update best accuracy if necessary
    if(avg_acc >= best_acc) {
    best_acc = avg_acc
    best_model = i
  }
}
  
saveRDS(list_of_models$rpart[[2]], "final_model1.rds")
