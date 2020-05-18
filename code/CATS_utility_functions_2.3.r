# UTILITY FUNCTIONS
library("neuralnet")
library("caret")


create_nested_list <- function(names, list_length){
  #' Creates a nested list of length n with objects in 'names'
  result <- lapply(names,
                   function(i){vector(mode = "list", length = list_length)})
  names(result) <- names
  return(result)
}


fisher_feature_selection <- function(data, n=20, p_cutoff=0.0005, gene_features) {
  #' Selects top n most important features from data based on Fisher's Exact test
  #' data = data used for feature selection (class label should be the first feature column)
  #' n = amount of features to be selected (default = 20)
  #' p_cutoff = P-value cutoff value for the features to be selected (default = 0.0005)
  #' gene_features = list of features that will be added unconditionally
  #' 
  #' Note: features will be selected based on the p_cutoff value only if n = 0 
  df_features <- data[-1]
  df_classes <- data[1]
  p_values <- c()
  
  for (i in names(data)[-1]) {
    # Create a contingency table
    cont_tbl <- table(data[[i]], df_classes[[1]])
    # Perform Fisher's exact test
    #set.seed(42)
    test_result <- fisher.test(cont_tbl)
    # Add the p-value of the test to the list
    p_values[i] <- test_result$p.value
  }
  
  if(n>0) {
    # Select top n features
    filtered <- sort(p_values)[1:n]
  }else{
    filtered <- which(sort(p_values) < p_cutoff)
    cat(length(filtered), "additional features selected using Fisher's Exact test with a P-value of", p_cutoff, "\n")
  }
  
  # Add 3 important features to the list
  features_selected <- names(filtered)
  gene_feature_intersection <- intersect(gene_features, features_selected)
  gene_features <- gene_features[!(gene_features %in% gene_feature_intersection)]
  features_selected <- features_selected[!(features_selected %in% gene_feature_intersection)]
  
  final_features <- c(gene_feature_intersection, gene_features, features_selected)
  
  # Combine classes and selected features into a single df
  result <- cbind(df_classes, data[final_features])
  
  return(result)
}

#function of stratified partition
partition <- function(data, n, label) {
  #' Perform stratified partition on data
  #' x = dataset, n = number of folds, f = name of the label column in data
  data$partition <- NA
  groupsizes <- table(data[[label]])
  groupnames <- names(groupsizes)
  for (name in groupnames) {
    # The first call to 'sample' randomly permutes the row index numbers in a group
    # The second call to 'sample' adds one random value from 1 to n to the index number.
    # This prevents the first groupsize %% n partitions to always get the remaining rows.
    #set.seed(42)
    labels <- (sample(1:groupsizes[name]) + sample(1:n,1))%%n + 1
    data$partition[data[[label]]==name] <- labels
  }
  return(data)
}


NN_model <- function(df, s, weight=NULL) {
  #' Function for building a neural network model (using neuralnet package)
  f <- as.formula(paste(paste(s,"~"),paste(names(df)[-1],collapse=" + ")))
  
  #model
  #set.seed(42)
  NNs <- neuralnet(f,df,hidden=c(8,8),
                   act.fct = "logistic",err.fct="ce",learningrate=0.001,
                   algorithm="backprop",linear.output = FALSE,
                   startweights = weight,threshold=0.1)
  return(NNs)
}


build_model <- function(data, name, n=10) {
  #' function for building building a model (using caret package)
  #' data = training data, model = name of the model used, n = number of folds in CV (default = 10)
  
  ctrl <- trainControl(method = "repeatedcv", repeats = n)
  train_data <- data[, 2:ncol(data)]
  train_subgroups <- data[, 1]
  
  # Set grid for grid search
  nr_feats = length(names(train_data))
  if (name == 'rpart')
    {model_fit <- train(train_data,
                        train_subgroups,
                        method = name,
                        trControl = ctrl,
                        tuneLength = 10)}
  else {
    if (name == 'cforest'){grid <- expand.grid(mtry=c(nr_feats, sqrt(nr_feats), log(nr_feats)))}
    if (name == 'svmLinear'){grid <- expand.grid(C=seq(0.05, 2, 0.15))}
    
    
    # train model
    #set.seed(42)
    model_fit <- train(train_data,
                       train_subgroups,
                       method = name,
                       trControl = ctrl,
                       tuneGrid = grid)}
  return(model_fit) 
}


depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)

nested_list_to_tibble <- function(l) {
  #' Transforms a nested list into a tibble object
  if (depth(l) == 2) {
    return(bind_rows(l))
  } else {
    l <- map_depth(l, depth(l) - 2, bind_rows)
    nested_list_to_tibble(l)
  }
}
