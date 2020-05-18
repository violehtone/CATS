# Author: Chao (Cico) Zhang
# Date: 31 Mar 2017
# Usage: Rscript run_model.R -i unlabelled_sample.txt -m model.pkl -o output.txt
# If you are using python, please use the Python script template instead.
# Set up R error handling to go to stderr
options(show.error.messages=F, error=function(){cat(geterrmessage(),file=stderr());q("no",1,F)})

# Import required libraries
# You might need to load other packages here.
suppressPackageStartupMessages({
  library('getopt')
  library('caret')
})

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get options using the spec as defined by the enclosed list
# Read the options from the default: commandArgs(TRUE)
option_specification <- matrix(c(
  'input', 'i', 2, 'character',
  'model', 'm', 2, 'character',
  'output', 'o', 2, 'character'
), byrow=TRUE, ncol=4);

# Parse options
options <- getopt(option_specification);

# Start your coding
## FOR DEBUGGING PURPOSES
#my_model <- readRDS("final_model1.rds")
#Train <- read.delim("Validation_call.txt",header=FALSE,sep = "\t")
#output <- "output.txt"

my_model <- readRDS(options$model)
Train <- read.delim(options$input,header=FALSE,sep = "\t")

Train <- t(Train[,-(1:4)])
Train <- as.data.frame(Train)

# Rename column names to keep chromosome + position information
head(Train[, 1:10])

feature_information <- read.delim(options$input, header=F, sep="\t")[, 1:3]
#feature_information <- read.delim("Validation_call.txt",header=FALSE,sep = "\t")[, 1:3]

# head(feature_information)

for(i in seq_along(names(Train[-1]))){
  new_col_name <- toString(paste("C", feature_information[i + 1, 1], "S", feature_information[i + 1, 2], "E", feature_information[i + 1, 3], sep = ""))
  colnames(Train)[i + 1] <- new_col_name
}

#print(names(Train))

Train2 <- Train[,-1]
rownames(Train2) <- Train[,1]
Train <- Train2

# Convert dataset values to numeric
for (i in 2:ncol(Train)) {
  Train[[i]] <- as.numeric(as.character(Train[[i]]))
}

#predictions <- predict(my_model, newdata = Train[,-1])

predictions <- predict(my_model, newdata = Train[,-1])

predictions_df <- data.frame(rownames(Train),predictions)

colnames(predictions_df) <- c("Sample", "Subgroup")

write.table(predictions_df, options$output, append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = TRUE)

# suggested steps
# Step 1: load the model from the model file (options$model)
# Step 2: apply the model to the input file (options$inout) to do the prediction
# Step 3: write the prediction into the desinated output file (options$output)

# End your coding
message ("Done!")
