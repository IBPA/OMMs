#!/usr/bin/env Rscript
#' Select a set of OMMs from candidate mixed meals, given a set of approved mixed meals.

suppressPackageStartupMessages(library("argparse"))
library(infotheo)
library(ggplot2)
library(reshape2)
library(scales)

source("./src/common.R")
source("./src/theme_util.R")
DEV_MODE <- TRUE

# Parse and return command line arguments of this script (return defaults if DEV_MODE is TRUE).
get_args <- function() {
  parser <- ArgumentParser(description = "Select a set of OMMs from candidate mixed meals, given a set of approved mixed meals.")
  parser$add_argument("--glycanDB", required=TRUE,
                      help="A .csv file containing glycan content profiles of individual foods when dried.")
  parser$add_argument("--moistureDB", required=TRUE,
                      help="A .xlsx file containing the moisture percentages of individual foods before drying.")
  parser$add_argument("--num-meals", type="integer", required=TRUE,
                      help="The number of mixed meals to select.")
  parser$add_argument("--candidate-meals", required=TRUE,
                      help="A .csv file containing the the set of candidate meals to select from.")
  parser$add_argument("--approved-meals", required=FALSE,
                      help="A .csv file containing the the set of approved meals.")
  parser$add_argument("--compositional", action='store_true', required=FALSE,
                      help="Transform glycans to compositional values.")
  parser$add_argument("--output", required=TRUE,
                      help="A .csv file containing the selected OMMs where each mixed meal is defined by the individual food proportions that it includes,")
  
  testargs <- c("--glycanDB", "./data/draft_mono_for_ameen_060321.csv", 
                "--moistureDB", "./data/% moisture content for all foods.xlsx", 
                "--candidate-meals", "./results/cand_OMMS.csv",
                "--approved-meals", "./results/appr_OMMS.csv", 
                "--num-meals", 10,
                "--compositional",
                "--output", "./results/sel_OMMs.csv")

  if (DEV_MODE){
    return(parser$parse_args(testargs))
  }
  return(parser$parse_args())
}

args <- get_args()

# 1) Load data
# 1.1) Load food data
df_food_vectors <- load_food_data(args$glycanDB, args$moistureDB)
rownames(df_food_vectors) <- df_food_vectors$uid
# 1.2) Transform food data
if(args$compositional){
  df_food_vectors[,MONO_COLUMNS] <- df_food_vectors[,MONO_COLUMNS]/rowSums(df_food_vectors[,MONO_COLUMNS])
  df_food_vectors[is.na(df_food_vectors)] <- 0
} else {
  df_food_vectors[,MONO_COLUMNS] <- minmax_normalize(df_food_vectors[,MONO_COLUMNS])
}
# 1.3) Load mixed-meal data
cand_MMs <- load_MMs(args$candidate_meals)
appr_MMs <- matrix(0,ncol = 1, nrow = 1, dimnames = list(1,df_food_vectors$uid[1])) # just an empty meal
if (length(args$approved_meals)){
  appr_MMs <- load_MMs(args$approved_meals)
}

# 2) Calculate corresponding glycan vectors
cand_MM_vectors <- cand_MMs %*% data.matrix(df_food_vectors[colnames(cand_MMs), MONO_COLUMNS])
appr_MM_vectors <- appr_MMs %*% data.matrix(df_food_vectors[colnames(appr_MMs), MONO_COLUMNS])

# 3) Pick from candidate vectors given the approved ones.
non_const_columns <- MONO_COLUMNS[(sapply(data.frame(cand_MM_vectors), var) != 0)]
D_new <- MaxProAugment(appr_MM_vectors[,non_const_columns], cand_MM_vectors[,non_const_columns], args$num_meals)
df_D_new <- data.frame(D_new$Design[(nrow(appr_MM_vectors)+1):nrow(D_new$Design),])
colnames(df_D_new) <- non_const_columns
df_design <- data.frame(matrix(0, nrow(df_D_new), length(MONO_COLUMNS)))
colnames(df_design) <- MONO_COLUMNS
df_design[,non_const_columns] <- df_D_new[,non_const_columns]

# 4) Identify the corresponding meals
design_vec <- unname(apply(df_design, 1, function(x){paste(x, collapse = ",")}))
cand_vec <- unname(apply(cand_MM_vectors, 1, function(x){paste(x, collapse = ",")}))
idx_selected_meals <- match(design_vec, cand_vec)
selected_MMs <- cand_MMs[idx_selected_meals,]

# 6) Save the selected mixed meals
save_MMs(selected_MMs, args$output)
