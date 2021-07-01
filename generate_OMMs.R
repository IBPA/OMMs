#!/usr/bin/env Rscript
#' Design mixed meals to maximize the information content of glycan profiles.

suppressPackageStartupMessages(library("argparse"))
library(MaxPro)
library(readxl)
library(stringr)
library(gurobi)

source("./src/common.R")
DEV_MODE <- TRUE

# Parse and return command line arguments of this script (return defaults if DEV_MODE is TRUE).
get_args <- function() {
  parser <- ArgumentParser(description = "Design mixed meals to maximize the information content of glycan profiles.")
  parser$add_argument("--glycanDB", required=TRUE,
                      help="A .csv file containing glycan content profiles of individual foods when dried.")
  parser$add_argument("--moistureDB", required=TRUE,
                      help="A .xlsx file containing the moisture percentages of individual foods before drying.")
  parser$add_argument("--max-meal-ingredients", type="integer",required=FALSE,
                      help="The maximum number of ingredients in each designed mixed meal.")
  parser$add_argument("--compositional", action='store_true', required=FALSE,
                      help="Transform glycans to compositional values.")
  parser$add_argument("--num-meals", type="integer", required=TRUE,
                      help="The number of mixed meals to design.")
  parser$add_argument("--output", required=TRUE,
                      help="The output .csv filename that will contain the designed OMMs where each mixed meal is defined by the individual food proportions that it includes.")
  
  testargs <- c("--glycanDB", "./data/draft_mono_for_ameen_060321.csv", 
                "--moistureDB", "./data/% moisture content for all foods.xlsx", 
                "--num-meals", 50,
                "--max-meal-ingredients", 5,
                "--compositional",
                "--output", "./results/gen_OMMS.csv")
  if (DEV_MODE){
    return(parser$parse_args(testargs))
  }
  return(parser$parse_args())
}

# Given the 'df_food_glycans', find a mixed meal (i.e. proportions of each food),
# with a glycan content that is most similar to the 'target_glycans' specified.
get_proportions <- function(df_food_glycans, target_glycans){
  # 1) initialize the dimensions
  n <- nrow(df_food_glycans) # Number of foods
  d <- ncol(df_food_glycans) # Number of glycans
  
  # 2) Define the linear programming (LP) problem
  # 2.1) Construct the 'a' vector used in the LP objective.
  a <- c(rep(0, n), rep(1, d), rep(1, d))
  
  # 2.2) Construct the 'G' matrix used in the LP equality constraints.
  A <- t(df_food_glycans)
  G <- cbind(A, diag(1, d, d), diag(-1, d, d))
  G <- rbind(G, c(rep(1, n), rep(0, d), rep(0, d)))
  
  # 2.3) Construct the 'g' vector used in the right-hand-side of LP equality constraints.
  g <- c(t(target_glycans), 1)
  
  # 2.4) Construct the 'l' vector used as the lower-bound for the variable constraints.
  l <- rep(0, n + 2*d)
  
  # 3) Setup the Gurobi model based on the LP that is defined above.
  model = list()
  model$A <- G
  model$rhs <- g
  model$sense <- "="
  model$obj <- a
  model$lb <- l
  model$modelsense <- "min"
  
  # 4) Solve the defined LP using Gurobi
  result <- gurobi(model)
  
  # 5) Return the LP solution. Note, we use the first portion (i.e. 1:n) as the rest represent the minimized errors.
  return(result$x[1:n])
}

args <- get_args()

# 1) Load food data
df <- load_food_data(args$glycanDB, args$moistureDB)
if(args$compositional){
  df[,MONO_COLUMNS] <- df[,MONO_COLUMNS]/rowSums(df[,MONO_COLUMNS])
  df[is.na(df)] <- 0
} else {
  df[,MONO_COLUMNS] <- minmax_normalize(df[,MONO_COLUMNS])
}
df_food_glycans <- df[,MONO_COLUMNS]
rownames(df_food_glycans) <- df$uid

# 2) If max-meal-ingredients is set, include only foods with maximum difference in their glycan content
if(length(args$max_meal_ingredients)){
  # 2.1) Reorder foods to maximize the difference of their glycan content
  D_new <- MaxProRunOrder(df_food_glycans)
  df_new_D <- data.frame(D_new$Design[1:args$max_meal_ingredients,]) # keep only the top foods
  df_new_D <- df_new_D[df_new_D$X1,] # Ensure the order is correct
  df_new_D$X1 <- NULL # Remove first column as it determines the order only
  colnames(df_new_D) <- MONO_COLUMNS
  
  # 2.2) Identify and assign the correct food uids as rownames
  design_vec <- unname(apply(df_new_D, 1, function(x){paste(x, collapse = ",")}))
  all_vec <- unname(apply(df_food_glycans, 1, function(x){paste(x, collapse = ",")}))
  idx_selected_meals <- match(design_vec, all_vec)
  df_food_glycans <- df_food_glycans[idx_selected_meals,]
}


# 3) Generate design using MaxPro. A design here is the hypothetical glycan vector.
#    The goal here is to generate a set of glycan vectors that are very different from each other.
# 3.1) Ignore columns that have 0 variance and reorder
non_const_columns <- MONO_COLUMNS[(sapply(df_food_glycans, var) != 0)]
# 3.2) Create an initial design (generate more than needed), ignore const columns.
D_init <- MaxProLHD(args$num_meals*2, length(non_const_columns))$Design
# 3.3) Optimize the initialized design.
D <- MaxPro(D_init)

# 4) Select best design
D_all <- MaxProRunOrder(D$Design) # Select N new designs given the data
df_new_D <- data.frame(D_all$Design)
df_new_D <- df_new_D[df_new_D$X1,] # Ensure the order is correct
df_new_D$X1 <- NULL # Remove first column as it determines the order only
colnames(df_new_D) <- non_const_columns
df_design <- data.frame(matrix(0, nrow(df_new_D), length(MONO_COLUMNS)))
colnames(df_design) <- MONO_COLUMNS
df_design[,non_const_columns] <- df_new_D[,non_const_columns]

# 5) Find a mixed-meal (i.e. food proportions) for each design considering that
#    the glycan content of a mixed-meal is additive with respect to the 
#    proportions of its constituent foods.
df_food_proportions <- data.frame(matrix(0, nrow(df_design), nrow(df_food_glycans)))
colnames(df_food_proportions) <- rownames(df_food_glycans)
for(i in 1:nrow(df_design)){
  df_food_proportions[i,] <- get_proportions(df_food_glycans, df_design[i,])
}

# 6) Calculate the glycan content of the designed mixed-meals above to select 
#    the best set of mixed-meals amongst them. Note that the glycan content of
#    a given mixed meal, may not exactly add up to the designed target glycan
#    content since the corresponding LP may have non-zero error. 
# 6.A) Calculate the glycan vectors of mixed meals, then reorder accordingly to select the top ones
MM_vectors <- data.matrix(df_food_proportions) %*% data.matrix(df_food_glycans[colnames(df_food_proportions), MONO_COLUMNS])
MM_vectors_reorder <- MaxProRunOrder(MM_vectors)
df_new_D <- data.frame(MM_vectors_reorder$Design)
df_new_D$X1 <- NULL # Remove first column as it determines the order only
colnames(df_new_D) <- non_const_columns
df_design <- data.frame(matrix(0, nrow(df_new_D), length(MONO_COLUMNS)))
colnames(df_design) <- MONO_COLUMNS
df_design[,non_const_columns] <- df_new_D[,non_const_columns]

# 6.B) Find food proportions for each design (considering additive model).
#      TODO: there is alternative way, it is done like this by resolving the LPs
#            for convenience of coding.
df_food_proportions <- data.frame(matrix(0, args$num_meals, nrow(df_food_glycans)))
colnames(df_food_proportions) <- rownames(df_food_glycans)
for(i in 1:args$num_meals){
  df_food_proportions[i,] <- get_proportions(df_food_glycans, df_design[i,])
}

# 7) Save the designed mixed meals
save_MMs(df_food_proportions, args$output)

