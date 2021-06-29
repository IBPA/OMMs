#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
library(MaxPro) # MaxPro, is considered a space filling design
library(readxl)
library(stringr)
library(gurobi)

source("./src/common.R")

get_args <- function() {
  parser <- ArgumentParser(description = "Design mixed meals to maximize the information content of glycan profiles.")
  parser$add_argument("--glycanDB", required=TRUE,
                      help="A .csv file containing glycan content profiles of individual foods when dried.")
  parser$add_argument("--moistureDB", required=TRUE,
                      help="A .xlsx file containing the moisture percentages of individual foods before drying.")
  parser$add_argument("--max-meal-ingredients", type="integer",required=FALSE,
                      help="The maximum number of ingredients in each designed mixed meal.")
  parser$add_argument("--num-meals", type="integer", required=TRUE,
                      help="The number of mixed meals to design.")
  parser$add_argument("--output", required=TRUE,
                      help="The output .csv filename that will contain the designed OMMs where each mixed meal is defined by the individual food proportions that it includes.")
  
  testargs <- c("--glycanDB", "./data/draft_mono_for_ameen_060321.csv", 
                "--moistureDB", "./data/% moisture content for all foods.xlsx", 
                "--num-meals", 50, 
                "--output", "./results/gen_OMMS.csv")
  args <- parser$parse_args(testargs)
  
  return(args)
}

get_proportions <- function(df_data, target){
  n <- nrow(df_data)
  d <- ncol(df_data)
  
  a <- c(rep(0, n), rep(1, d), rep(1, d))
  
  A <- t(df_data)
  G <- cbind(A, diag(1, d, d), diag(-1, d, d))
  G <- rbind(G, c(rep(1, n), rep(0, d), rep(0, d)))
  
  g <- c(t(target), 1)
  
  l <- rep(0, n + 2*d)
  
  model = list()
  model$A <- G
  model$rhs <- g
  model$sense <- "="
  model$obj <- a
  model$lb <- l
  model$modelsense <- "min"
  
  result <- gurobi(model)
  
  return(result$x[1:n])
}

args <- get_args()

# 1) Load data
df <- load_data(args$glycanDB, args$moistureDB)
df[,MONO_COLUMNS] <- minmax_normalize(df[,MONO_COLUMNS])
df_data <- df[,MONO_COLUMNS]
rownames(df_data) <- df$uid

# 2) Generate design
non_const_columns <- MONO_COLUMNS[(sapply(df_data, var) != 0)]
dim <- length(non_const_columns) # ignore columns that have 0 variance
D_init <- MaxProLHD(args$num_meals*2, dim)$Design # Generate more than needed
D <- MaxPro(D_init)

# 3) Select best design
D_all <- MaxProRunOrder(D$Design) # Select N new designs given the data
df_new_D <- data.frame(D_all$Design)
df_new_D <- df_new_D[df_new_D$X1,] # Ensure the order is correct
df_new_D$X1 <- NULL # Remove first column as it determines the order only
colnames(df_new_D) <- non_const_columns
df_design <- data.frame(matrix(0, nrow(df_new_D), length(MONO_COLUMNS)))
colnames(df_design) <- MONO_COLUMNS
df_design[,non_const_columns] <- df_new_D[,non_const_columns]

# 4) Find food proportions for each design (considering additive model)
df_food_proportions <- data.frame(matrix(0, nrow(df_design), nrow(df)))
colnames(df_food_proportions) <- df$uid
for(i in 1:nrow(df_design)){
  df_food_proportions[i,] <- get_proportions(df_data, df_design[i,])
}

# 5.A) Calculate the glycan vectors of mixed meals, then reorder accordingly to select the top ones
MM_vectors <- data.matrix(df_food_proportions) %*% data.matrix(df_data[colnames(df_food_proportions), MONO_COLUMNS])
MM_vectors_reorder <- MaxProRunOrder(MM_vectors)
df_new_D <- data.frame(MM_vectors_reorder$Design)
df_new_D$X1 <- NULL # Remove first column as it determines the order only
colnames(df_new_D) <- non_const_columns
df_design <- data.frame(matrix(0, nrow(df_new_D), length(MONO_COLUMNS)))
colnames(df_design) <- MONO_COLUMNS
df_design[,non_const_columns] <- df_new_D[,non_const_columns]

# 5.b) Find food proportions for each design (considering additive model)
df_food_proportions <- data.frame(matrix(0, args$num_meals, nrow(df)))
colnames(df_food_proportions) <- df$uid
for(i in 1:args$num_meals){
  df_food_proportions[i,] <- get_proportions(df_data, df_design[i,])
}

# 6) Save
save_MMs(df_food_proportions, args$output)
  
