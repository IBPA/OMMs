#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
library(MaxPro) # MaxPro, is considered a space filling design
library(readxl)
library(stringr)
library(gurobi)

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
                "--num-meals", 10, 
                "--output", "./results/gen_OMMS.csv")
  args <- parser$parse_args(testargs)
  
  return(args)
}

load_data <- function(data_file, metadata_file){
  # A) Load data
  df <- read.table(data_file, sep = ",", quote = '"', header = TRUE)
  df_data <- df[,MONO_COLUMNS]
  df_data[is.na(df_data)] <- 0
  df[,MONO_COLUMNS] <- df_data
  colnames(df)[4] <- "uid"
  
  # B) Load water content metadata
  df_meta <- read_excel(metadata_file)
  df_meta <- df_meta[,c("Unique Identifier Code","% water content")]
  colnames(df_meta) <- c("uid", "water")
  df_meta$uid <- as.numeric(str_remove(df_meta$uid, "^0+"))
  df_meta$dry <- (100-df_meta$water)/100
  df_meta$water <- NULL
  
  # C) Incorporate dry content % of each food into the corresponding glycan content vector
  df_comb <- merge(df, df_meta, by = "uid")
  df_comb[,MONO_COLUMNS] <- df_comb[,MONO_COLUMNS] * df_comb$dry
  #  ***** Skip this for now, due to issue in water % (not between 0-100)
  
  # D) Min-Max Normalization
  scale_minmax <- function(x){
    (x- min(x)) /(max(x)-min(x))
  }
  df_data_scaled <- as.data.frame(lapply(df_data, scale_minmax))
  df[,MONO_COLUMNS] <- df_data_scaled
  df[is.na(df)] <- 0
  
  
  return(df[,c("uid", MONO_COLUMNS)])
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

save_MMs <- function(df_food_proportions, filename) {
  # Write the non-zero proportions to the file for each meal.
  nz_proportions <- apply(df_food_proportions, 1, FUN = function(x){x[x!=0]})
  
  lines <- c()
  for (meal_idx in 1:length(nz_proportions)) {
    cur_meal_props <- nz_proportions[[meal_idx]]
    line_str <- paste(sprintf("%s:%.2f", names(cur_meal_props), unname(cur_meal_props)), collapse = ",")
    lines <- c(lines, line_str)
  }
  writeLines(lines, filename)
}


args <- get_args()

# Constants
MONO_COLUMNS <- c('Glucose', 'Galactose', 'Fructose', 'Xylose', 'Arabinose', 'Fucose', 'Rhamnose', 'GlcA', 'GalA', 'GlcNAc', 'GalNAc', 'Mannose', 'Allose', 'Ribose')

# 1) Load data
df <- load_data(args$glycanDB, args$moistureDB)
df_data <- df[,MONO_COLUMNS]

# 2) Generate design
non_const_columns <- MONO_COLUMNS[(sapply(df_data, var) != 0)]
dim <- length(non_const_columns) # ignore columns that have 0 variance
D_init <- MaxProLHD(args$num_meals*4, dim)$Design # Generate more than needed
D <- MaxPro(D_init)

# 3) Select best design
D_all <- MaxProAugment(df_data[,non_const_columns], D$Design, args$num_meals) # Select N new designs given the data
df_new_D <- data.frame(D_all$Design[(nrow(df_data)+1):nrow(D_all$Design),])
df_design <- data.frame(matrix(0, nrow(df_new_D), length(MONO_COLUMNS)))
colnames(df_design) <- MONO_COLUMNS
df_design[,non_const_columns] <- df_new_D

# 4) Find food proportions for each design (considering additive model)
df_food_proportions <- data.frame(matrix(0, nrow(df_design), nrow(df)))
colnames(df_food_proportions) <- df$uid
for(i in 1:nrow(df_design)){
  df_food_proportions[i,] <- get_proportions(df_data, df_design[i,])
}

# 5) Save
save_MMs(df_food_proportions, args$output)
  
