#!/usr/bin/env Rscript
#' Visualize the expected glycan content of OMMs.

suppressPackageStartupMessages(library("argparse"))
library(infotheo)
library(readxl)
library(stringr)
library(ggplot2)
library(reshape2)
library(scales)

DEV_MODE <- FALSE

# Parse and return command line arguments of this script (return defaults if DEV_MODE is TRUE).
get_args <- function() {
  parser <- ArgumentParser(description = "Visualize the expected glycan content of OMMs.")
  parser$add_argument("--glycanDB", required=TRUE,
                      help="A .csv file containing glycan content profiles of individual foods when dried.")
  parser$add_argument("--moistureDB", required=TRUE,
                      help="A .xlsx file containing the moisture percentages of individual foods before drying.")
  parser$add_argument("--OMMs", required=TRUE,
                      help="A .csv file containing the the OMMs to visualize.")
  parser$add_argument("--output-dir", required=TRUE,
                      help="The output directory for saving the generated figures.")
  
  testargs <- c("--glycanDB", "./data/draft_mono_for_ameen_060321.csv", 
                "--moistureDB", "./data/% moisture content for all foods.xlsx", 
                "--OMMs", "./results/gen_OMMS.csv", 
                "--output-dir", "./results/")
  if (DEV_MODE){
    return(parser$parse_args(testargs))
  }
  return(parser$parse_args())
}

# source the filename from this script's directory
source_from_here <- function(filename) {
  if (DEV_MODE){
    source(file.path("R", filename))
  } else {
    initial.options <- commandArgs(trailingOnly = FALSE)
    script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
    source(file.path(dirname(script.name), filename))
  }
}
source_from_here("common_OMMs.R")
source_from_here("theme_util_OMMs.R")

args <- get_args()

# 1) Load data
MMs <- load_MMs(args$OMMs)
df_food_vectors <- load_food_data(args$glycanDB, args$moistureDB)
rownames(df_food_vectors) <- df_food_vectors$uid

# 2) Calculate the glycan vectors of mixed meals
MM_vectors <- MMs %*% data.matrix(df_food_vectors[colnames(MMs), MONO_COLUMNS])

# 3.A) Plot stacked bar-charts of mixed meals
df_MM <- data.frame(MM_vectors)
df_MM$id <- 1:nrow(df_MM)
df_MM_melt <- melt(df_MM, measure.vars=MONO_COLUMNS, value.name = "Quantity", variable.name = "Glycan")
gPlot_bars <- ggplot(df_MM_melt, aes(x=id, y=Quantity, fill=Glycan))+
          geom_bar(stat = "identity")+
          scale_x_continuous("Meal ids")+
          scale_y_continuous("Glycan weight per mg of wet sample")+
          my_base_theme

print(gPlot_bars)
ggsave(file.path(args$output_dir,"Fig1_A.pdf"), gPlot_bars, dpi = 400, width=8 , height = 4)
ggsave(file.path(args$output_dir,"Fig1_A.png"), gPlot_bars, dpi = 400, width=8 , height = 4)

# 3.B) Plot stacked bar-charts of mixed meals (proportions)
df_MM <- data.frame(MM_vectors)
df_MM$id <- 1:nrow(df_MM)
df_MM_melt <- melt(df_MM, measure.vars=MONO_COLUMNS, value.name = "Quantity", variable.name = "Glycan")
gPlot_bars <- ggplot(df_MM_melt, aes(x=id, y=Quantity, fill=Glycan))+
  geom_bar(stat = "identity", position = "fill")+
  scale_x_continuous("Meal ids")+
  scale_y_continuous("Glycan %weight per mg of wet sample", labels = scales::percent)+
  my_base_theme
print(gPlot_bars)
ggsave(file.path(args$output_dir,"Fig1_B.pdf"), gPlot_bars, dpi = 400, width=8 , height = 4)
ggsave(file.path(args$output_dir,"Fig1_B.png"), gPlot_bars, dpi = 400, width=8 , height = 4)

# 3.C) Plot stacked bar-charts of foods (proportions)
df_food_vectors_vis <- data.frame(df_food_vectors)
df_food_vectors_vis[is.na(df_food_vectors_vis)] <- 0
df_food_vectors_vis[df_food_vectors_vis < 0] <- 0
df_food_vectors_vis$id <- 1:nrow(df_food_vectors_vis)
df_food_vectors_melt <- melt(df_food_vectors_vis, measure.vars=MONO_COLUMNS, value.name = "Quantity", variable.name = "Glycan")
gPlot_bars <- ggplot(df_food_vectors_melt, aes(x=id, y=Quantity, fill=Glycan))+
  geom_bar(stat = "identity", position = "fill")+
  scale_x_continuous("Food Index")+
  scale_y_continuous("Glycan %weight per mg of wet sample", labels = scales::percent)+
  my_base_theme
print(gPlot_bars)
ggsave(file.path(args$output_dir,"Fig1_C.pdf"), gPlot_bars, dpi = 400, width=8 , height = 4)
ggsave(file.path(args$output_dir,"Fig1_C.png"), gPlot_bars, dpi = 400, width=8 , height = 4)

# 4) Calculate information content iteratively
df_vectors_all <- data.frame(MM_vectors)
df_vectors_all <- discretize(df_vectors_all, disc = "equalwidth", nbins = 20)

infos <- c()
i_numbers <- 1:nrow(df_vectors_all)
for (i in i_numbers){
  df_vectors_cur <- df_vectors_all[1:i,]
  info_cur <- entropy(df_vectors_cur)
  infos <- c(infos, info_cur)
}

# 5) Plot the information content of inidividual foods and mixed meals
df_infos <- data.frame(iter=i_numbers, info=infos)
gPlot <- ggplot(df_infos, aes(x=iter, y=infos))+
          geom_point(size=.7)+
          geom_line()+
          scale_x_continuous("Meal ids")+
          scale_y_continuous("Glycan information content (entropy)", breaks = 0:round(max(df_infos$info)))+
          my_base_theme
print(gPlot)
ggsave(file.path(args$output_dir,"Fig2.pdf"), gPlot, dpi = 400, width=8 , height = 4)
ggsave(file.path(args$output_dir,"Fig2.png"), gPlot, dpi = 400, width=8 , height = 4)

# 6) Histogram for number of foods used in the mixed-meals
df_hist <- data.frame(id=1:nrow(MMs), count = rowSums(MMs!=0))
gPlot_hist <- ggplot(df_hist, aes(count))+
                geom_histogram( colour="black", fill=unname(colors_assigned["A"]), binwidth = 1)+
                scale_x_continuous("Number of foods in the mixed meal", limits = c(1, max(df_hist$count)+1))+
                scale_y_continuous("Number of mixed meals")+
                my_base_theme
sprintf("Mixed meals use an average %.2f foods.", mean(df_hist$count))
print(gPlot_hist)
ggsave(file.path(args$output_dir,"Fig3.pdf"), gPlot_hist, dpi = 400, width=8 , height = 4)
ggsave(file.path(args$output_dir,"Fig3.png"), gPlot_hist, dpi = 400, width=8 , height = 4)
print("Completed Successfully!")
