#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
library(infotheo)
library(ggplot2)
library(reshape2)
library(scales)

source("./src/common.R")
source("./src/theme_util.R")

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
  args <- parser$parse_args(testargs)
  
  return(args)
}

args <- get_args()

# 1) Load data
MMs <- load_MMs(args$OMMs)
df_food_vectors <- load_data(args$glycanDB, args$moistureDB)
rownames(df_food_vectors) <- df_food_vectors$uid

# 2) Calculate the glycan vectors of mixed meals
MM_vectors <- MMs %*% data.matrix(df_food_vectors[colnames(MMs), MONO_COLUMNS])

# 3) Plot stacked bar-charts of mixed meals
df_MM <- minmax_normalize(data.frame(MM_vectors))
df_MM <- data.frame(t(apply(df_MM, 1, FUN = function(x) {x/sum(x)})))
df_MM$id <- 1:nrow(df_MM)
df_MM_melt <- melt(df_MM, measure.vars=MONO_COLUMNS, value.name = "Quantity", variable.name = "Glycan")
gPlot_bars <- ggplot(df_MM_melt, aes(x=id, y=Quantity, fill=Glycan))+
          geom_bar(stat = "identity")+
          scale_x_continuous("Meal ID")+
          scale_y_continuous("Glycan Proportions")+
          my_base_theme

print(gPlot_bars)
ggsave("./results/Fig1.pdf", gPlot_bars, dpi = 400, width=8 , height = 4)
ggsave("./results/Fig1.png", gPlot_bars, dpi = 400, width=8 , height = 4)


# 4) Calculate information content iteratively
df_vectors_all <- rbind(df_food_vectors[,MONO_COLUMNS], data.frame(MM_vectors))
df_vectors_all <-minmax_normalize(df_vectors_all)
df_vectors_all <- discretize(df_vectors_all, disc = "equalwidth")

infos <- c()
i_numbers <- 1:nrow(df_vectors_all)
for (i in i_numbers){
  df_vectors_cur <- df_vectors_all[1:i,]
  info_cur <- entropy(df_vectors_cur)
  print(info_cur)
  infos <- c(infos, info_cur)
}

# 5) Plot the information content of inidividual foods and mixed meals
df_infos <- data.frame(iter=i_numbers, info=infos)
gPlot <- ggplot(df_infos, aes(x=iter-nrow(df_food_vectors), y=infos))+
          geom_point(size=.7)+
          geom_line()+
          geom_text(x=nrow(MMs)/2, y=0, label=sprintf("%d mixed meals",nrow(MMs)))+
          geom_text(x=-nrow(df_food_vectors)/2, y=0, label=sprintf("%d individual foods",nrow(df_food_vectors)))+
          geom_vline(xintercept=0, linetype="dashed", color = "red")+
          scale_x_continuous("", breaks = c())+
          scale_y_continuous("Glycan information content (entropy)")+
          my_base_theme
print(gPlot)
ggsave("./results/Fig2.pdf", gPlot, dpi = 400, width=8 , height = 4)
ggsave("./results/Fig2.png", gPlot, dpi = 400, width=8 , height = 4)
