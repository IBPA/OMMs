#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
library(stringr)
source("./src/common.R")

get_args <- function() {
  parser <- ArgumentParser(description = "Use friendly food names instead of uids for the OMMs.")
  parser$add_argument("--glycanDB", required=TRUE,
                      help="A .csv file containing glycan content profiles of individual foods when dried.")
  parser$add_argument("--OMMs", required=TRUE,
                      help="A .csv file containing the the OMMs to visualize.")
  parser$add_argument("--output-filename", required=TRUE,
                      help="The name of the .csv file with friendly names.")

  args <- parser$parse_args()
  
  return(args)
}

args <- get_args()

MMs <- load_MMs(args$OMMs)
colnames(MMs) <- get_friendly_food_names(colnames(MMs), args$glycanDB)
save_MMs(MMs, args$output_filename)
