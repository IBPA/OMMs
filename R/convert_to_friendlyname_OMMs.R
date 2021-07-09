#!/usr/bin/env Rscript
#' Replace the uids in a mixed-meals file, with the corresponding friendly food names.

suppressPackageStartupMessages(library("argparse"))
library(stringr)

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

# source the filename from this script's directory
source_from_here <- function(filename) {
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  source(file.path(dirname(script.name), filename))
}
source_from_here("common_OMMs.R")

args <- get_args()

# 1) Load mixed-meals
MMs <- load_MMs(args$OMMs)

# 2) Replace the column names (uids) with corresponding friendly names
colnames(MMs) <- get_friendly_food_names(colnames(MMs), args$glycanDB)

# 3) Save
save_MMs(MMs, args$output_filename)
print("Completed Successfully!")
