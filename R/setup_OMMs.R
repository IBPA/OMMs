#!/usr/bin/env Rscript
#' One-time setup script

suppressPackageStartupMessages(library("argparse"))
library(stringr)

# Parse and return command line arguments of this script
get_args <- function() {
  parser <- ArgumentParser(description = "One-time setup for OMMs package.")
  parser$add_argument("--gurobi-R-package", required=TRUE,
                      help="The path to the gurobi R package.")
  args <- parser$parse_args()
  return(args)
}

args <- get_args()
install.packages(args$gurobi_R_package, repos=NULL)
