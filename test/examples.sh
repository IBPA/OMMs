#!/bin/bash

# Change to this script's directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
pushd $SCRIPT_DIR

mkdir results

# Generate 50 OMMs using 5 ingredients
generate_OMMs.R --glycanDB ./data/draft_mono_for_ameen_060321.csv \
                  --moistureDB "./data/% moisture content for all foods.xlsx" \
                  --num-meals 50 \
                  --max-meal-ingredients 5 \
                  --compositional \
                  --output ./results/gen_OMMS.csv

# Visualize generated OMMs
visualize_OMMs.R --glycanDB ./data/draft_mono_for_ameen_060321.csv \
                  --moistureDB "./data/% moisture content for all foods.xlsx" \
                  --OMMs ./results/gen_OMMS.csv \
                  --output-dir ./results/


# Convert generated OMMs to use friendly food neames
convert_to_friendlyname_OMMs.R --glycanDB ./data/draft_mono_for_ameen_060321.csv \
                               --OMMs ./results/gen_OMMS.csv \
                               --output-filename ./results/gen_OMMS_friendly.csv

# Select 10 OMMs from the generated OMMs
select_OMMs.R --glycanDB ./data/draft_mono_for_ameen_060321.csv \
              --moistureDB "./data/% moisture content for all foods.xlsx" \
              --candidate-meals ./results/gen_OMMS.csv \
              --num-meals 10 \
              --compositional \
              --output ./results/sel_OMMS.csv


# Back to original directory
popd