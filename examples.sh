#!/bin/bash

./generate_OMMs.R --glycanDB ./data/draft_mono_for_ameen_060321.csv \
                  --moistureDB "./data/% moisture content for all foods.xlsx" \
                  --num-meals 50 \
                  --max-meal-ingredients 5 \
                  --compositional \
                  --output ./results/gen_OMMS.csv

./visualize_OMMs.R --glycanDB ./data/draft_mono_for_ameen_060321.csv \
                   --moistureDB "./data/% moisture content for all foods.xlsx" \
                   --OMMs ./results/gen_OMMS.csv \
                   --output-dir ./results/

./select_OMMs.R --glycanDB ./data/draft_mono_for_ameen_060321.csv \
                --moistureDB "./data/% moisture content for all foods.xlsx" \
                --candidate-meals ./results/cand_OMMS.csv \
                --approved-meals ./results/appr_OMMS.csv \
                --num-meals 10 \
                --compositional \
                --output ./results/sel_OMMS.csv

./convert_to_friendlyname_OMMs.R --glycanDB ./data/draft_mono_for_ameen_060321.csv \
                                 --OMMs ./results/gen_OMMS.csv \
                                 --output-filename ./results/gen_OMMS_friendly.csv