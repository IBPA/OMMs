# Meal Design for Heterogeneous Glycan Content
Design mixed meals automatically for maximizing the information content of their glycan compositions.

## Requirements
* **[Gurobi](https://www.gurobi.com/)** (version ≥ 9.1)
* **[R](https://www.r-project.org/)** (version ≥ 4.0.2)
* **R Packages** gurobi, MaxPro, ggplot2, readxl, stringr, argparse, infotheo, reshape2

## Installation
### Step1) Install Gurobi
Install [Gurobi](https://www.gurobi.com/) (version ≥ 9.1). The Gurobi academic licence is free. Then run the following command (ATTENTION: use the location of the Gurobi R package on your system which varies depending on your installation):<br>
`export GUROBI_R_PATH=$HOME/local/gurobi912/linux64/R/gurobi_9.1-2_R_4.0.2.tar.gz`

### Step2) Install OMMs using conda
Ensure that you have conda installed (see [miniconda](https://docs.conda.io/en/latest/miniconda.html)), then run:<br>
`conda create -n omms -c conda-forge -c ameenetemady r-omms`<br>
`conda activate omms`<br>


### Step3) Setup OMMs to use Gurobi
Run: `setup_OMMs.R  --gurobi-R-package $GUROBI_R_PATH`

## 1. OMM Generation
Design optimal mixed meals (OMMs) to maximize the information content of glycan profiles. In other words, in terms of the glycan content, we want the meals to be most different from each other while also being most different from the individual foods. To run, use the following command (change the argument values as needed):

`generate_OMMs.R --glycanDB glycanDB.csv --moistureDB moistureDB.xlsx --num-meals N [--max-meal-ingredients M] --output generated_OMMs.csv`
### Arguments:
* `--glycanDB filename.csv`: The glycan content profiles of individual foods when dried. The format is described in [./data/glycanDB_format.txt](./data/glycanDB_format.txt).
* `--moistureDB filename.xlsx`: The moisture percentages of individual foods before drying. The format is described in [./data/moistureDB_format.txt](./data/moistureDB_format.txt).
* `--max-meal-ingredients M` (optional): The maximum number of ingredients in each designed mixed meal.
* `--compositional` (optional): Transform glycan content of each food to compositional values. If not set, glyans will be transformed using minimax normalization for each glycan independently.
* `--num-meals N`: The number of mixed meals to design.
* `--output filename.csv`: The designed OMMs where each mixed meal is defined by the individual food proportions that it includes. The format is described in [./data/MMs_format.txt](./data/MMs_format.txt).


## 2. OMM Selection
Select a set of OMMs from candidate mixed meals, given a set of approved mixed meals. To run, use the following command (change the argument values as needed):

`select_OMMs.R --glycanDB glycanDB.csv --moistureDB moistureDB.xlsx [--approved-meals approved_meals.csv] --candidate-meals candidate_meals.csv --num-meals N --output selected_OMMs.csv`

### Arguments:
* `--glycanDB filename.csv`: The glycan content profiles of individual foods when dried. The format is described in [./data/glycanDB_format.txt](./data/glycanDB_format.txt).
* `--moistureDB filename.xlsx`: The moisture percentages of individual foods before drying. The format is described in [./data/moistureDB_format.txt](./data/moistureDB_format.txt).
* `--num-meals N`: The number of mixed meals to select.
* `--approved-meals filename.csv` (optional): The set of approved meals. The format is described in [./data/MMs_format.txt](./data/MMs_format.txt).
* `--compositional` (optional): Transform glycan content of each food to compositional values. If not set, glyans will be transformed using minimax normalization for each glycan independently.
* `--candidate-meals filename.csv`: The set of candidate meals to select from. The format is described in [./data/MMs_format.txt](./data/MMs_format.txt).
* `--output filename.csv`: The selected OMMs where each mixed meal is defined by the individual food proportions that it includes. The format is described in [./data/MMs_format.txt](./data/MMs_format.txt).

## 3. OMM Visualization
Visualize the expected glycan content of OMMs. *Figure1* will include a stacked barchart visualizing the glycan content of each OMM. *Figure2* will show the information content of the OMMs' expected glycan profiles when measured sequentially. To run, use the following command (change the argument values as needed):

`visualize_OMMs.R --glycanDB glycanDB.csv --moistureDB moistureDB.xlsx --OMMs OMMs.csv --output-dir ./resutls`

### Arguments:
* `--glycanDB filename.csv`: The glycan content profiles of individual foods when dried. The format is described in [./data/glycanDB_format.txt](./data/glycanDB_format.txt)
* `--moistureDB filename.xlsx`: The moisture percentages of individual foods before drying. The format is described in [./data/moistureDB_format.txt](./data/moistureDB_format.txt)
* `--OMMs filename.csv`: The OMMs to visualize. The format is described in [./data/MMs_format.txt](./data/MMs_format.txt).
* `--output-dir dir_path`: The output directory for saving the generated figures.

## Examples
See [examples.sh](./test/examples.sh). Following a successful installation, you can run `./test/examples.sh`.

## Support
For any questions contact Ameen Eetemadi (eetemadi@ucdavis.edu).

## Licence
See the [LICENSE](./LICENSE) file for license rights and limitations (MIT License).

## Acknowledgement
This work was supported by the United States Department of
Agriculture (USDA)/NSF AI Institute for Next Generation Food Systems (AIFS), USDA award number 2020-67021-32855.

