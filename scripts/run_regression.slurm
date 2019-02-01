#!/bin/bash
#SBATCH --mail-user=jakejhughey@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=01:00:00
module restore R343
Rscript scripts/run_regression.R params/exome/params_test2.yaml
