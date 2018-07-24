#!/bin/bash
#SBATCH --mail-user=jakejhughey@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
module restore R343
Rscript scripts/split_regression.R params/exome/params_test4.yaml
