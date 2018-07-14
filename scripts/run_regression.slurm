#!/bin/bash

#SBATCH --mail-user=jakejhughey@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1

# test1
# #SBATCH --cpus-per-task=6
# #SBATCH --time=00:10:00
# #SBATCH --mem=16G

# test2
# #SBATCH --cpus-per-task=6
# #SBATCH --time=01:00:00
# #SBATCH --mem=32G

# test3
#SBATCH --cpus-per-task=6
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00

module restore R343

# Rscript scripts/run_regression.R params/exome_params_test1.yaml
# Rscript scripts/run_regression.R params/exome_params_test2.yaml
Rscript scripts/run_regression.R params/exome_params_test3.yaml
