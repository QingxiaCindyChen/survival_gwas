#!/bin/bash
#SBATCH --mail-user=jakejhughey@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
module restore R343
Rscript scripts/gather_regression.R results/exome/20180715_204406_
