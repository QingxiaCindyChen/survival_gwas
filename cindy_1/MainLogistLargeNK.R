# rm(list = ls())

args = commandArgs(TRUE)
njob = as.integer(args[1])
wd = args[2]

library(survival)
library(rms)
# load("./Data/summary_demographics.Rdata")
# load("./Data/summary_geno_pheno.Rdata")
# source("SIMTrunc.R")
load("~/JakeHughey/summary_demographics.Rdata")
load("~/JakeHughey/summary_geno_pheno.Rdata")
source("~/JakeHughey/SIMLogistLargeNK.R")

NSIM <- 2
set.seed(12345+5000*njob)
file.out <- paste0("Result", njob, ".Rdata")
setwd(wd)

n <- 50000; p <- 100; q <- 800000-p

Result <- matrix(NA, NSIM, 13)
for(sim in 1:NSIM) {
   print(Sys.time())
   cat("sim=", sim, "\t")
   Result[sim, ] <- SIMLogistLargeNK(n, p, q, sim, njob)
   save(Result, file=file.out)
 }

save(Result, file=file.out)
