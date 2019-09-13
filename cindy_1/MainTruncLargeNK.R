# rm(list = ls())

args = commandArgs(TRUE)
njob = as.integer(args[1])
wd = args[2]
# n =args[3]
# K = args[4]
# prate = args[5]

library(survival)
library(rms)
# load("./Data/summary_demographics.Rdata")
# load("./Data/summary_geno_pheno.Rdata")
# source("SIMTrunc.R")
load("~/JakeHughey/summary_demographics.Rdata")
load("~/JakeHughey/summary_geno_pheno.Rdata")
source("~/JakeHughey/SIMTruncLargeNK2.R")

# ResultCox
 NSIM <- 1
 n <- 50000; p <- 100; q <- 800000-p

# ResultCox0
# NSIM <- 10
# n <- 50000; p <- 100; q <- 101-p

set.seed(12345+5000*njob)
file.out <- paste0("Result", njob, ".txt")
setwd(wd)

Result <- matrix(NA, NSIM, 14)
for(sim in 1:NSIM) {
   print(Sys.time())
   cat("sim=", sim, "\t")
   Result[sim, ] <- SIMTruncLargeNK2(n, p, q, sim, njob)
   write.table(Result, file=file.out, row.names=FALSE, col.names=FALSE)
 }
write.table(Result, file=file.out, row.names=FALSE, col.names=FALSE)

# apply(Result, 2, mean)
