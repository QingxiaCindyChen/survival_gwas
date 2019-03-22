
rm(list = ls())
library(survival)
load("./Data/summary_demographics.Rdata")
load("./Data/summary_geno_pheno.Rdata")

source("SIMTrunc.R")

NSIM <- 1000
set.seed(123)
n <- 29000; p <- 100; q <- 30000-p

Result <- matrix(NA, NSIM, 8)
for(sim in 1:NSIM) {
   print(Sys.time())
   cat("sim=", sim, "\t")
   Result[sim, ] <- SIMTrunc(n, p, q)
 }
write.csv(Result, file="Result.csv")

# apply(Result, 2, mean)
