SIMTrunc <- function(n, p, q) {
  # seed to generate random number
  # n is the total sample size
  # p is the number of true risk alles
  # q is the number of false risk alles
  
  Genv <- matrix(NA, n, p+q)
  pval <- matrix(NA, p+q, 6)
  
  cat("1");  print(Sys.time())
  
  # generate the data
  # generate risk alleles
  prob <- as.vector(quantile(minor_allele_fraction, runif(p+q)))
  for(k in 1:(p+q)) Genv[,k] <- rbinom(n, size=1, prob=prob[k])
  # the first p is true risk alles with effect randomly simulated from Unif(0.5,1) 
  # the rest q alles are false risk alles
  beta <- runif(p, 0.3,0.5)
  # beta <- rep(0.5,p)
  # generate the survival time S(t)=exp(-Lamda0(t)*exp(beta'Genv))~Unif(0,1) 
  # Lambda0(t)=10000t below to get reasonable event rate
  survt <- -log(runif(n))*exp(-Genv[,1:p]%*%beta)*10000
  
  # get the censoring
  censt <- pmin(rgamma(n,1,1), 2)
  delta <- ifelse(survt <= censt, 1, 0)
  survt <- pmin(survt, censt)
  
  cat("2"); print(Sys.time())
  
  # simulate truncation time and drop off the rows with survt<=trunct
  trunct <- runif(n, 0, 0.1)
  keep <- (trunct<survt)
  trunct <- trunct[keep]
  survt <- survt[keep]
  delta <- delta[keep]
  Genv <- Genv[keep,]
  
  ## test on time
  cat("3"); print(Sys.time())
  system.time(agreg.fit(x=Genv[,k], y=cbind(trunct, survt, delta),strata=NULL, init=NULL, method='breslow', 
                        control=coxph.control(), rownames=c('x1')))
  testfun <- function(survt, trunct) {  
    sort.end <- order(-survt) - 1L
    sort.start <- order(-trunct) - 1L
  }
    system.time(testfun(survt, trunct))
  
  cat("4"); print(Sys.time())
  ## sorting time 
  sort.end <- order(-survt) - 1L
  sort.start <- order(-trunct) - 1L
  cat("5"); print(Sys.time())
  
  # Fit the Cox model ignoring left truncation, Cox model with left truncation, and logistic regression model
  # use survt, delta, and Genv from now on
  for(k in 1:(p+q)) {
    sf <- tryCatch(coxph(Surv(survt, delta) ~ Genv[,k]), error = function(e) NA)
    sft <- tryCatch(coxph(Surv(trunct, survt, delta) ~ Genv[,k]), error = function(e) NA)
    sl <- tryCatch(glm(delta ~ Genv[,k], family=binomial(link = "logit")), error = function(e) NA)
    if(is.list(sf) & is.list(sft) & is.list(sl)) {
      pval[k,1] <- 2*pnorm(-abs(sf$coefficients[1]/sqrt(sf$var[1,1])))
      pval[k,2] <- coef(summary(sl))[2,4]
      pval[k,3] <- 2*pnorm(-abs(sft$coefficients[1]/sqrt(sf$var[1,1])))
    }
  }  
  cat("4"); print(Sys.time())
  
  pval[,4] <- pval[,1]<(0.05/(p+q))
  pval[,5] <- pval[,2]<0.05/(p+q)
  pval[,6] <- pval[,3]<0.05/(p+q)
  
  cat("5"); print(Sys.time())
  
  # calculate TPR and TNR based on Bonforroni correction
  # pval has two columns, 1st one from cox model and 2nd one from Logistic regression model 
  CoxTPR <- mean(pval[1:p,1]<0.05/(p+q), na.rm=T) 
  CoxTNR <- mean(pval[(p+1):(p+q),1]>0.05/(p+q), na.rm=T) 
  LogitTPR <- mean(pval[1:p,2]<0.05/(p+q), na.rm=T) 
  LogitTNR <- mean(pval[(p+1):(p+q),2]>0.05/(p+q), na.rm=T) 
  TCoxTPR <- mean(pval[1:p,3]<0.05/(p+q), na.rm=T) 
  TCoxTNR <- mean(pval[(p+1):(p+q),3]>0.05/(p+q), na.rm=T) 
  eventr <- mean(delta)
  n1 <- sum(keep)
  cat("n1=", n1, "event rate", mean(delta), "\n")
  
  cat("6");  print(Sys.time())
  
  rm(list=ls()[!ls() %in% c("n1", "eventr", "CoxTPR", "CoxTNR", "LogitTPR", "LogitTNR", "TCoxTPR", "TCoxTNR")])

  return(c(n1, eventr, CoxTPR, CoxTNR, LogitTPR, LogitTNR, TCoxTPR, TCoxTNR))
}
