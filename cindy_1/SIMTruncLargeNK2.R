SIMTruncLargeNK2 <- function(n, p, q, sim, njob) {
  # seed to generate random number
  # n is the total sample size
  # p is the number of true risk alles
  # q is the number of false risk alles
  # compare logistic regression with SNP and rcs(last.age, 5) and Cox model with right censoring and left truncation
  
  # Genv <- matrix(NA, n, p+q)
  GenvP <- matrix(NA, n, p)
  GenvQ <- rep(0,n)
  pval <- matrix(9.99, p+q, 2) # create double matrix with prefilled impossible value
  nominal <- c(0.01, 0.05, 0.1)
  L <- length(nominal)
  # CoxTPR <- CoxTNR <- 
  LogitTPR <- LogitTNR <- TCoxTPR <- TCoxTNR <- rep(NA, L)
 
  cat("starting");  print(Sys.time())
  
  # generate the data
  # generate risk alleles
  prob <- as.vector(quantile(minor_allele_fraction, runif(p+q)))
  # for(k in 1:(p+q)) Genv[,k] <- rbinom(n, size=1, prob=prob[k])
  for(k in 1:p) GenvP[,k] <- rbinom(n, size=1, prob=prob[k])
  # the first p is true risk alles with effect randomly simulated from Unif(0.5,1) 
  # the rest q alles are false risk alles
  beta <- runif(p, 0.3,0.5)
  # beta <- rep(0.5,p)
  # generate the survival time S(t)=exp(-Lamda0(t)*exp(beta'Genv))~Unif(0,1) 
  # Lambda0(t)=10000t below to get reasonable event rate
  survt <- -log(runif(n))*exp(-GenvP[,1:p]%*%beta)*10000
  
  # get the censoring
  censt <- pmin(rgamma(n,1,1), 2)
  delta <- ifelse(survt <= censt, 1, 0)
  survt <- pmin(survt, censt)
  
  # simulate truncation time and drop off the rows with survt<=trunct
  trunct <- runif(n, 0, 0.1)
  keep <- (trunct<survt)
  trunct <- trunct[keep]
  survt <- survt[keep]
  delta <- delta[keep]
  GenvP <- GenvP[keep,]
  diffage <- survt-trunct
  cat("finish data simulation"); print(Sys.time())
  
  eventr <- mean(delta)
  n1 <- sum(keep)
  cat("n1=", n1, "event rate", mean(delta), "\n")

  # Fit the Cox model ignoring left truncation, Cox model with left truncation, and logistic regression model
  # use survt, delta, and Genv from now on
  # for k in 1:p use simulated GenvP, otherwise simulate new null SNPs
  for(k in 1:p) {
    gc()
    sft <- tryCatch(agreg.fit(x=as.matrix(as.numeric(GenvP[,k])), y=Surv(trunct, survt, delta),strata=NULL, init=NULL, method='breslow', 
           control=coxph.control(), rownames=c('x1')), error = function(e) NA)
    sl <- tryCatch(glm(delta ~ GenvP[,k] + rcs(survt, 5) + rcs(diffage, 5), family=binomial(link = "logit")), error = function(e) NA)
    if(is.list(sft) & is.list(sl)) {
      pval[k,1] <- coef(summary(sl))[2,4]
      pval[k,2] <- 2*pnorm(-abs(sft$coefficients[1]/sqrt(sft$var[1,1])))
    }
  }  
  cat("Finish True SNPs model"); print(Sys.time())
  rm(list=c('censt', 'GenvP', 'keep'))

  file.out <- paste0("Pval_",sim, "_", njob,".Rdata")
  save(pval, file=file.out)
  
  for(k in (p+1):(p+q)) {
    gc()
    GenvQ <- rbinom(n1, size=1, prob=prob[k])
    sft <- tryCatch(agreg.fit(x=as.matrix(as.numeric(GenvP[,k])), y=Surv(trunct, survt, delta),strata=NULL, init=NULL, method='breslow', 
           control=coxph.control(), rownames=c('x1')), error = function(e) NA)
    sl <- tryCatch(glm(delta ~ GenvQ + rcs(survt, 5) + rcs(diffage, 5), family=binomial(link = "logit")), error = function(e) NA)
    if(is.list(sft) & is.list(sl)) {
      pval[k,1] <- coef(summary(sl))[2,4]
      pval[k,2] <- 2*pnorm(-abs(sft$coefficients[1]/sqrt(sft$var[1,1])))
    }
    if(k %% 1000 == 0) {
     save(pval, file=file.out)
     cat("SNP=",k,"\n"); 
     print(Sys.time())
    }
  }  
  cat("Finish Null effect SNPs model"); print(Sys.time())
  save(pval, file=file.out)
  
  # calculate TPR and TNR based on Bonforroni correction
  # pval has three columns, 1st one from cox without truncation, 2nd one from Logistic regression model, 3rd one is cox with truncation 
  for (l in 1:L) {
    LogitTPR[l] <- mean(pval[1:p,1]<nominal[l]/(p+q), na.rm=T) 
    LogitTNR[l] <- mean(pval[(p+1):(p+q),1]>nominal[l]/(p+q), na.rm=T) 
    TCoxTPR[l] <- mean(pval[1:p,2]<nominal[l]/(p+q), na.rm=T) 
    TCoxTNR[l] <- mean(pval[(p+1):(p+q),2]>nominal[l]/(p+q), na.rm=T) 
  }

  rm(list=ls()[!ls() %in% c("n1", "eventr", "LogitTPR", "LogitTNR", "TCoxTPR", "TCoxTNR")])

  return(c(n1, eventr, LogitTPR, LogitTNR, TCoxTPR, TCoxTNR))
}
