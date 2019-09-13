SIMLogistLargeNK <- function(n, p, q, sim, njob) {
  # seed to generate random number
  # n is the total sample size
  # p is the number of true risk alles
  # q is the number of false risk alles
  # compare logistic regression with SNP and rcs(last.age, 5) and Cox model with right censoring
  # the right model is logistic regression model
  
  # Genv <- matrix(NA, n, p+q)
  GenvP <- matrix(NA, n, p)
  GenvQ <- rep(0,n)
  pval <- matrix(9.99, p+q, 2) # create double matrix with prefilled impossible value
  nominal <- c(0.01, 0.05, 0.1)
  L <- length(nominal)
  LogitTPR <- LogitTNR <- CoxTPR <- CoxTNR <- rep(NA, L)
  
  cat("starting");  print(Sys.time())
  
  # generate the data
  # generate risk alleles
  prob <- as.vector(quantile(minor_allele_fraction, runif(p+q)))
  for(k in 1:p) GenvP[,k] <- rbinom(n, size=1, prob=prob[k])
  # generate the age of last record
  age <- rnorm(n, 60, 5)
  # the first p is true risk alles with effect randomly simulated from Unif(0.5,1) 
  # the rest q alles are false risk alles
  beta <- runif(p, 0.3,0.7)
  beta0 <- -11.5
  beta1 <- 0.001
  # event indicator 
  EventP <- 1/(1+exp(-beta0-GenvP[,1:p]%*%beta-beta1*age))
  delta <- rbinom(n, size=1, prob=EventP)
  
  # get the censoring and observed time
  censt <- runif(n,50,85)
  survt <- ifelse(delta==1, age, censt)
    
  cat("finish data simulation"); print(Sys.time())
  
  eventr <- mean(delta)
  cat("event rate", mean(delta), "\n")
  
  # Fit the Cox model without truncation, and logistic regression model
  # use survt, delta, and Genv from now on
  # for k in 1:p use simulated GenvP, otherwise simulate new null SNPs
  for(k in 1:p) {
    gc()
    sf <- tryCatch(coxph.fit(x=as.matrix(as.numeric(GenvP[,k])), y=Surv(survt, delta),strata=NULL, init=NULL, method='breslow', 
            control=coxph.control(), rownames=c('x1')), error = function(e) NA)
    sl <- tryCatch(glm(delta ~ GenvP[,k] + rcs(survt, 5), family=binomial(link = "logit")), error = function(e) NA)
    if(is.list(sf) & is.list(sl)) {
      pval[k,1] <- coef(summary(sl))[2,4]
      pval[k,2] <- 2*pnorm(-abs(sf$coefficients[1]/sqrt(sf$var[1,1])))
    }
  }  
  cat("Finish True SNPs model"); print(Sys.time())
  rm(list=c('GenvP'))
  
  file.out <- paste0("Pval_",sim, "_", njob,".Rdata")
  save(pval, file=file.out)
  
  for(k in (p+1):(p+q)) {
    gc()
    GenvQ <- rbinom(n, size=1, prob=prob[k])
    sf <- tryCatch(coxph.fit(x=as.matrix(as.numeric(GenvQ)), y=Surv(survt, delta),strata=NULL, init=NULL, method='breslow', 
                             control=coxph.control(), rownames=c('x1')), error = function(e) NA)
    sl <- tryCatch(glm(delta ~ GenvQ + rcs(survt, 5), family=binomial(link = "logit")), error = function(e) NA)
    if(is.list(sf) & is.list(sl)) {
      pval[k,2] <- 2*pnorm(-abs(sf$coefficients[1]/sqrt(sf$var[1,1])))
      pval[k,1] <- coef(summary(sl))[2,4]
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
    CoxTPR[l] <- mean(pval[1:p,2]<nominal[l]/(p+q), na.rm=T) 
    CoxTNR[l] <- mean(pval[(p+1):(p+q),2]>nominal[l]/(p+q), na.rm=T) 
  }
  
  rm(list=ls()[!ls() %in% c("eventr", "LogitTPR", "LogitTNR", "CoxTPR", "CoxTNR")])
  
  return(c(eventr, LogitTPR, LogitTNR, CoxTPR, CoxTNR))
}

