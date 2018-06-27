#ROC plot/calculation setup

#Calculate TPR and FPR from pvalues and binary GWAS catalogue calls
ROCLobster = function(SurvData, threshold){
  TP = nrow(SurvData[SurvData$GWAS==1 & SurvData$pval < threshold, ]) #snps in both significant lists below our pval
  FP = nrow(SurvData[SurvData$GWAS==0 & SurvData$pval < threshold, ]) #snps not found in the gwas catalogue but below our pval
  TN = nrow(SurvData[SurvData$GWAS==0 & SurvData$pval > threshold, ]) #snps not found in the gwas catalogue but above pval
  FN = nrow(SurvData[SurvData$GWAS==1 & SurvData$pval > threshold, ])
  TPR = TP/(TP+FN) #y axis
  FPR = FP/(FP+TN) #x axis
  Out = data.frame(TPR = TPR, FPR = FPR)
  return(Out)
}

#ROC Plot
ROCPlobster = function(PlotDF, PlotTitle){
  ROC = ggplot(PlotDF, aes(x = FPR, y = TPR, group = test, colour = test))+
    scale_x_continuous(expand = c(0.0075, 0.0075)) + 
    scale_y_continuous(expand = c(0.005, 0.005)) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), lty = 2, color = 'black')+
    geom_line(size = 1)+
    theme_bw()+
    ggtitle(paste(PlotTitle))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11, face = "bold"), 
          #panel.grid.major = element_line(colour = 'grey98'), panel.grid.minor = element_line(colour = 'grey98'),
          legend.title =  element_blank(), legend.text = element_text(size = 10), 
          legend.box.margin = margin(-2, -2, -2, -2), title = element_text(size = 11, face = "bold"), 
          plot.title = element_text(hjust = 0.5))
  return(ROC)
}

PhePhetchROC = function(FileNames, CallMatrix){ #Build ROC data for a given list of phecode file names
  
  ROCMat = NULL
  for(code in FileNames){
    CoxDF = setDF(read_tsv(file = file.path(filepath, 'GWASCatalog', 'exome_full', sprintf('exome_%s_cox.tsv.gz', code)), col_types = cols()))
    LogMatch = Logfiles[c(grep(sprintf('exome_%s_logistic.tsv.gz', code), Logfiles))]
    LogDF = setDF(read_tsv(file = file.path(filepath, 'GWASCatalog', 'exome_full', LogMatch), col_types = cols()))
    LogDF = rename(LogDF, pval = P, snp = SNP)
    
    #Get 'raw' phecode ID
    Phecode = gsub('phe', '', code)
    Phecode = gsub('p', '.', Phecode)
    
    #Only take what you need
    CoxDF = select(CoxDF, snp, pval)
    CoxDF$phecode = Phecode
    CoxDF$method = 'Cox'
    CoxDF = merge(CoxDF, SnpIDs, by.x = 'snp', by.y = 'SNP')
    CoxDF = select(CoxDF, -snp)
    CoxDF = CoxDF[CoxDF$rsID %in% CallMatrix$SNP, ]
    temp1 = inner_join(CoxDF, CallMatrix, by = c('rsID' = 'SNP', 'phecode' = 'jd_code'))
    if(nrow(temp1)>0){
      temp1$GWAS = 1
    }
    temp2 = anti_join(CoxDF, CallMatrix, by = c('rsID' = 'SNP', 'phecode' = 'jd_code'))
    if(nrow(temp2)>0){
      temp2$GWAS = 0
    }
    
    CoxDF = rbind(temp1, temp2)
    
    LogDF = select(LogDF, snp, pval)
    LogDF$phecode = Phecode
    LogDF$method = 'Logistic'
    LogDF = merge(LogDF, SnpIDs, by.x = 'snp', by.y = 'SNP')
    LogDF = select(LogDF, -snp)
    LogDF = LogDF[LogDF$rsID %in% CallMatrix$SNP, ]
    temp1 = inner_join(LogDF, CallMatrix, by = c('rsID' = 'SNP', 'phecode' = 'jd_code'))
    if(nrow(temp1)>0){
      temp1$GWAS = 1
    }
    temp2 = anti_join(LogDF, CallMatrix, by = c('rsID' = 'SNP', 'phecode' = 'jd_code'))
    if(nrow(temp2)>0){
      temp2$GWAS = 0
    }
    
    LogDF = rbind(temp1, temp2)
    
    #ID match
    ComboDF = rbind(CoxDF, LogDF)
    
    ROCtemp = select(ComboDF, pval, method, GWAS)
    ROCMat = rbind(ROCMat, ROCtemp)
  }
  
  ROCMat$pval = as.numeric(ROCMat$pval)
  return(ROCMat)
}

#Calculate AUC - found at https://mbq.me/blog/augh-roc/
auroc <- function(score, bool) {
  n1 <- sum(!bool)
  n2 <- sum(bool)
  U  <- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

BuildROCCurve = function(ROCDF, PlotTitle){ #Make the roc curve
  
  #Loop through thresholds yourself, and consider the very low pvalues (down to 1e-20)
  thresholds = c(1:100 %o% 10^-(30:2)) #fancy sequence, outer binary product operator - 2900 values here
  
  fillercox = as.data.frame(matrix(NA, nrow = length(thresholds), ncol = 2))
  colnames(fillercox) = c('TPR', 'FPR')
  fillerlog = as.data.frame(matrix(NA, nrow = length(thresholds), ncol = 2))
  colnames(fillerlog) = c('TPR', 'FPR')
  
  SurvDataCox = ROCDF[ROCDF$method=='Cox', ]
  #NAs in Logistic
  SurvDataLog = ROCDF[ROCDF$method=='Logistic', ]
  SurvDataLog = SurvDataLog[!is.na(SurvDataLog$pval), ]
  
  #the number of thresholds specified above, for 80+ phecodes, can take a few minutes
  for(thresh in 1:length(thresholds)){
    tempcox = ROCLobster(SurvDataCox, thresholds[thresh])
    templog = ROCLobster(SurvDataLog, thresholds[thresh])
    fillercox[thresh, ] = tempcox
    fillerlog[thresh, ] = templog
  }
  
  fillercox$test = c('CoxPH')
  fillerlog$test = c('Logistic')
  filler = rbind(fillercox, fillerlog)
  
  #Build plot
  Output = ROCPlobster(filler, PlotTitle)
  
  #Add AUC calc
  #Flip pval to threshold
  FlipCox = 1 - SurvDataCox$pval
  CoxAUC = auroc(FlipCox, SurvDataCox$GWAS)
  
  FlipLog = 1 - SurvDataLog$pval
  LogAUC = auroc(FlipLog, SurvDataLog$GWAS)
  
  OutROC = Output + annotate("text", x = 0.79, y = 0.1, label = paste("Cox AUC = ", round(CoxAUC, 3)), size = 3)+
    annotate("text", x = 0.763, y = 0.05, label = paste("Logistic AUC = ", round(LogAUC, 3)), size = 3)
  
  return(OutROC)
}
