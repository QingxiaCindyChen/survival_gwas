#Updated GWAS catalog v1.0.2 - includes a snomed icd match too
#Inputs - GWASCatalog_ICD9_OMIM_SNOMEDCT_Map.csv, phecode_icd9_rolled.csv, ExomeIDs.csv, snomed_ICD9, gwas_catalog-ancestry.tsv
#Outputs - MultiPheROC_ICD_GWASv102.pdf, MultiPheROC_SNO_GWASv102.pdf, MultiPheROC_ICD_SNO_GWASv102.pdf
library('readr')
library('cowplot')
library('data.table')
library('survival')
library('snpStats')
library('survminer')
library('dplyr')
library('plotROC')
options(warn = -1)

filepath = file.path('/Users', 'srhoades', 'Dropbox (VUMC)/', 'JakeColab_CountingGWAS', 'counting_gwas', 'analyze')

#New with omim and snomedct
GWASCat = as.data.frame(read_csv(file = file.path(filepath, 'GWASCatalog', 'GWASCatalog_ICD9_OMIM_SNOMEDCT_Map.csv')))

#Exome
SnpIDs = setDT(read_csv(file = file.path(filepath, 'GWASCatalog', 'ExomeIDs.csv'), col_types = cols()))

#SNOMED
snoIDs = setDT(read_csv(file = file.path(filepath, 'GWASCatalog', 'snomed_ICD9'), col_types = 'cc'))

#Ancestry filter
Anc = as.data.frame(read_tsv(file = file.path(filepath, 'GWASCatalog', 'gwas_catalog-ancestry.tsv')))
#Messy names, just anything with European in it? many are mixed
AncList = as.character(Anc$PUBMEDID[c(grep('Euro', Anc$`BROAD ANCESTRAL CATEGORY`))])

phemap = read_csv(file = file.path(filepath, 'GWASCatalog', 'phecode_icd9_rolled.csv'))

#Trim rows without an ICD9
GWASCatICD = GWASCat[!(is.na(GWASCat$ICD)), ] #roughly half of the rows removed, but only 180 unique ICD9s here

#Keep NAs - since we want the rows which have SNOMED but not ICD (if we think of ICD as gold standard)
GWASCatSNO = GWASCat[(is.na(GWASCat$ICD)), ] #roughly half of the rows removed, but only 180 unique ICD9s here

#European participants
GWASCatICD = GWASCatICD[c(GWASCatICD$PUBMEDID %in% AncList), ] #161 left
GWASCatSNO = GWASCatSNO[c(GWASCatSNO$PUBMEDID %in% AncList), ] #161 left

#Snp format, split by semi-colon
GWASCatICD = data.table(GWASCatICD)[, list(SNPS = unlist(strsplit(SNPS, '; '))), by=c('P-VALUE', 'ICD', 'SNOMEDCT')]
GWASCatICD = GWASCatICD[c(GWASCatICD$SNPS %in% SnpIDs$rsID), ] #only exome IDs
GWASCatICD$ICD = as.character(GWASCatICD$ICD)

GWASCatSNO = data.table(GWASCatSNO)[, list(SNPS = unlist(strsplit(SNPS, '; '))), by=c('P-VALUE', 'ICD', 'SNOMEDCT')]
GWASCatSNO = GWASCatSNO[c(GWASCatSNO$SNPS %in% SnpIDs$rsID), ]
GWASCatSNO$SNOMEDCT = as.character(GWASCatSNO$SNOMEDCT)

#Lots of duplicates
GWASCatSNO = unique(GWASCatSNO)
GWASCatICD = unique(GWASCatICD)

#Catch formatting
#if length of string before period is 2, add a 0 in front, if the length is 1, add two zeros in front
for(phe in 1:nrow(GWASCatICD)){
  if(nchar(gsub('\\.\\d?\\d?', '', GWASCatICD$ICD[phe]))==2){
    GWASCatICD$ICD[phe] = paste0('0', GWASCatICD$ICD[phe])
  }
  if(nchar(gsub('\\.\\d?\\d?', '', GWASCatICD$ICD[phe]))==1){
    GWASCatICD$ICD[phe] = paste0('00', GWASCatICD$ICD[phe])
  }
} 


GWASCatICD = merge(GWASCatICD, phemap, by.x = 'ICD', by.y = 'ICD9')
GWASCatSNO = merge(GWASCatSNO, snoIDs, by.x = 'SNOMEDCT', by.y = 'SNOMED_CID')
#GWASCatSNO = unique(GWASCatSNO)
GWASCatSNO = merge(GWASCatSNO, phemap, by.x = 'ICD_CODE', by.y = 'ICD9') #73 unique phecodes

#GWASCatICD  = dplyr::select(GWASCatICD, ICD, 'P-VALUE', SNOMEDCT, SNPS, PheCode)
#For now leave as is.. but keep in mind p-values go down to 9e-06
GWASCatICD = GWASCatICD %>% rename(SNP = SNPS, jd_code = PheCode) %>% dplyr::select(SNP, jd_code)
GWASCatSNO = GWASCatSNO %>% rename(SNP = SNPS, jd_code = PheCode) %>% dplyr::select(SNP, jd_code)

GWASCatSNO = GWASCatSNO[!(is.na(GWASCatSNO$jd_code)), ]
GWASCatICD = GWASCatICD[!(is.na(GWASCatICD$jd_code)), ]

#Then compare to above
GWASCatICD = unique(GWASCatICD)
GWASCatSNO = unique(GWASCatSNO)

#Bring matches together
GWASJoin = anti_join(GWASCatSNO, GWASCatICD, by = c('SNP', 'jd_code')) #anti-join -> adding unique elements of GWASCatSNO
GWASBoth = rbind(GWASCatICD, GWASJoin)

#Adapt from Exome_GWASCat_AltROCK.R
Coxfiles = list.files(file.path(filepath, 'GWASCatalog', 'exome_full'))[c(grep('_cox.tsv.gz', list.files(file.path(filepath, 'GWASCatalog', 'exome_full'))))]
Logfiles = list.files(file.path(filepath, 'GWASCatalog', 'exome_full'))[c(grep('_logistic.tsv.gz', list.files(file.path(filepath, 'GWASCatalog', 'exome_full'))))]

#Get file for phecode of interest, and make sure its matched between cox and logistic
FilePhecodes = gsub('exome_|_cox.tsv.gz', '', Coxfiles)

MatchPhestestICD = unique(gsub('\\.', 'p', GWASCatICD$jd_code))
for(item in 1:length(MatchPhestestICD)){
  #Matches to phecode filenames require removing starting zeros
  if(substr(MatchPhestestICD[item], start = 1, stop = 1)=='0'){
    MatchPhestestICD[item] = substr(MatchPhestestICD[item], start = 2, stop = nchar(MatchPhestestICD[item]))
  }
}

MatchPhestestSNO = unique(gsub('\\.', 'p', GWASCatSNO$jd_code))
for(item in 1:length(MatchPhestestSNO)){
  #Matches to phecode filenames require removing starting zeros
  if(substr(MatchPhestestSNO[item], start = 1, stop = 1)=='0'){
    MatchPhestestSNO[item] = substr(MatchPhestestSNO[item], start = 2, stop = nchar(MatchPhestestSNO[item]))
  }
}

MatchPhestestJoin = unique(gsub('\\.', 'p', GWASBoth$jd_code))
for(item in 1:length(MatchPhestestJoin)){
  #Matches to phecode filenames require removing starting zeros
  if(substr(MatchPhestestJoin[item], start = 1, stop = 1)=='0'){
    MatchPhestestJoin[item] = substr(MatchPhestestJoin[item], start = 2, stop = nchar(MatchPhestestJoin[item]))
  }
}

MatchPhestestICD = paste0('phe', MatchPhestestICD)
MatchPhestestSNO = paste0('phe', MatchPhestestSNO)
MatchPhestestJoin = paste0('phe', MatchPhestestJoin)

FilePhecodesICD = FilePhecodes[FilePhecodes %in% MatchPhestestICD] 
FilePhecodesSNO = FilePhecodes[FilePhecodes %in% MatchPhestestSNO] 
FilePhecodesBoth = FilePhecodes[FilePhecodes %in% MatchPhestestJoin] 
FilePhecodesBoth = unique(FilePhecodesBoth)

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

#Build ROC DFs
ICDROC = PhePhetchROC(FilePhecodesICD, GWASCatICD)
SNOROC = PhePhetchROC(FilePhecodesSNO, GWASCatSNO)
BothROC = PhePhetchROC(FilePhecodesBoth, GWASBoth)

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

#GGplots
OutputICD = BuildROCCurve(ICDROC, 'GWAS Catalog v1.0.2 - ICD Map')
OutputSNO = BuildROCCurve(SNOROC, 'GWAS Catalog v1.0.2 - SNOMEDCT Map')
OutputBoth = BuildROCCurve(BothROC, 'GWAS Catalog v1.0.2 - ICD and SNOMEDCT Map')

ggsave(filename = 'MultiPheROC_ICD_GWASv102.pdf', path = file.path(filepath, 'GWASCatalog'), plot = OutputICD, width = 5, height = 3.5, units = c('in'))
ggsave(filename = 'MultiPheROC_ICD_GWASv102.png', path = file.path(filepath, 'GWASCatalog'), plot = OutputICD, width = 5, height = 3.5, units = c('in'), dpi = 300)
ggsave(filename = 'MultiPheROC_SNO_GWASv102.pdf', path = file.path(filepath, 'GWASCatalog'), plot = OutputSNO, width = 5, height = 3.5, units = c('in'))
ggsave(filename = 'MultiPheROC_SNO_GWASv102.png', path = file.path(filepath, 'GWASCatalog'), plot = OutputSNO, width = 5, height = 3.5, units = c('in'), dpi = 300)
ggsave(filename = 'MultiPheROC_ICD_SNO_GWASv102.pdf', path = file.path(filepath, 'GWASCatalog'), plot = OutputBoth, width = 5, height = 3.5, units = c('in'))
ggsave(filename = 'MultiPheROC_ICD_SNO_GWASv102.png', path = file.path(filepath, 'GWASCatalog'), plot = OutputBoth, width = 5, height = 3.5, units = c('in'), dpi = 300)
