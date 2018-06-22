#Comparisons of coxph and gwas catalogue on all phecodes for a collective ROC
#Inputs: catalog_phecode_map.txt, ExomeIDs.csv, all exome_[phecode].csv.gz results
library('tidyverse')
library('cowplot')
library('lubridate')
library('survival')
library('snpStats')
library('survminer')
options(warn = -1)


filepath = file.path('/Users', 'srhoades', 'Dropbox (VUMC)/', 'JakeColab_CountingGWAS', 'counting_gwas', 'analyze', 'GWASCatalog')
#filepath = file.path('..')

GWASCat = read_delim(file = file.path(filepath, 'catalog_phecode_map.txt'), delim = '\t')

#Exome IDs
SnpIDs = read_csv(file = file.path(filepath, 'ExomeIDs.csv'))

#Only keep the GWAS Cat SNPs - we can only compare up to 589 SNPs to the GWAS Cat with our coxph approach
GWASCat = GWASCat[c(GWASCat$SNP %in% SnpIDs$rsID), ]
#Count number of SNPs by jd_code
#GWASCatTally = GWASCat %>% group_by(jd_code) %>% tally()

Coxfiles = list.files(file.path(filepath))[c(grep('_cox.tsv.gz', list.files(file.path(filepath))))]
Logfiles = list.files(file.path(filepath))[c(grep('_logistic.tsv.gz', list.files(file.path(filepath))))]

#Get file for phecode of interest, and make sure its matched between cox and logistic
FilePhecodes = gsub('exome_|_cox.tsv.gz', '', Coxfiles)

#Save yourself some memory/time by only considering the phecodes which exist on the GWASCat
MatchPhes = gsub('\\.', 'p', GWASCat$jd_code)
MatchPhes = paste0('phe', MatchPhes)
FilePhecodes = FilePhecodes[FilePhecodes %in% MatchPhes]

#Now we're interested in recall of the GWAS catalog
#How many Exome SNPs via Coxph are in the gold standard?
#True positives = hits on both; false positives = hits only on CoxPH; true negatives = not sig on both; false negatives = hits only on GWAS Cat

ROCLobster = function(SurvData, GWASCatPhenoSpecific, threshold, AllGWASCat){
  TP = length(SurvData$rsID[SurvData$rsID %in% GWASCatPhenoSpecific$SNP & SurvData$pval < threshold]) #snps in both significant lists below our pval
  #FP = length(SurvData$rsID[!(SurvData$rsID %in% GWASCatPhenoSpecific$SNP) & SurvData$pval < threshold]) #snps not found in the gwas catalogue but below our pval
  FP = length(SurvData$rsID[!(SurvData$rsID %in% GWASCatPhenoSpecific$SNP) & SurvData$rsID %in% GWASCat$SNP & SurvData$pval < threshold]) #snps not found in the gwas catalogue but below our pval
  #TN = length(SurvData$rsID[!(SurvData$rsID %in% GWASCatPhenoSpecific$SNP) & SurvData$pval > threshold]) #snps not found in the gwas catalogue but above pval
  TN = length(SurvData$rsID[!(SurvData$rsID %in% GWASCatPhenoSpecific$SNP) & SurvData$rsID %in% GWASCat$SNP & SurvData$pval > threshold]) #snps not found in the gwas catalogue but above pval
  FN = length(SurvData$rsID[(SurvData$rsID %in% GWASCatPhenoSpecific$SNP) & SurvData$pval > threshold])
  #TPR = TP/(TP+FN) #y axis - orig
  #FPR = FP/(FP+TN) #x axis - orig
  #Out = data.frame(TPR = TPR, FPR = FPR) - orig
  Out = data.frame(TP = TP, FP = FP, TN = TN, FN = FN) #new - for looping through phecodes
  return(Out)
}

ROCPlobster = function(PlotDF, PlotTitle){
  ROC = ggplot(PlotDF, aes(x = FPR, y = TPR, group = test, colour = test))+
    scale_x_continuous(expand = c(0.0075, 0.0075)) + 
    scale_y_continuous(expand = c(0.005, 0.005)) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), lty = 2, color = 'black')+
    geom_line(size = 1)+
    theme_bw()+
    #ggtitle(paste('Phecode', PlotTitle))+
    ggtitle(NULL)+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11, face = "bold"), 
          #panel.grid.major = element_line(colour = 'grey98'), panel.grid.minor = element_line(colour = 'grey98'),
          legend.title =  element_blank(), legend.text = element_text(size = 10), 
          legend.box.margin = margin(-2, -2, -2, -2), title = element_text(size = 11, face = "bold"), 
          plot.title = element_text(hjust = 0.5))
  return(ROC)
}

#Generating ROC requires looping through thresholds, but within each threshold, we need to compile every phecodes' result
thresholds = seq(0, 1, 0.01) #This is the major determinant of speed - how high of plotting resolution do you want?
fillercox = as.data.frame(matrix(NA, nrow = length(thresholds), ncol = 2))
colnames(fillercox) = c('TPR', 'FPR')
fillerlog = as.data.frame(matrix(NA, nrow = length(thresholds), ncol = 2))
colnames(fillerlog) = c('TPR', 'FPR')

for(thresh in 1:length(thresholds)){
  tempconfusioncox = as.data.frame(matrix(NA, nrow = length(FilePhecodes), ncol = 4))
  colnames(tempconfusioncox) = c('TP', 'FP', 'TN', 'FN')
  tempconfusionlog = as.data.frame(matrix(NA, nrow = length(FilePhecodes), ncol = 4))
  colnames(tempconfusionlog) = c('TP', 'FP', 'TN', 'FN')
  
  for(phe in 1:length(FilePhecodes)){
    CoxDF = read_tsv(file = file.path(filepath, sprintf('exome_%s_cox.tsv.gz', FilePhecodes[phe])), col_types = cols())
    LogMatch = Logfiles[c(grep(sprintf('exome_%s_logistic.tsv.gz', FilePhecodes[phe]), Logfiles))]
    LogDF = read_tsv(file = file.path(filepath, LogMatch), col_types = cols())
    
    #ID match
    CoxDF = merge(CoxDF, SnpIDs, by.x = 'snp', by.y = 'SNP')
    LogDF = merge(LogDF, SnpIDs, by.x = 'SNP', by.y = 'SNP')
    #rename
    LogDF = rename(LogDF, pval = P)
    
    #Get 'raw' phecode ID
    Phecode = gsub('phe', '', FilePhecodes[phe])
    Phecode = gsub('p', '.', Phecode)
    
    #Select SNPs from GWAS Catalog
    GWASCatMatch = GWASCat[GWASCat$jd_code==Phecode, ]
    
    tempphecox = ROCLobster(CoxDF, GWASCatMatch, thresholds[thresh], GWASCat)
    tempphelog = ROCLobster(LogDF, GWASCatMatch, thresholds[thresh], GWASCat)
    tempconfusioncox[phe, ] = tempphecox
    tempconfusionlog[phe, ] = tempphelog
  }
  
  #Now get collective TPR and FPR
  TPRcox = sum(tempconfusioncox$TP)/(sum(tempconfusioncox$TP) + sum(tempconfusioncox$FN))
  FPRcox = sum(tempconfusioncox$FP)/(sum(tempconfusioncox$FP) + sum(tempconfusioncox$TN)) #y axis - orig
  CoxAdd = data.frame(TPR = TPRcox, FPR = FPRcox)
  fillercox[thresh, ] = CoxAdd
  
  TPRlog = sum(tempconfusionlog$TP)/(sum(tempconfusionlog$TP) + sum(tempconfusionlog$FN))
  FPRlog = sum(tempconfusionlog$FP)/(sum(tempconfusionlog$FP) + sum(tempconfusionlog$TN)) #y axis - orig
  LogAdd = data.frame(TPR = TPRlog, FPR = FPRlog)
  fillerlog[thresh, ] = LogAdd
}

fillercox$test = c('CoxPH')
fillerlog$test = c('Logistic')
filler = rbind(fillercox, fillerlog)

Outplot = ROCPlobster(filler, PlotTitle = '')
ggsave(filename = 'MultiPhecodesROC.pdf', path = file.path(filepath), plot = Outplot, width = 5, height = 3.5, units = c('in'))

