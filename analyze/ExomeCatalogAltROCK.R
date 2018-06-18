#Alternative format - using data.table, and building a table to use plotROC
#Inputs: catalog_phecode_map.txt, ExomeIDs.csv, all exome_[phecode].csv.gz results
library('readr')
library('cowplot')
library('data.table')
library('survival')
library('snpStats')
library('survminer')
library('dplyr')
library('plotROC')
options(warn = -1)


filepath = file.path('/Users', 'srhoades', 'Dropbox (VUMC)/', 'JakeColab_CountingGWAS', 'counting_gwas', 'analyze', 'GWASCatalog', 'exome_full')
#filepath = file.path('/home', 'rhoadesd', 'counting_gwas', 'analyze', 'results', 'GWASCatalog', 'exome_full')

#GWASCat = read_delim(file = file.path(filepath, 'catalog_phecode_map.txt'), delim = '\t')
GWASCat = setDT(read_delim(file = file.path(filepath, 'catalog_phecode_map.txt'), delim = '\t', col_types = cols()))

#Exome IDs
#SnpIDs = read_csv(file = file.path(filepath, 'ExomeIDs.csv'))
SnpIDs = setDT(read_csv(file = file.path(filepath, 'ExomeIDs.csv'), col_types = cols()))

#Only keep the GWAS Cat SNPs - we can only compare up to 589 SNPs to the GWAS Cat with our coxph approach
GWASCat = GWASCat[c(GWASCat$SNP %in% SnpIDs$rsID), ]
#Count number of SNPs by jd_code
GWASCatTally = GWASCat %>% group_by(jd_code) %>% tally()

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
#Example
#RefPhecode = FilePhecodes[3]
ROC = NULL
for(code in FilePhecodes){
  CoxDF = setDF(read_tsv(file = file.path(filepath, sprintf('exome_%s_cox.tsv.gz', code)), col_types = cols()))
  LogMatch = Logfiles[c(grep(sprintf('exome_%s_logistic.tsv.gz', code), Logfiles))]
  LogDF = setDF(read_tsv(file = file.path(filepath, LogMatch), col_types = cols()))
  LogDF = rename(LogDF, pval = P)
  LogDF = rename(LogDF, snp = SNP)

  #Get 'raw' phecode ID
  Phecode = gsub('phe', '', code)
  Phecode = gsub('p', '.', Phecode)

  #Only take what you need
  CoxDF = select(CoxDF, snp, pval)
  CoxDF$phecode = Phecode
  CoxDF$method = 'Cox'
  CoxDF = merge(CoxDF, SnpIDs, by.x = 'snp', by.y = 'SNP')
  CoxDF = select(CoxDF, -snp)
  CoxDF = CoxDF[CoxDF$rsID %in% GWASCat$SNP, ]
  temp1 = inner_join(CoxDF, GWASCat, by = c('rsID' = 'SNP', 'phecode' = 'jd_code'))
  temp1$GWAS = 1
  temp2 = anti_join(CoxDF, GWASCat, by = c('rsID' = 'SNP', 'phecode' = 'jd_code'))
  temp2$GWAS = 0

  CoxDF = rbind(temp1, temp2)

  LogDF = select(LogDF, snp, pval)
  LogDF$phecode = Phecode
  LogDF$method = 'Logistic'
  LogDF = merge(LogDF, SnpIDs, by.x = 'snp', by.y = 'SNP')
  LogDF = select(LogDF, -snp)
  LogDF = LogDF[LogDF$rsID %in% GWASCat$SNP, ]
  temp1 = inner_join(LogDF, GWASCat, by = c('rsID' = 'SNP', 'phecode' = 'jd_code'))
  temp1$GWAS = 1
  temp2 = anti_join(LogDF, GWASCat, by = c('rsID' = 'SNP', 'phecode' = 'jd_code'))
  temp2$GWAS = 0

  LogDF = rbind(temp1, temp2)

  #ID match
  ComboDF = rbind(CoxDF, LogDF)

  ROCtemp = select(ComboDF, pval, method, GWAS)
  ROC = rbind(ROC, ROCtemp)
}

ROC$pval = as.numeric(ROC$pval)
#Opposite thresholding since we're using pvals
ROC$pval = 1 - ROC$pval
ROCPlot = ggplot(ROC, aes(d = GWAS, m = pval, colour = method)) + geom_roc(n.cuts = 0)+
  scale_x_continuous(expand = c(0.0075, 0.0075)) + 
  scale_y_continuous(expand = c(0.005, 0.005)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), lty = 2, color = 'black')+
  theme_bw()+
  xlab('FPF')+
  ylab('TPF')+
  #ggtitle(paste('Phecode', PlotTitle))+
  ggtitle(NULL)+
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11, face = "bold"), 
        #panel.grid.major = element_line(colour = 'grey98'), panel.grid.minor = element_line(colour = 'grey98'),
        legend.title =  element_blank(), legend.text = element_text(size = 10), 
        legend.box.margin = margin(-5, -5, -5, -5), title = element_text(size = 11, face = "bold"), 
        plot.title = element_text(hjust = 0.5))
ROCPlotOut = ROCPlot + annotate("text", x = 0.8, y = 0.1, label = paste("Cox AUC = ", round(calc_auc(ROCPlot)["AUC"][1, ], 2)))+
  annotate("text", x = 0.765, y = 0.05, label = paste("Logistic AUC = ", round(calc_auc(ROCPlot)["AUC"][2, ], 2)))

ggsave(filename = 'MultiPhecodesROC.pdf', path = file.path(filepath), plot = ROCPlotOut, width = 5, height = 3.5, units = c('in'))

