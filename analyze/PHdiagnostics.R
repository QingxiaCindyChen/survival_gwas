#CoxPH diagnostics - run through the SNPs for a given phecode, and also calculate schoenfelds and see if anything correlates
#Inputs: catalog_phecode_map.txt, exome_phe555p1.csv.gz, ExomeIDs.csv
#library('tidyverse')
library('cowplot')
library('lubridate')
library('survival')
library('snpStats')
library('survminer')
library('readr')
library('dplyr')

#filepath = file.path('/Users', 'srhoades', 'Dropbox (VUMC)/', 'JakeColab_CountingGWAS', 'counting_gwas')
filepath = file.path('/home', 'rhoadesd', 'counting_gwas')
resultDir = 'results'

source(file.path(filepath, 'analyze', 'analyze_setup.R'))

filePrefix = 'exome'

minRecLen = 0 # years
minEvents = 2 # to be defined as a case; controls have zero events
maxAgeAtEvent = 90 # years; dates in the SD after this age are untrustworthy
minCases = 50 # to analyze the phecode

nPC = 2 # for both cox and logistic regression
splineDf = 4 # for last_age in logistic regression (ns::splines)
buffer = 1 # years; for calculating age1 in cox regression

minMaf = 0.01
minCallRate = 0.95
minHwePval = 0.001

############################################################
# load snp data

snpData = setDT(read_csv(file.path(filepath, 'analyze', procDir, 'exome_map.csv.gz'), col_types = cols()))
snpData = snpData[,.(snp = snp.name, chr = chromosome, pos = position)]

genoData = readRDS(file.path(filepath, 'analyze', procDir, sprintf('%s_genotype_data.rds', filePrefix)))

idx = (genoData$genoSummary$MAF >= minMaf) &
  (genoData$genoSummary$Call.rate >= minCallRate) &
  (2 * pnorm(-abs(genoData$genoSummary$z.HWE)) >= minHwePval)
snps = intersect(colnames(genoData$genoFull$genotypes)[idx], snpData[chr <= 22, snp])

############################################################
# load grid data

gridData = setDT(read_csv(file.path(filepath, 'analyze', procDir, sprintf('%s_grid_data.csv.gz', filePrefix)),
                          col_types = 'ccDDD'))
gridData[, first_age := time_length(first_entry_date - dob, 'years')]
gridData[, last_age := time_length(last_entry_date - dob, 'years')]
gridData[, rec_len := last_age - first_age]
gridData = gridData[first_age >= 0 & rec_len >= minRecLen,]

gridData = cbind(gridData, data.table(splines::ns(gridData$last_age, df = splineDf)))
colnames(gridData)[(ncol(gridData) - splineDf + 1):ncol(gridData)] = paste0('last_age', 1:splineDf)

# load PC data
pcData = setDT(read_csv(file.path(filepath, 'analyze', procDir, sprintf('%s_pc_data.csv.gz', filePrefix)), col_types = cols()))
if (nPC > 0) {
  gridData = merge(gridData, pcData[, 1:(1 + nPC)], by = 'grid')}

genoTmp = data.table(grid = genoData$genoFull$fam$pedigree, sex = genoData$genoFull$fam$sex)
gridData = merge(gridData, genoTmp, by = 'grid')

############################################################
# load phenotype data

phenoData = setDT(read_csv(file.path(filepath, 'analyze', procDir, sprintf('%s_phenotype_data.csv.gz', filePrefix)),
                           col_types = 'ccD'))

phenoData = merge(phenoData, gridData[, .(grid, dob, sex)], by = 'grid')
phenoData[, age := time_length(entry_date - dob, 'years')]
phenoData = phenoData[age <= maxAgeAtEvent]

phenoData = merge(phenoData, phecodeData[, .(phecode, whichSex)], by = 'phecode')
phenoData = phenoData[(whichSex == 'both') | (whichSex == 'male' & sex == 1) |
                        (whichSex == 'female' & sex == 2)]

phenoTmp = phenoData[, .N, by = .(grid, phecode)]
phenoTmp = phenoTmp[, .(nCases = sum(N >= minEvents)), by = phecode]
phenoData = merge(phenoData, phenoTmp[nCases >= minCases], by = 'phecode')

phenoData = phenoData[, .(grid, phecode, age)]

# testing
# phenoData = phenoData[phecode %in% c('274.1', '714.1', '335', '185', '427.21', '290.11')]

############################################################

#Testing
#phecodeDataKeep = phecodeData[phecode %in% c('290.11'), .(phecode, whichSex)]
phecodeDataKeep = phecodeData[phecode %in% c('274.1', '714.1', '335', '185', '427.21', '290.11'), .(phecode, whichSex)]
phecodeDataKeep[, filename := sprintf('%s_phe%s.csv.gz', filePrefix, gsub('.', 'p', phecode, fixed = TRUE))]

# genoFlat = data.table(grid = gridData$grid, as(genoData$genoFull$genotypes[gridData$grid, snps], 'numeric'))
# genoFlat = melt(genoFlat, id.vars = 'grid', variable.name = 'snpId', value.name = 'snp')
# genoFlat = genoFlat[!is.na(snp)]
# genoFlat[, snp := as.integer(snp)]

methodVec = rep.int(c('logistic', 'cox'), times = length(snps))
snpVec = rep(snps, each = 2)

registerDoParallel(cores = 20)

logFilepath = file.path(filepath, 'analyze', resultDir, sprintf('%s_progress.txt', filePrefix))
timeStarted = Sys.time()
cat(sprintf('%s started analysis\n', timeStarted), file = logFilepath)

# microbenchmark::microbenchmark({

#test on 100 snps
#snps = snps[1:50]
schoens = data.frame(snp = c(snps), method = c(rep('cox', length(snps))), schoenfeld = c(rep(NA, length(snps))))

done = foreach(ii = 1:nrow(phecodeDataKeep)) %dopar% {
  whichSex = phecodeDataKeep$whichSex[ii]
  phenoDataNow = phenoData[phecode == phecodeDataKeep$phecode[ii], .(grid, age)]
  
  inputBase = makeInput(phenoDataNow, gridData, whichSex, minEvents, buffer)
  glmStr = getGlmStr(whichSex, nPC, splineDf)
  coxStr = getCoxStr(whichSex, nPC)
  
  resultNow = foreach(snp = snps, .combine = rbind) %do% {
    inputNow = addSnpToInput(inputBase, genoData$genoFull, snp)
    # inputNow = addSnpToInput1(inputBase, genoFlat, snp)
    glmFit = runGlm(glmStr, inputNow) #don't need but won't mess with dt structure, it will add time though
    coxFit = runCox(coxStr, inputNow)
    #Schoenfelds
    testPH = cox.zph(coxFit)
    schoens$schoenfeld[schoens$snp==snp] = testPH$table['snp', 'p']

    data.table(rbind(coef(summary(glmFit))[2, 1:3], coef(summary(coxFit))[1, c(1, 3:4)]))}

  colnames(resultNow) = c('coef', 'se', 'z')
  
  # speedglm encodes pvals as factors, and coxph sets small pvals to zero.
  resultNow[, pval := 2 * pnorm(-abs(z))]
  resultNow[, method := methodVec]
  resultNow[, snp := snpVec]
  resultNow[, phecode := phecodeDataKeep$phecode[ii]]
  
  SchoenTest = merge(resultNow, schoens, by = c('snp', 'method'))
  
  SchoenTest = SchoenTest[c(SchoenTest$pval < 0.05 | SchoenTest$schoenfeld < 0.05), ]
  SchoenTest$log10pval = -log10(SchoenTest$pval)
  SchoenTest$log10schoen = -log10(SchoenTest$schoenfeld)
  
  maxval = max(max(SchoenTest$log10pval), max(SchoenTest$log10schoen))
  
  SchoenCompare = ggplot(data = SchoenTest, aes(x = log10pval, y = log10schoen))+
    geom_point(size = 0.5)+
    geom_abline(slope = 1, intercept = 0, color = 'black') +
    theme_bw()+
    xlim(c(0, maxval)) +
    ylim(c(0, maxval)) + 
    #xlab('-log10(CoxPH p')+
    xlab(expression('CoxPH -log'[10]~'(p)')) +
    #ylab('Schoenfeld p')+
    ylab(expression('Schoenfeld -log'[10]~'(p)')) +
    ggtitle(paste('Phecode', phecodeDataKeep$phecode[ii]))+
    theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11, face = "bold"), 
          panel.grid.major = element_line(colour = 'grey98'), panel.grid.minor = element_line(colour = 'grey98'),
          legend.position = 'none', title = element_text(size = 11, face = "bold"), 
          plot.title = element_text(hjust = 0.5))
  ggsave(filename = sprintf('Phe%sPHAssume.pdf', gsub('\\.', 'p', phecodeDataKeep$phecode[ii])), path = file.path(filepath, 'analyze', resultDir), plot = SchoenCompare, width = 5, height = 4, units = c('in'))
}

