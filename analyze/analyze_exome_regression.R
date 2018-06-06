source(file.path('analyze', 'analyze_setup.R'))

testing = FALSE
if (testing) {
  resultDir = 'results_test'
}

registerDoParallel(cores = 20)
filePrefix = 'exome'

plinkFilepath = '../plink-1.90b6/plink'
plinkMemSize = 100 # MB
plinkVif = 1e7 # reduce plink's tendency to output NA due to multicollinearity of the spline bases

genoDir = '../genotype_data/exome'
genoPrefix = 'Exome_GRID_Euro'

minRecLen = 0 # years
minEvents = 2 # to be defined as a case; controls have zero events
maxAgeAtEvent = 90 # years; dates in the SD after this age are untrustworthy
minCases = 50 # to analyze the phecode

nPC = 2 # for both cox and logistic regression
splineDf = 3 # for last_age in logistic regression (ns::splines)
buffer = 1 # years; for calculating age1 in cox regression

minMaf = 0.01
minCallRate = 0.95
minHwePval = 0.001

############################################################
# load snp data

snpData = setDT(read_csv(file.path(procDir, 'exome_map.csv.gz'), col_types = cols()))
snpData = snpData[, .(snp = snp.name, chr = chromosome, pos = position)]

genoData = readRDS(file.path(procDir, sprintf('%s_genotype_data.rds', filePrefix)))

idx = (genoData$genoSummary$MAF >= minMaf) &
  (genoData$genoSummary$Call.rate >= minCallRate) &
  (2 * pnorm(-abs(genoData$genoSummary$z.HWE)) >= minHwePval)
snps = intersect(colnames(genoData$genoFull$genotypes)[idx], snpData[chr <= 22, snp])

if (testing) {
  snps = unique(read_csv('results/exome_pilot_top20.csv', col_types = cols())$snp)
  # set.seed(17)
  # snps = snps[sample.int(length(snps), 1000)]
}

snpFilepath = tempfile('snp_', fileext = '.tsv')
write_tsv(data.table(snps), snpFilepath, col_names = FALSE)

############################################################
# load grid data

gridData = setDT(read_csv(file.path(procDir, sprintf('%s_grid_data.csv.gz', filePrefix)),
                          col_types = 'ccDDD'))
gridData[, first_age := time_length(first_entry_date - dob, 'years')]
gridData[, last_age := time_length(last_entry_date - dob, 'years')]
gridData[, rec_len := last_age - first_age]
gridData = gridData[first_age >= 0 & rec_len >= minRecLen]

gridData = cbind(gridData, data.table(splines::ns(gridData$last_age, df = splineDf)))
colnames(gridData)[(ncol(gridData) - splineDf + 1):ncol(gridData)] = paste0('last_age', 1:splineDf)

covarColnames = c('rec_len', paste0('last_age', 1:splineDf))

# load PC data
if (nPC > 0) {
  pcData = setDT(read_csv(file.path(procDir, sprintf('%s_pc_data.csv.gz', filePrefix)), col_types = cols()))
  gridData = merge(gridData, pcData[, 1:(1 + nPC)], by = 'grid')
  covarColnames = c(covarColnames, colnames(pcData)[2:(1 + nPC)])}

genoTmp = data.table(grid = genoData$genoFull$fam$pedigree, sex = genoData$genoFull$fam$sex)
gridData = merge(gridData, genoTmp, by = 'grid')

covarData = gridData[, .(FID = grid, IID = grid)]
covarData = cbind(covarData, gridData[, c(covarColnames, 'sex'), with = FALSE])

covarFilepath = tempfile('covar_', fileext = '.tsv')
write_tsv(covarData, covarFilepath)

############################################################
# load phenotype data

phenoData = setDT(read_csv(file.path(procDir, sprintf('%s_phenotype_data.csv.gz', filePrefix)),
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

if (testing) {
  phenoData = phenoData[phecode %in% c('274.1', '714.1', '335', '185', '427.21', '290.11')]
}

############################################################

phecodeDataKeep = phecodeData[phecode %in% unique(phenoData$phecode), .(phecode, whichSex)]
phecodeDataKeep[, phecodeStr := paste0('phe', gsub('.', 'p', phecode, fixed = TRUE))]
phecodeDataKeep[, coxFilename := sprintf('%s_%s_cox.tsv.gz', filePrefix, phecodeStr)]
phecodeDataKeep[, logisticFilename := sprintf('%s_%s_logistic.tsv', filePrefix, phecodeStr)]
phecodeDataKeep[, covarNum := paste0('1-', ncol(covarData) - ifelse(whichSex == 'both', 2, 3))]
phecodeDataKeep = merge(phecodeDataKeep, phenoTmp, by = 'phecode')

# split analysis across servers
# phecodeDataKeep = phecodeDataKeep[1:round(nrow(phecodeDataKeep) / 2)]
# phecodeDataKeep = phecodeDataKeep[(1 + round(nrow(phecodeDataKeep) / 2)):nrow(phecodeDataKeep)]
# change names of progress files and workspace file

############################################################
# run cox regression

logFilepath = file.path(resultDir, sprintf('%s_progress_cox.txt', filePrefix))
timeStarted = Sys.time()
cat(sprintf('%s started analysis\n', timeStarted), file = logFilepath)

phenoLogisticList = foreach(ii = 1:nrow(phecodeDataKeep)) %dopar% {
  whichSex = phecodeDataKeep$whichSex[ii]
  phenoDataNow = phenoData[phecode == phecodeDataKeep$phecode[ii], .(grid, age)]

  inputBase = makeInput(phenoDataNow, gridData, whichSex, minEvents, buffer)
  coxStr = makeCoxStr(whichSex, nPC)

  resultNow = foreach(snp = snps, .combine = rbind) %do% {
    inputNow = addSnpToInput(inputBase, genoData$genoFull, snp)
    coxFit = runCox(coxStr, inputNow)
    data.table(coef(summary(coxFit))[1, c(1, 3:4), drop = FALSE])}

  colnames(resultNow) = c('beta', 'se', 'z')
  resultNow[, pval := 2 * pnorm(-abs(z))] # coxph sets small pvals to zero
  resultNow[, snp := snps]
  write_tsv(resultNow, gzfile(file.path(resultDir, phecodeDataKeep$coxFilename[ii])))

  cat(sprintf('%s completed phecode %s (%d of %d)\n', Sys.time(), phecodeDataKeep$phecode[ii],
              ii, nrow(phecodeDataKeep)), file = logFilepath, append = TRUE)

  # prepare phenotype data for plink
  phenoLogisticNow = inputBase[, .(grid, status)]
  colnames(phenoLogisticNow)[2] = phecodeDataKeep$phecodeStr[ii]
  setkey(phenoLogisticNow, 'grid')
  phenoLogisticNow}

timeElapsed = Sys.time() - timeStarted
cat(sprintf('Time elapsed of %.2f %s\n', timeElapsed, attr(timeElapsed, 'units')),
    file = logFilepath, append = TRUE)

############################################################
# run logistic regression in plink

# make phenotype file for plink
phenoPlink = Reduce(function(...) merge(..., all = TRUE), phenoLogisticList)
phenoPlink[is.na(phenoPlink)] = -9

colnames(phenoPlink)[1] = 'FID'
phenoPlink[, IID := FID]
setcolorder(phenoPlink, c(1, ncol(phenoPlink)))

phenoFilepath = tempfile('pheno_', fileext = '.tsv')
write_tsv(phenoPlink, phenoFilepath)

plinkArgs = sprintf('%s --bfile %s --extract %s --covar %s --pheno %s --1 --memory %d --vif %d',
                    '--logistic hide-covar beta --ci 0.95', file.path(genoDir, genoPrefix),
                    snpFilepath, covarFilepath, phenoFilepath, plinkMemSize, plinkVif)

logFilepath = file.path(resultDir, sprintf('%s_progress_logistic.txt', filePrefix))
timeStarted = Sys.time()
cat(sprintf('%s started analysis\n', timeStarted), file = logFilepath)

done = foreach(ii = 1:nrow(phecodeDataKeep)) %dopar% {
  # run plink
  plinkArgsNow = sprintf('%s --pheno-name %s --covar-number %s --out %s_%s',
                         plinkArgs, phecodeDataKeep$phecodeStr[ii], phecodeDataKeep$covarNum[ii],
                         file.path(resultDir, filePrefix), phecodeDataKeep$phecodeStr[ii])
  system2(plinkFilepath, plinkArgsNow)

  # fix plink's stupid output spacing
  outputFilename = sprintf('%s_%s.assoc.logistic', filePrefix, phecodeDataKeep$phecodeStr[ii])
  outputFilepath = file.path(resultDir, outputFilename)
  tmpFilepath = tempfile(paste0(phecodeDataKeep$phecodeStr[ii], '_'))
  system(sprintf("cat %s | tr -s ' ' '\t' > %s", outputFilepath, tmpFilepath))
  system(sprintf("cat %s | sed 's/^[[:space:]]*//g' | sed 's/[[:space:]]*$//g' > %s",
                 tmpFilepath, outputFilepath))
  unlink(tmpFilepath)

  # rename and compress
  file.rename(outputFilepath, file.path(resultDir, phecodeDataKeep$logisticFilename[ii]))
  system2('gzip', paste('-f', file.path(resultDir, phecodeDataKeep$logisticFilename[ii])))

  cat(sprintf('%s completed phecode %s (%d of %d)\n', Sys.time(), phecodeDataKeep$phecode[ii],
              ii, nrow(phecodeDataKeep)), file = logFilepath, append = TRUE)}

phecodeDataKeep[, logisticFilename := paste0(logisticFilename, '.gz')]
# unlink(c(snpFilepath, covarFilepath, phenoFilepath))

timeElapsed = Sys.time() - timeStarted
cat(sprintf('Time elapsed of %.2f %s\n', timeElapsed, attr(timeElapsed, 'units')),
    file = logFilepath, append = TRUE)

############################################################

d = c('genoData', 'genoTmp', 'idx', 'pcData', 'phenoTmp', 'phenoLogisticList',
      'done', 'timeStarted', 'timeElapsed', 'phenoPlink')
save(list = setdiff(ls(), c(d, 'd')),
     file = file.path(resultDir, sprintf('%s_workspace.Rdata', filePrefix)))
