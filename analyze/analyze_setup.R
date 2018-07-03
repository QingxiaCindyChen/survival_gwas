library('data.table')
# library('speedglm')
library('readr')
library('lubridate')
library('doParallel')
library('snpStats')
library('cowplot')
library('qqman')
library('yaml')

procDir = 'processed'
# resultDir = 'results'

phecodeData = setDT(read_csv(file.path(procDir, 'phecode_data.csv.gz'), col_types = 'ccc??????'))

eb = element_blank()
theme_set(theme_light() +
            theme(axis.text = element_text(color = 'black'), strip.text = element_text(color = 'black'),
                  panel.grid.minor = eb, legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))

############################################################
# functions for loading data

loadGeno = function(procDir, filePrefix, p, snpSubsetPath = NULL) {
  genoData = readRDS(file.path(procDir, sprintf('%s_genotype_data.rds', filePrefix)))

  idx = (genoData$genoSummary$MAF >= p$minMaf) &
    (genoData$genoSummary$Call.rate >= p$minCallRate) &
    (2 * pnorm(-abs(genoData$genoSummary$z.HWE)) >= p$minHwePval) &
    (genoData$genoFull$map$chromosome <= 22)

  genoFull = genoData$genoFull
  genoFull$genotypes = genoFull$genotypes[, idx]
  genoFull$map = genoFull$map[idx,]

  if (!is.null(snpSubsetPath) & length(snpSubsetPath) > 0) {
    snps = unique(read_tsv(snpSubsetPath, col_types = 'c', col_names = FALSE)$X1)
    genoFull$genotypes = genoFull$genotypes[, colnames(genoFull$genotypes) %in% snps]
    genoFull$map = genoFull$map[genoFull$map$snp.name %in% snps,]}

  return(genoFull)}


loadGrid = function(procDir, filePrefix, minRecLen, splineDf, nPC, fam) {
  gridData = read_csv(file.path(procDir, sprintf('%s_grid_data.csv.gz', filePrefix)), col_types = 'ccDDD')
  setDT(gridData)
  gridData[, first_age := time_length(first_entry_date - dob, 'years')]
  gridData[, last_age := time_length(last_entry_date - dob, 'years')]
  gridData[, rec_len := last_age - first_age]
  gridData = gridData[first_age >= 0 & rec_len >= minRecLen]

  gridData = cbind(gridData, data.table(splines::ns(gridData$last_age, df = splineDf)))
  colnames(gridData)[(ncol(gridData) - splineDf + 1):ncol(gridData)] = paste0('last_age', 1:splineDf)
  covarColnames = c('rec_len', paste0('last_age', 1:splineDf))

  if (nPC > 0) {
    pcData = read_csv(file.path(procDir, sprintf('%s_pc_data.csv.gz', filePrefix)), col_types = cols())
    setDT(pcData)
    gridData = merge(gridData, pcData[, 1:(1 + nPC)], by = 'grid')
    covarColnames = c(covarColnames, colnames(pcData)[2:(1 + nPC)])}

  genoTmp = data.table(fam)[, .(grid = pedigree, sex)]
  gridData = merge(gridData, genoTmp, by = 'grid')

  return(list(gridData, covarColnames))}


loadPheno = function(procDir, filePrefix, p, gridData, phecodeSubsetPath) {
  phenoData = read_csv(file.path(procDir, sprintf('%s_phenotype_data.csv.gz', filePrefix)),
                       col_types = 'ccD')
  setDT(phenoData)

  if (!is.null(phecodeSubsetPath) & length(phecodeSubsetPath) > 0) {
    phecodes = unique(read_tsv(phecodeSubsetPath, col_types = 'c', col_names = FALSE)$X1)
    phenoData = phenoData[phecode %in% phecodes]}

  phenoData = merge(phenoData, gridData[, .(grid, dob, sex)], by = 'grid')
  phenoData[, age := time_length(entry_date - dob, 'years')]
  phenoData = phenoData[age <= p$maxAgeAtEvent]

  phenoData = merge(phenoData, phecodeData[, .(phecode, whichSex)], by = 'phecode')
  phenoData = phenoData[(whichSex == 'both') | (whichSex == 'male' & sex == 1) |
                          (whichSex == 'female' & sex == 2)]

  phenoSummary = phenoData[, .N, by = .(grid, phecode)]
  phenoSummary = phenoSummary[, .(nCases = sum(N >= p$minEvents)), by = phecode]
  phenoData = merge(phenoData, phenoSummary[nCases >= p$minCases], by = 'phecode')
  phenoData = phenoData[, .(grid, phecode, age)]
  return(list(phenoData, phenoSummary))}


############################################################
# functions for cox regression

makeGwasMetadata = function(filePrefix, phecodeData, phenoData, phenoSummary) {
  d = phecodeData[phecode %in% unique(phenoData$phecode), .(phecode, whichSex)]
  d[, phecodeStr := paste0('phe', gsub('.', 'p', phecode, fixed = TRUE))]
  d[, coxFilename := sprintf('%s_%s_cox.tsv.gz', filePrefix, phecodeStr)]
  d[, logisticFilename := sprintf('%s_%s_logistic.tsv', filePrefix, phecodeStr)]
  # d[, covarNum := paste0('1-', ncol(covarData) - ifelse(whichSex == 'both', 2, 3))]
  return(merge(d, phenoSummary, by = 'phecode'))}


# expected colnames in phenoData: grid, age
# expected colnames in gridData: grid, first_age, last_age
makeInput = function(phenoData, gridData, whichSex, minEvents, ageBuffer) {
  phenoCase = phenoData
  setkeyv(phenoCase, c('grid', 'age'))
  phenoCase = phenoCase[, if (.N >= minEvents) .SD[minEvents,], by = grid]
  phenoControl = fsetdiff(gridData[, .(grid)], phenoData[, .(grid)])

  input = rbind(phenoCase, phenoControl, fill = TRUE)
  input = merge(input, gridData, by = 'grid')

  if (whichSex == 'male') {
    input = input[sex == 1]
  } else if (whichSex == 'female') {
    input = input[sex == 2]}

  input[, status := ifelse(is.na(age), 0, 1)]
  input[, age2 := ifelse(status, age, last_age)]
  input[, age1 := min(first_age, max(0, age2 - ageBuffer)), by = grid]
  input}


addSnpToInput = function(input, genoFull, snp) {
  input[, genotype := as(genoFull$genotypes[grid, snp], 'numeric')[,1]]
  input[!is.na(genotype),]}


makeCoxStr = function(whichSex, nPC) {
  formStr = 'Surv(age1, age2, status) ~ genotype'
  if (whichSex == 'both') {
    formStr = paste(formStr, '+ sex')}
  if (nPC > 0) {
    formStr = paste(formStr, '+', paste0('PC', 1:nPC, collapse = ' + '))}
  formStr}


runCox = function(formulaStr, input) {
  coxph(formula(formulaStr), data = input)}



runGwasCox = function(inputBase, genoFull, formulaStr) {
  snps = colnames(genoFull$genotypes)
  d = foreach(snp = snps, .combine = rbind) %do% {
    inputNow = addSnpToInput(inputBase, genoFull, snp)
    coxFit = runCox(formulaStr, inputNow)
    data.table(coef(summary(coxFit))[1, c(1, 3:4), drop = FALSE])}
  colnames(d) = c('beta', 'se', 'z')
  d[, pval := 2 * pnorm(-abs(z))] # coxph sets small pvals to zero
  d[, snp := snps]
  return(d)}

############################################################
# functions for log files

createLogFile = function(resultDir, filePrefix, fileSuffix) {
  path = file.path(resultDir, sprintf('%s_progress_%s.txt', filePrefix, fileSuffix))
  timeStarted = Sys.time()
  cat(sprintf('%s started analysis\n', timeStarted), file = path)
  return(list(path = path, timeStarted = timeStarted))}


appendLogFile = function(logFile, gwasMetadata, ii) {
  cat(sprintf('%s completed phecode %s (%d of %d)\n', Sys.time(), gwasMetadata$phecode[ii],
              ii, nrow(gwasMetadata)), file = logFile$path, append = TRUE)}


finishLogFile = function(logFile) {
  timeElapsed = Sys.time() - logFile$timeStarted
  cat(sprintf('Time elapsed of %.2f %s\n', timeElapsed, attr(timeElapsed, 'units')),
      file = logFile$path, append = TRUE)}

############################################################
# functions for logistic regression using plink

makePhenoPlink = function(inputBase, phecodeStr) {
  phenoPlinkNow = inputBase[, .(grid, status)]
  colnames(phenoPlinkNow)[2] = phecodeStr
  setkey(phenoPlinkNow, 'grid')
  return(phenoPlinkNow)}


prepForPlink = function(snps, gridData, covarColnames, gwasMetadata, phenoPlinkList) {
  snpPath = tempfile('snp_', fileext = '.tsv')
  write_tsv(data.table(snps), snpPath, col_names = FALSE)

  covarData = gridData[, .(FID = grid, IID = grid)]
  covarData = cbind(covarData, gridData[, c(covarColnames, 'sex'), with = FALSE])
  gwasMetadata[, covarNum := paste0('1-', ncol(covarData) - ifelse(whichSex == 'both', 2, 3))]
  covarPath = tempfile('covar_', fileext = '.tsv')
  write_tsv(covarData, covarPath)

  phenoPlink = Reduce(function(...) merge(..., all = TRUE), phenoPlinkList)
  phenoPlink[is.na(phenoPlink)] = -9
  colnames(phenoPlink)[1] = 'FID'
  phenoPlink[, IID := FID]
  setcolorder(phenoPlink, c(1, ncol(phenoPlink)))
  phenoPath = tempfile('pheno_', fileext = '.tsv')
  write_tsv(phenoPlink, phenoPath)

  return(list(gwasMetadata, list(snp = snpPath, covar = covarPath, pheno = phenoPath)))}


makePlinkArgs = function(p, paths) {
  sprintf('%s --bfile %s --extract %s --covar %s --pheno %s --1 --memory %d --vif %d',
          '--logistic hide-covar beta --ci 0.95',
          file.path(p$dataPath, p$dataPrefix),paths$snp, paths$covar,
          paths$pheno, as.numeric(p$memSize), as.numeric(p$maxVif))}


runGwasPlink = function(resultDir, filePrefix, phecodeStr, covarNum, plinkArgs, execPath) {
  argsNow = sprintf('%s --pheno-name %s --covar-number %s --out %s_%s',
                    plinkArgs, phecodeStr, covarNum,
                    file.path(resultDir, filePrefix), phecodeStr)
  system2(execPath, argsNow)}


cleanPlinkOutput = function(resultDir, filePrefix, phecodeStr) {
  outputFilename = sprintf('%s_%s.assoc.logistic', filePrefix, phecodeStr)
  outputFilepath = file.path(resultDir, outputFilename)
  tmpFilepath = tempfile(paste0(phecodeStr, '_'))
  system(sprintf("cat %s | tr -s ' ' '\t' > %s", outputFilepath, tmpFilepath))
  system(sprintf("cat %s | sed 's/^[[:space:]]*//g' | sed 's/[[:space:]]*$//g' > %s",
                 tmpFilepath, outputFilepath))
  unlink(tmpFilepath)
  return(outputFilepath)}


# makeGlmStr = function(whichSex, nPC, splineDf) {
#   formStr = sprintf('status ~ genotype + rec_len + %s', paste0('last_age', 1:splineDf, collapse = ' + '))
#   if (whichSex == 'both') {
#     formStr = paste(formStr, '+ sex')}
#   if (nPC > 0) {
#     formStr = paste(formStr, '+', paste0('PC', 1:nPC, collapse = ' + '))}
#   formStr}
#
#
# runGlm = function(formulaStr, input) {
#   speedglm(formula(formulaStr), family = binomial(), data = input)}
