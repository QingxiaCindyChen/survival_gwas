library('data.table')
library('readr')
library('lubridate')
library('doParallel')
library('snpStats')
library('cowplot')
library('qqman')
library('yaml')

procParent = 'processed'
resultParent = 'results'

phecodeData = read_csv(file.path(procParent, 'phecode_data.csv.gz'),
                       col_types = 'ccc??????')
setDT(phecodeData)

theme_set(theme_light() +
            theme(axis.text = element_text(color = 'black'),
                  strip.text = element_text(color = 'black'),
                  panel.grid.minor = element_blank(),
                  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))

############################################################
# functions for splitting and gathering

writeTaskDirs = function(gwasMetadata, params, paramDir, resultDir) {
  gwasMetadata[, taskId := rep_len(1:params$slurm$nTasks, length.out = .N)]
  gwasMetadataList = split(gwasMetadata, by = 'taskId')
  params$phecodeSubsetFile = 'phecodes.tsv'

  for(ii in 1:length(gwasMetadataList)) {
    resultDirNow = sprintf('%s_%d', resultDir, ii)
    dir.create(resultDirNow, recursive = TRUE)
    write_tsv(gwasMetadataList[[ii]][, .(phecode)],
              file.path(resultDirNow, 'phecodes.tsv'), col_names = FALSE)
    if (!is.null(params$snpSubsetFile)) {
      file.copy(file.path(paramDir, params$snpSubsetFile), resultDirNow)}
    write_yaml(params, file.path(resultDirNow, 'params.yaml'))}}


writeResultDir = function(params, paramDir, resultDir) {
  dir.create(resultDir, recursive = TRUE)
  if (!is.null(params$snpSubsetFile)) {
    file.copy(file.path(paramDir, params$snpSubsetFile), resultDir)}
  if (!is.null(params$phecodeSubsetFile)) {
    file.copy(file.path(paramDir, params$phecodeSubsetFile), resultDir)}
  write_yaml(params, file.path(resultDir, 'params.yaml'))}


writeSlurmRun = function(p, resultDir) {
  txt1 = c('#!/bin/bash',
           '#SBATCH --mail-user=%s',
           '#SBATCH --mail-type=ALL',
           '#SBATCH --nodes=1',
           '#SBATCH --ntasks=1',
           '#SBATCH --cpus-per-task=%d',
           '#SBATCH --mem=%s',
           '#SBATCH --time=%s',
           '#SBATCH --array=1-%d',
           'module restore %s')
  if (p$nTasks <= 1) {
    txt2 = 'Rscript scripts/run_regression.R %s/params.yaml'
  } else {
    txt2 = 'Rscript scripts/run_regression.R %s_${SLURM_ARRAY_TASK_ID}/params.yaml'}

  txt = sprintf(paste0(c(txt1, txt2), collapse = '\n'),
                p$email, p$cpusPerTask, p$mem, p$time,
                p$nTasks, p$lmodCollection, resultDir)
  filename = sprintf('%s_run.slurm', basename(resultDir))
  con = file(file.path('scripts', filename))
  writeLines(txt, con)
  close(con)}


writeSlurmGather = function(p, resultDir) {
  txt = c('#!/bin/bash',
          '#SBATCH --mail-user=%s',
          '#SBATCH --mail-type=ALL',
          '#SBATCH --nodes=1',
          '#SBATCH --ntasks=1',
          '#SBATCH --cpus-per-task=4',
          '#SBATCH --mem=16G',
          '#SBATCH --time=01:00:00',
          'module restore %s',
          'Rscript scripts/gather_regression.R %s')
  txt = sprintf(paste0(txt, collapse = '\n'),
                p$email, p$lmodCollection, resultDir)
  filename = sprintf('%s_gather.slurm', basename(resultDir))
  con = file(file.path('scripts', filename))
  writeLines(txt, con)
  close(con)}


gatherTaskResults = function(taskDirs, taskFiles, resultDir, subDir) {
  gmList = foreach(taskDir = taskDirs) %do% {
    gwasFile = file.path(taskDir, 'gwas_metadata.tsv')
    if (file.exists(gwasFile)) {
      gm = setDT(read_tsv(gwasFile, col_types = 'cccccdc'))
      taskIdNow = strsplit(basename(taskDir), '_')[[1]][3]
      gm[, taskId := as.integer(taskIdNow)]

      file.copy(file.path(taskDir, gm$coxFilename), resultDir)
      file.copy(file.path(taskDir, gm$logisticFilename), resultDir)
      file.copy(file.path(taskDir, paste0(gm$phecodeStr, '.log')), resultDir)

      for (taskFile in taskFiles) {
        taskFileNew = sprintf('%s_%s.%s', tools::file_path_sans_ext(taskFile),
                              taskIdNow, tools::file_ext(taskFile))
        file.copy(file.path(taskDir, taskFile),
                  file.path(resultDir, subDir, taskFileNew))}
      gm
    } else {
      'incomplete'}}
  return(gmList)}

############################################################
# functions for loading data

filterGenoData = function(genoData, idx) {
  genoData$genoFull$genotypes = genoData$genoFull$genotypes[, idx]
  genoData$genoFull$map = genoData$genoFull$map[idx,]
  genoData$genoSummary = genoData$genoSummary[idx,]
  return(genoData)}


loadGeno = function(procDir, p, snpSubsetPath = NULL) {
  genoData = readRDS(file.path(procDir, 'genotype_data.rds'))

  idx = (genoData$genoSummary$MAF >= p$minMaf) &
    (genoData$genoSummary$Call.rate >= p$minCallRate) &
    (2 * pnorm(-abs(genoData$genoSummary$z.HWE)) >= p$minHwePval) &
    (genoData$genoFull$map$chromosome <= 22)
  genoData = filterGenoData(genoData, idx)

  if (!is.null(snpSubsetPath) & length(snpSubsetPath) > 0) {
    snps = unique(read_tsv(snpSubsetPath, col_types = 'c', col_names = FALSE)$X1)
    idx2 = colnames(genoData$genoFull$genotypes) %in% snps
    genoData = filterGenoData(genoData, idx2)}
  return(genoData)}


loadGrid = function(procDir, fam, minRecLen, p, paramDir = NULL) {
  splineDf = p$splineDf
  nPC = p$nPC

  gridData = read_csv(file.path(procDir, 'grid_data.csv.gz'), col_types = 'ccDDD')
  setDT(gridData)
  gridData[, first_age := time_length(first_entry_date - dob, 'years')]
  gridData[, last_age := time_length(last_entry_date - dob, 'years')]
  gridData[, rec_len := last_age - first_age]
  gridData = gridData[first_age >= 0 & rec_len >= minRecLen]

  gridData = cbind(gridData, data.table(splines::ns(gridData$last_age, df = splineDf)))
  colnames(gridData)[(ncol(gridData) - splineDf + 1):ncol(gridData)] = paste0('last_age', 1:splineDf)
  covarColnames = c('rec_len', paste0('last_age', 1:splineDf))

  if (nPC > 0) {
    pcData = read_csv(file.path(procDir, 'pc_data.csv.gz'), col_types = cols())
    setDT(pcData)
    gridData = merge(gridData, pcData[, 1:(1 + nPC)], by = 'grid')
    covarColnames = c(covarColnames, colnames(pcData)[2:(1 + nPC)])}

  if (!is.null(p$covarFile)) {
    covarData = read_delim(file.path(paramDir, p$covarFile),
                           delim = p$covarFileDelim, col_types = cols())
    setDT(covarData)
    gridData = merge(gridData,
                     covarData[, c('IID', p$covarsFromFile), with = FALSE],
                     by.x = 'grid', by.y = 'IID') # not checking for name clashes
    covarColnames = c(covarColnames, p$covarsFromFile)}

  genoTmp = data.table(fam)[, .(grid = pedigree, sex)]
  gridData = merge(gridData, genoTmp, by = 'grid')
  return(list(gridData, covarColnames))}


loadPheno = function(procDir, p, gridData, phecodeSubsetPath) {
  phenoData = read_csv(file.path(procDir, 'phenotype_data.csv.gz'), col_types = 'ccD')
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
# functions for log files

createLogFile = function(resultDir, fileSuffix) {
  path = file.path(resultDir, sprintf('progress_%s.txt', fileSuffix))
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
# functions for cox regression

makeGwasMetadata = function(phecodeData, phenoData, phenoSummary) {
  d = phecodeData[phecode %in% unique(phenoData$phecode), .(phecode, whichSex)]
  d[, phecodeStr := paste0('phe', gsub('.', 'p', phecode, fixed = TRUE))]
  d[, coxFilename := paste0(phecodeStr, '_cox.tsv.gz')]
  d[, logisticFilename := paste0(phecodeStr, '_logistic.tsv')]
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
  return(input)}


getColnamesKeep = function(whichSex, nPC, covarsFromFile) {
  colnamesKeep = c('grid', 'age1', 'age2', 'status')
  if (whichSex == 'both') {
    colnamesKeep = c(colnamesKeep, 'sex')}
  if (nPC > 0) {
    colnamesKeep = c(colnamesKeep, paste0('PC', 1:nPC))}
  colnamesKeep = c(colnamesKeep, covarsFromFile) # ok if covarsFromFile is null
  return(colnamesKeep)}


makeAgregInput = function(input, genoFull, snp) {
  input[, genotype := as(genoFull$genotypes[grid, snp], 'numeric')[,1]]
  setcolorder(input, c('grid', 'age1', 'age2', 'status', 'genotype'))
  idx = !is.na(input$genotype)
  x = as.matrix(input[idx, 5:ncol(input)])
  y = with(input[idx], Surv(age1, age2, status))
  return(list(x = x, y = y))}


runAgreg = function(x, y, control) {
  fit = agreg.fit(x, y, strata = NULL, init = NULL, control = control,
                  method = 'efron', resid = FALSE, concordance = FALSE)
  return(fit)}


runGwasCox = function(inputBase, genoFull, whichSex, nPC, covarsFromFile) {
  colnamesKeep = getColnamesKeep(whichSex, nPC, covarsFromFile)
  inputKeep = inputBase[, colnamesKeep, with = FALSE]
  snps = colnames(genoFull$genotypes)
  control = coxph.control()

  dVec = foreach(snp = snps, .combine = c) %do% {
    agInput = makeAgregInput(inputKeep, genoFull, snp)
    agFit = runAgreg(agInput$x, agInput$y, control)
    c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))}

  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
  d[, snp := snps]
  return(d)}

############################################################
# functions for logistic regression using plink

makePhenoPlink = function(inputBase, phecodeStr) {
  phenoPlinkNow = inputBase[, .(grid, status)]
  colnames(phenoPlinkNow)[2] = phecodeStr
  setkeyv(phenoPlinkNow, 'grid')
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
  sprintf('%s --bfile %s --extract %s --covar %s --pheno %s --memory %d --vif %d',
          '--1 --logistic hide-covar beta --ci 0.95',
          p$dataPathPrefix, paths$snp, paths$covar,
          paths$pheno, as.numeric(p$memSize), as.numeric(p$maxVif))}


runGwasPlink = function(resultDir, phecodeStr, covarNum, plinkArgs, execPath) {
  argsNow = sprintf('%s --pheno-name %s --covar-number %s --out %s',
                    plinkArgs, phecodeStr, covarNum,
                    file.path(resultDir, phecodeStr))
  system2(execPath, argsNow)}


cleanPlinkOutput = function(resultDir, phecodeStr) {
  outputFile = file.path(resultDir, paste0(phecodeStr, '.assoc.logistic'))
  tmpFile = tempfile(paste0(phecodeStr, '_'))
  system(sprintf("cat %s | tr -s ' ' '\t' > %s", outputFile, tmpFile))
  system(sprintf("cat %s | sed 's/^[[:space:]]*//g' | sed 's/[[:space:]]*$//g' > %s",
                 tmpFile, outputFile))
  unlink(tmpFile)
  return(outputFile)}

############################################################

loadGwasCox = function(path) {
  # col_type 'n' chokes on exponential notation
  d = setDT(read_tsv(path, col_types = 'dddc'))
  d[, beta := -beta] # coefficients from cox are for the major allele
  return(d)}


loadGwasLogistic = function(path) {
  d = setDT(read_tsv(path, col_types = 'dcdccddddddd'))
  colnames(d) = tolower(colnames(d))
  d = d[, .(beta, se, pval = p, snp = snp)]
  return(d)}


loadGwasPerPhecode = function(resultDir, coxFilename, logisticFilename) {
  coxFilepath = file.path(resultDir, coxFilename)
  dCox = loadGwasCox(coxFilepath)
  dCox[, method := 'cox']

  logisticFilepath = file.path(resultDir, logisticFilename)
  dLogistic = loadGwasLogistic(logisticFilepath)
  dLogistic[, method := 'logistic']
  return(rbind(dCox, dLogistic))}


loadGwas = function(resultDir, gwasMetadata, maxPvalLoad) {
  gdList = foreach(ii = 1:nrow(gwasMetadata), .combine = rbind) %dopar% {
    d = loadGwasPerPhecode(resultDir, gwasMetadata$coxFilename[ii],
                           gwasMetadata$logisticFilename[ii])
    d[, phecode := gwasMetadata$phecode[ii]]

    dLambda = d[, .(lambdaMed = median((beta / se)^2, na.rm = TRUE) /
                      qchisq(0.5, 1)), by = .(phecode, method)]

    d = d[, if (mean(log(pval), na.rm = TRUE) <= log(maxPvalLoad)) .SD, by = snp]
    list(d, dLambda)}

  gd = rbindlist(gdList[,1], use.names = TRUE)
  gdLambda = rbindlist(gdList[,2], use.names = TRUE)
  gdLambda = dcast(gdLambda, phecode ~ method, value.var = 'lambdaMed')
  return(list(gd, gdLambda))}


mergeAll = function(gwasData, phecodeData, gwasMetadata, genoData) {
  mapData = setDT(genoData$genoFull$map)
  mapData = mapData[, .(snp = snp.name, chr = chromosome, pos = position)]
  genoSummary = setDT(genoData$genoSummary, keep.rownames = TRUE)
  genoSummary = genoSummary[, .(snp = rn, maf = MAF)]
  gData = merge(gwasData, phecodeData[, .(phecode, phenotype, group)],
                by = 'phecode')
  gData = merge(gData, gwasMetadata[, .(phecode, nCases)], by = 'phecode')
  gData = merge(gData, mapData, by = 'snp')
  gData = merge(gData, genoSummary, by = 'snp')
  return(gData)}


plotManhattan = function(byList, dtSubset, plotDir, main = NULL) {
  filename = sprintf('%s_%s_man.pdf', byList$phecodeStr, byList$method)
  pdf(file.path(plotDir, filename), width = 6, height = 4)
  manhattan(dtSubset, p = 'pval', snp = 'snp', chr = 'chr', bp = 'pos', main = main)
  dev.off()}


plotQq = function(byList, dtSubset, plotDir, main = NULL) {
  filename = sprintf('%s_%s_qq.pdf', byList$phecodeStr, byList$method)
  pdf(file.path(plotDir, filename), width = 6, height = 4)
  qq(dtSubset$pval, main = main)
  dev.off()}


plotManhattanAndQq = function(byList, dtSubset, plotDir) {
  main = sprintf('%s (%s), %s regression', byList$phenotype,
                 byList$phecode, byList$method)
  plotManhattan(byList, dtSubset, plotDir, main)
  plotQq(byList, dtSubset, plotDir, main)}


filterForSignificance = function(gData, maxPval) {
  d = gData[, if (mean(-log(pval), na.rm = TRUE) >= -log(maxPval)) .SD,
            by = .(phecode, snp)]
  d[, logRatio := log2(exp(beta))]
  d[, negLogPval := -log10(pval)]
  return(d)}


plotEffectSize = function(gData, lnCol, lnSz, ptShp, ptSz, ptAlph) {
  d = dcast(gData, phecode + snp ~ method, value.var = 'logRatio')

  pTmp = ggplot(d) +
    geom_abline(slope = 1, intercept = 0, color = lnCol, size = lnSz) +
    geom_point(aes(x = logistic, y = cox), shape = ptShp, size = ptSz, alpha = ptAlph) +
    labs(title = 'Effect size', x = 'log2(hazard ratio)', y = 'log2(odds ratio)')

  paramList = list(col = 'white', fill = 'darkgray', size = 0.25)
  p = ggExtra::ggMarginal(pTmp, type = 'histogram', binwidth = 0.15, boundary = 0,
                          xparams = paramList, yparams = paramList)
  return(list(d, p))}


plotPval = function(gData, lnCol, lnSz, ptShp, ptSz, ptAlph) {
  d = dcast(gData, phecode + snp ~ method, value.var = 'negLogPval')
  p = ggplot(d) +
    geom_hline(yintercept = 0, color = lnCol, size = lnSz) +
    geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic),
               shape = ptShp, size = ptSz, alpha = ptAlph) +
    geom_smooth(aes(x = (logistic + cox) / 2, y = cox - logistic),
                size = 0.5, method = 'loess', span = 0.5) +
    labs(title = '-log10(p)')
  return(list(d, p))}


plotSe = function(gData, binwidth = 0.0005, limits = c(-0.015, 0.005)) {
  d = dcast(gData, phecode + snp ~ method, value.var = 'se')
  p = ggplot(d) +
    geom_histogram(aes(x = cox - logistic), binwidth = binwidth, boundary = 0,
                   size = 0.25, fill = 'darkgray', color = 'white') +
    labs(title = 'Standard error') +
    scale_x_continuous(limits = limits)
  return(list(d, p))}


plotLambda = function(gwasLambdaData, lnCol, lnSz, ptShp, ptSz, ptAlph,
                      xlims = NULL, ylims = NULL) {
  p = ggplot(gwasLambdaData) +
    geom_hline(yintercept = 0, color = lnCol, size = lnSz) +
    geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic),
               shape = ptShp, size = ptSz, alpha = ptAlph) +
    labs(title = 'Lambda median') +
    scale_x_continuous(limits = xlims) +
    scale_y_continuous(limits = ylims)
  return(p)}
