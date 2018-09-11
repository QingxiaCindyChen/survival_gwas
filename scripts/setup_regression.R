library('BEDMatrix')
library('cowplot')
library('data.table')
library('doParallel')
library('readr')
library('lubridate')
library('qqman')
library('survival')
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
  filenames = c(params$snpSubsetFile, params$phecodeSubsetFile,
                params$gwas$covarFile)
  for (filename in filenames) {
    if (!is.null(filename)) {
      file.copy(file.path(paramDir, filename), resultDir)}}
  write_yaml(params, file.path(resultDir, 'params.yaml'))}


writeSlurmRun = function(p, resultDir) {
  txt1 = c('#!/bin/bash',
           '#SBATCH --mail-user=%s',
           '#SBATCH --mail-type=ALL',
           '#SBATCH --constraint=[sandybridge|haswell|skylake]',
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
          '#SBATCH --constraint=[sandybridge|haswell|skylake]',
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

loadGenoData = function(plinkDataPathPrefix, editSampleNames = TRUE,
                        editVariantNames = TRUE) {
  genoData = BEDMatrix(plinkDataPathPrefix)
  if (editSampleNames) {
    rownames(genoData) = sapply(strsplit(rownames(genoData), '_'),
                                function(r) r[1])}
  if (editVariantNames) {
    colnames(genoData) = substr(colnames(genoData), 1,
                                nchar(colnames(genoData)) - 2)}
  return(genoData)}


loadSnpGenoData = function(p, plinkDataPathPrefix, snpSubsetPath = NULL) {
  warnOld = getOption('warn')
  options(warn = -1) # suppress harmless warning about NULL in read_tokens_
  snpsAll = read_table2(paste0(plinkDataPathPrefix, '.bim'), na = c('0', '-9'),
                        col_names = FALSE, col_types = cols_only(X2 = 'c'))$X2
  options(warn = warnOld)

  if (length(snpSubsetPath) == 0) {
    snpsIdx = 1:length(snpsAll)
  } else {
    snpsPref = unique(read_tsv(snpSubsetPath, col_names = FALSE, col_types = 'c')$X1)
    snpsIdx = match(snpsPref, snpsAll)
    snpsIdx = sort(snpsIdx[!is.na(snpsIdx)])}

  if (is.null(p$maxSnpsPerChunk)) {
    p$maxSnpsPerChunk = 1e4}
  chunkIdx = sort(rep_len(1:ceiling(length(snpsIdx) / as.integer(p$maxSnpsPerChunk)),
                          length(snpsIdx)))

  genoData = loadGenoData(plinkDataPathPrefix)

  if (p$qc) { # TODO: needs to be rewritten
    snpData = foreach(chunkIdxNow = unique(chunkIdx), .combine = rbind) %dopar% {
      snpsIdxNow = snpsIdx[chunkIdxNow == chunkIdx]
      # genoFull = read.plink(plinkDataPathPrefix, select.snps = snpsIdxNow)

      genoSummary = col.summary(genoFull$genotypes)
      idx = (genoSummary$MAF >= p$minMaf) &
        (genoSummary$Call.rate >= p$minCallRate) &
        (2 * pnorm(-abs(genoSummary$z.HWE)) >= p$minHwePval) &
        (genoFull$map$chromosome <= 22)

      data.table(snpName = genoFull$map$snp.name[idx],
                 snpIdx = snpsIdxNow[idx],
                 chunkIdx = chunkIdxNow)}
  } else {
    snpData = data.table(snpName = snpsAll[snpsIdx],
                         snpIdx = snpsIdx,
                         chunkIdx = chunkIdx)}

  return(list(snpData = snpData, genoData = genoData))}


loadGrid = function(procDir, plinkDataPathPrefix, minRecLen, paramsGwas, paramDir = NULL) {
  splineDf = paramsGwas$splineDf
  nPC = paramsGwas$nPC

  gridData = read_csv(file.path(procDir, 'grid_data.csv.gz'), col_types = 'ccDDD')
  setDT(gridData)
  gridData[, first_age := time_length(first_entry_date - dob, 'years')]
  gridData[, last_age := time_length(last_entry_date - dob, 'years')]
  gridData[, rec_len := last_age - first_age]
  gridData = gridData[first_age >= 0 & rec_len >= minRecLen]

  gridData = cbind(gridData, data.table(splines::ns(gridData$last_age, df = splineDf)))
  colIdx = (ncol(gridData) - splineDf + 1):ncol(gridData)
  colnames(gridData)[colIdx] = paste0('last_age', 1:splineDf)
  covarColnames = c('rec_len', colnames(gridData)[colIdx])

  if (nPC > 0) {
    pcData = read_csv(file.path(procDir, 'pc_data.csv.gz'), col_types = cols())
    setDT(pcData)
    gridData = merge(gridData, pcData[, 1:(1 + nPC)], by = 'grid')
    covarColnames = c(covarColnames, colnames(pcData)[2:(1 + nPC)])}

  if (!is.null(paramsGwas$covarFile)) {
    warnOld = getOption('warn')
    options(warn = -1) # suppress harmless warning about missing colunm name
    covarData = read_delim(file.path(paramDir, paramsGwas$covarFile),
                           delim = paramsGwas$covarFileDelim, col_types = cols())
    options(warn = warnOld)
    setDT(covarData)
    gridData = merge(gridData,
                     covarData[, c('IID', paramsGwas$covarsFromFile), with = FALSE],
                     by.x = 'grid', by.y = 'IID') # not checking for name clashes
    covarColnames = c(covarColnames, paramsGwas$covarsFromFile)}

  fam = read_table2(paste0(plinkDataPathPrefix, '.fam'),
                    col_names = FALSE, col_types = cols())
  setDT(fam)
  gridData = merge(gridData, fam[, .(grid = X1, sex = X5)], by = 'grid')
  return(list(gridData, covarColnames))}


loadPheno = function(procDir, p, gridData, phecodeSubsetPath) {
  phenoData = read_csv(file.path(procDir, 'phenotype_data.csv.gz'), col_types = 'ccD')
  setDT(phenoData)

  if (length(phecodeSubsetPath) > 0) {
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

createLogFile = function(resultDir, fileSuffix, nPhecodes, nChunks = NULL) {
  txt = c('datetime', 'phecode', 'phecodeIdx')
  d = data.table(datetimeStarted = Sys.time(), nPhecodes = nPhecodes)
  if (!is.null(nChunks)) {
    txt = c(txt, 'chunkIdx')
    d$nChunks = nChunks}

  path = file.path(resultDir, sprintf('progress_%s.tsv', fileSuffix))
  writeLines(paste(txt, collapse = '\t'), con = path)

  metaPath = file.path(resultDir, sprintf('progress_%s_meta.tsv', fileSuffix))
  write_tsv(d, metaPath)

  return(list(path = path, metaPath = metaPath,
              datetimeStarted = d$datetimeStarted,
              nPhecodes = nPhecodes, nChunks = nChunks))}


appendLogFile = function(logFile, gwasMetadata, phecodeIdx, chunkIdx = NULL) {
  if (is.null(logFile$nChunks)) {
    # completed an unchunked phecode
    d = data.table(datetime = Sys.time(),
                   phecode = gwasMetadata$phecode[phecodeIdx],
                   phecodeIdx = phecodeIdx)
  } else {
    # loaded genotypes for a chunk or completed a chunk for a phecode
    d = data.table(datetime = Sys.time(),
                   phecode = ifelse(phecodeIdx == 0, NA,
                                    gwasMetadata$phecode[phecodeIdx]),
                   phecodeIdx = phecodeIdx,
                   chunkIdx = chunkIdx)}
  return(write_tsv(d, logFile$path, append = TRUE))}


finishLogFile = function(logFile) {
  d = read_tsv(logFile$metaPath, col_types = cols())
  d$timeElapsed = Sys.time() - d$datetimeStarted
  d$timeElapsedUnits = attr(d$timeElapsed, 'units')
  invisible(write_tsv(d, logFile$metaPath))}


getProgress = function(resultDir) {
  progFilenames = c('progress_cox_meta.tsv', 'progress_cox.tsv',
                    'progress_logistic_meta.tsv', 'progress_logistic.tsv')
  names(progFilenames) = c('coxMeta', 'cox', 'logisticMeta', 'logistic')

  warnOld = getOption('warn')
  options(warn = -1) # suppress warning about parsers not matching colnames

  progList = foreach(progFilename = progFilenames) %do% {
    progFilepath = file.path(resultDir, progFilename)
    if (file.exists(progFilepath)) {
      d = setDT(read_tsv(progFilepath, col_types = cols(phecode = 'c')))
    } else {
      d = NA}}
  options(warn = warnOld)

  names(progList) = names(progFilenames)
  return(progList)}

############################################################
# functions for cox regression

makeGwasMetadata = function(phecodeData, phenoData, phenoSummary, p) {
  d = phecodeData[phecode %in% unique(phenoData$phecode), .(phecode, whichSex)]
  d[, phecodeStr := paste0('phe', gsub('.', 'p', phecode, fixed = TRUE))]
  if (p$cox) {
    d[, coxFilename := paste0(phecodeStr, '_cox.tsv')]
  } else {
    d[, coxFilename := NA]}
  if (p$logistic) {
    d[, logisticFilename := paste0(phecodeStr, '_logistic.tsv')]
  } else {
    d[, logisticFilename := NA]}
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


prepPhenoDataForGwas = function(resultDir, gwasMetadata, phenoData, gridData,
                                minEvents, ageBuffer) {
  nCores = getDoParWorkers()
  registerDoParallel(cores = 4)
  phenoList = foreach(phenoIdx = 1:nrow(gwasMetadata), .combine = rbind) %dopar% {
    whichSex = gwasMetadata$whichSex[phenoIdx]
    phenoDataNow = phenoData[phecode == gwasMetadata$phecode[phenoIdx], .(grid, age)]
    inputBase = makeInput(phenoDataNow, gridData, whichSex, minEvents, ageBuffer)
    phenoFilename = tempfile(sprintf('pheno_%s_', gwasMetadata$phecodeStr[phenoIdx]),
                             tmpdir = '', fileext = '.rds')
    saveRDS(inputBase, file.path(resultDir, phenoFilename), compress = FALSE)
    phenoPlink = makePhenoPlink(inputBase, gwasMetadata$phecodeStr[phenoIdx])
    list(phenoPlink, phenoFilename)}
  registerDoParallel(cores = nCores)
  return(phenoList)}


getColnamesKeep = function(whichSex, nPC, covarsFromFile) {
  colnamesKeep = c('grid', 'age1', 'age2', 'status')
  if (whichSex == 'both') {
    colnamesKeep = c(colnamesKeep, 'sex')}
  if (nPC > 0) {
    colnamesKeep = c(colnamesKeep, paste0('PC', 1:nPC))}
  colnamesKeep = c(colnamesKeep, covarsFromFile) # ok if covarsFromFile is null
  return(colnamesKeep)}


makeAgregInput = function(input, genoMat, snp) {
  input[, genotype := genoMat[grid, snp]]
  setcolorder(input, c('grid', 'age1', 'age2', 'status', 'genotype'))
  idx = !is.na(input$genotype)
  x = as.matrix(input[idx, 5:ncol(input)])
  y = with(input[idx], Surv(age1, age2, status))
  return(list(x = x, y = y))}


runAgreg = function(x, y, control) {
  fit = agreg.fit(x, y, strata = NULL, init = NULL, control = control,
                  method = 'efron', resid = FALSE, concordance = FALSE)
  return(fit)}


runGwasCox = function(inputBase, genoMat, whichSex, nPC, covarsFromFile) {
  colnamesKeep = getColnamesKeep(whichSex, nPC, covarsFromFile)
  inputKeep = inputBase[, colnamesKeep, with = FALSE]
  snps = colnames(genoMat)
  control = coxph.control()

  dVec = foreach(snp = snps, .combine = c) %do% {
    agInput = makeAgregInput(inputKeep, genoMat, snp)
    agFit = runAgreg(agInput$x, agInput$y, control)
    c(agFit$coefficients[1], sqrt(diag(agFit$var)[1]))}

  idx = seq(1, length(dVec), 2)
  d = data.table(beta = dVec[idx], se = dVec[idx + 1])
  d[, pval := 2 * pnorm(-abs(beta / se))]
  d[, snp := snps]
  return(d)}


runGwasPhewasChunkCox = function(byList, snpDataSubset, genoData, gwasMetadata,
                                 phenoFilenames, params, resultDir, coxLog) {
  genoMat = genoData[, snpDataSubset$snpIdx]
  appendLogFile(coxLog, gwasMetadata, 0, byList$chunkIdx)

  gwasChunkMetadata = foreach(phenoIdx = 1:nrow(gwasMetadata), .combine = rbind) %do% {
    whichSex = gwasMetadata$whichSex[phenoIdx]
    inputBase = readRDS(file.path(resultDir, phenoFilenames[phenoIdx]))

    gwasResult = runGwasCox(inputBase, genoMat, whichSex, params$gwas$nPC,
                            params$gwas$covarsFromFile)

    filePre = sprintf('%s_%.4d_', gwasMetadata$phecodeStr[phenoIdx], byList$chunkIdx)
    filename = tempfile(filePre, tmpdir = '', fileext = '.tsv')
    write_tsv(gwasResult, file.path(resultDir, filename), col_names = FALSE)
    appendLogFile(coxLog, gwasMetadata, phenoIdx, byList$chunkIdx)

    data.table(phecode = gwasMetadata$phecode[phenoIdx],
               chunkIdx = byList$chunkIdx,
               filename = filename)}
  return(gwasChunkMetadata)}


createGwasFiles = function(resultDir, filenames) {
  done = foreach(filename = filenames, .combine = c) %dopar% {
    con = file(file.path(resultDir, filename))
    writeLines('beta\tse\tpval\tsnp', con)
    close(con)}
  invisible(done)}


gatherGwasChunks = function(gwasChunkMetadata, gwasMetadata, resultDir) {
  createGwasFiles(resultDir, gwasMetadata$coxFilename)
  done = foreach(phenoIdx = 1:nrow(gwasMetadata), .combine = c) %dopar% {
    gcmNow = gwasChunkMetadata[phecode == gwasMetadata$phecode[phenoIdx]]
    catArgs = paste(file.path(resultDir, gcmNow$filename), '>>',
                    file.path(resultDir, gwasMetadata$coxFilename[phenoIdx]),
                    collapse = ' ')
    system(paste('cat', catArgs))}
  invisible(done)}


compressFiles = function(filepaths) {
  done = foreach(filepath = filepaths, .combine = c) %dopar% {
    system2('gzip', paste('-f', filepath))}
  invisible(done)}

############################################################
# functions for logistic regression using plink

makePhenoPlink = function(inputBase, phecodeStr) {
  phenoPlinkNow = inputBase[, .(grid, status)]
  colnames(phenoPlinkNow)[2] = phecodeStr
  setkeyv(phenoPlinkNow, 'grid')
  return(phenoPlinkNow)}


prepForPlink = function(snpData, gridData, covarColnames, gwasMetadata,
                        phenoPlinkList, resultDir) {
  snpFile = tempfile('snp_', tmpdir = '', fileext = '.tsv')
  write_tsv(snpData[, .(snpName)], file.path(resultDir, snpFile), col_names = FALSE)

  covarData = gridData[, .(FID = grid, IID = grid)]
  covarData = cbind(covarData, gridData[, c(covarColnames, 'sex'), with = FALSE])
  gwasMetadata[, covarNum := paste0('3-', ncol(covarData) - ifelse(whichSex == 'both', 0, 1))]
  covarFile = tempfile('covar_', tmpdir = '', fileext = '.tsv')
  write_tsv(covarData, file.path(resultDir, covarFile))

  phenoPlink = Reduce(function(...) merge(..., all = TRUE), phenoPlinkList)
  phenoPlink[is.na(phenoPlink)] = -9
  colnames(phenoPlink)[1] = 'FID'
  phenoPlink[, IID := FID]
  setcolorder(phenoPlink, c(1, ncol(phenoPlink)))
  phenoFile = tempfile('pheno_', tmpdir = '', fileext = '.tsv')
  write_tsv(phenoPlink, file.path(resultDir, phenoFile))

  return(list(gwasMetadata, list(snp = snpFile, covar = covarFile, pheno = phenoFile)))}


makePlinkArgs = function(p, paths) {
  sprintf('%s --bfile %s --extract %s --covar %s --pheno %s --memory %d --vif %d --max-corr %.10g',
          '--threads 1 --covar-variance-standardize --1 --glm hide-covar',
          p$dataPathPrefix, paths$snp, paths$covar, paths$pheno,
          as.numeric(p$memSize), as.numeric(p$maxVif), as.numeric(p$maxCorr))}


runGwasPlink = function(resultDir, phecodeStr, covarNum, plinkArgs, execPath) {
  argsNow = sprintf('%s --pheno-name %s --covar-col-nums %s --out %s',
                    plinkArgs, phecodeStr, covarNum,
                    file.path(resultDir, phecodeStr))
  system2(execPath, argsNow)}


cleanPlinkOutput = function(resultDir, phecodeStr, filenameNew, compress = TRUE) {
  filenameOrig = sprintf('%s.%s.glm.logistic', phecodeStr, phecodeStr)
  unlink(file.path(resultDir, paste0(filenameOrig, '.id')))
  file.rename(file.path(resultDir, filenameOrig),
              file.path(resultDir, filenameNew))
  if (compress) {
    system2('gzip', paste('-f', file.path(resultDir, filenameNew)))
    filenameNew = paste0(filenameNew, '.gz')}
  invisible(filenameNew)}

############################################################

loadGwasCox = function(path) {
  # col_type 'n' chokes on exponential notation
  d = setDT(read_tsv(path, col_types = 'dddc'))
  # d[, beta := -beta] # coefficients from cox are for the major allele
  return(d)}


loadGwasLogistic = function(path) {
  #CHROM	POS	ID	REF	ALT	A1	TEST	OBS_CT	OR	SE	Z_STAT	P
  d = setDT(read_tsv(path, col_types = 'ddcccccddddd'))
  d = d[, .(beta = log(OR), se = SE, pval = P, snp = ID)]
  return(d)}


loadGwasPerPhecode = function(resultDir, coxFilename, logisticFilename) {
  if (!is.na(coxFilename)) {
    coxFilepath = file.path(resultDir, coxFilename)
    dCox = loadGwasCox(coxFilepath)
    dCox[, method := 'cox']
  } else {
    dCox = data.table()}

  if (!is.na(logisticFilename)) {
    logisticFilepath = file.path(resultDir, logisticFilename)
    dLogistic = loadGwasLogistic(logisticFilepath)
    dLogistic[, method := 'logistic']
  } else {
    dLogistic = data.table()}
  return(rbind(dCox, dLogistic))}


loadGwas = function(resultDir, gwasMetadata, maxPvalLoad) {
  gdList = foreach(ii = 1:nrow(gwasMetadata), .combine = rbind) %dopar% {
    d = loadGwasPerPhecode(resultDir, gwasMetadata$coxFilename[ii],
                           gwasMetadata$logisticFilename[ii])
    d[, phecode := gwasMetadata$phecode[ii]]

    dLambda = d[, .(lambdaMed = median((beta / se)^2, na.rm = TRUE) /
                      qchisq(0.5, 1)), by = .(phecode, method)]

    d = d[, if (any(!is.na(pval)) && mean(log(pval), na.rm = TRUE) <= log(maxPvalLoad)) .SD,
          by = snp]
    list(d, dLambda)}

  gd = rbindlist(gdList[,1], use.names = TRUE)
  gdLambda = rbindlist(gdList[,2], use.names = TRUE)
  gdLambda = dcast(gdLambda, phecode ~ method, value.var = 'lambdaMed')
  return(list(gd, gdLambda))}


mergeAll = function(gwasData, phecodeData, gwasMetadata, mapData) {
  gwasData = merge(gwasData, phecodeData[, .(phecode, phenotype, group)],
                   by = 'phecode')
  gwasData = merge(gwasData, gwasMetadata[, .(phecode, phecodeStr, nCases)],
                   by = 'phecode')
  gwasData = merge(gwasData, mapData[, .(snp, chr, pos)], by = 'snp')
  # genoSummary = setDT(genoData$genoSummary, keep.rownames = TRUE)
  # genoSummary = genoSummary[, .(snp = rn, maf = MAF)]
  # gwasData = merge(gwasData, genoSummary, by = 'snp')
  return(gwasData)}


plotManhattan = function(byList, dtSubset, plotDir, main = NULL, ...) {
  filename = sprintf('%s_%s_man.pdf', byList$phecodeStr, byList$method)
  pdf(file.path(plotDir, filename), width = 6, height = 4)
  manhattan(dtSubset, p = 'pval', snp = 'snp', chr = 'chr', bp = 'pos',
            main = main, ...)
  dev.off()}


plotQq = function(byList, dtSubset, plotDir, main = NULL) {
  filename = sprintf('%s_%s_qq.pdf', byList$phecodeStr, byList$method)
  pdf(file.path(plotDir, filename), width = 6, height = 4)
  qq(dtSubset$pval, main = main)
  dev.off()}


plotManhattanAndQq = function(byList, dtSubset, plotDir, ...) {
  main = sprintf('%s (%s), %s regression', byList$phenotype,
                 byList$phecode, byList$method)
  dt = dtSubset[!is.na(pval)]
  plotManhattan(byList, dt, plotDir, main, ...)
  plotQq(byList, dt, plotDir, main)}


filterForSignificance = function(gwasData, maxPval) {
  d = gwasData[, if (mean(-log(pval), na.rm = TRUE) >= -log(maxPval)) .SD,
               by = .(phecode, snp)]
  d[, logRatio := log2(exp(beta))]
  d[, negLogPval := -log10(pval)]
  return(d)}


plotEffectSize = function(gwasData, lnCol, lnSz, ptShp, ptSz, ptAlph, md = TRUE) {
  d = dcast(gwasData, phecode + snp ~ method, value.var = 'logRatio')
  if (md) {
    p = ggplot(d) +
      geom_hline(yintercept = 0, color = lnCol, size = lnSz) +
      # geom_point(aes(x = (logistic - cox)/2, y = cox + logistic),
      geom_point(aes(x = (logistic + cox) / 2, y = logistic - cox),
                 shape = ptShp, size = ptSz, alpha = ptAlph) +
      labs(title = 'Effect size', x = 'Mean of hazard ratio and odds ratio',
           y = 'Odds ratio - hazard ratio')
  } else {
    pTmp = ggplot(d) +
      geom_abline(slope = 1, intercept = 0, color = lnCol, size = lnSz) +
      geom_point(aes(x = cox, y = logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
      labs(title = 'Effect size', x = 'log2(hazard ratio)', y = 'log2(odds ratio)')

    paramList = list(col = 'white', fill = 'darkgray', size = 0.25)
    p = ggExtra::ggMarginal(pTmp, type = 'histogram', binwidth = 0.15, boundary = 0,
                            xparams = paramList, yparams = paramList)}
  return(list(d, p))}


plotPval = function(gwasData, lnCol, lnSz, ptShp, ptSz, ptAlph) {
  d = dcast(gwasData, phecode + snp ~ method, value.var = 'negLogPval')
  p = ggplot(d) +
    geom_hline(yintercept = 0, color = lnCol, size = lnSz) +
    geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic),
               shape = ptShp, size = ptSz, alpha = ptAlph) +
    geom_smooth(aes(x = (logistic + cox) / 2, y = cox - logistic),
                size = 0.5, method = 'loess', span = 0.5) +
    labs(title = '-log10(p)')
  return(list(d, p))}


plotSe = function(gwasData, binwidth = 0.0005, limits = c(-0.015, 0.005)) {
  d = dcast(gwasData, phecode + snp ~ method, value.var = 'se')
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
