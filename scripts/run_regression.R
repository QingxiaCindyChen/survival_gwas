source(file.path('scripts', 'setup_regression.R'))

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) == 0) {
  cmdArgs = 'params/mega/params_test1.yaml'}
paramDir = dirname(cmdArgs[1])
paramFile = basename(cmdArgs[1])

params = read_yaml(file.path(paramDir, paramFile))
procDir = file.path(procParent, params$datasetName)

############################################################

if (Sys.getenv('SLURM_ARRAY_TASK_ID') == '') {
  resultDir = file.path(resultParent, params$datasetName,
                        format(Sys.time(), '%Y%m%d_%H%M%S'))
  writeResultDir(params, paramDir, resultDir)
} else {
  resultDir = paramDir}

registerDoParallel(cores = params$slurm$cpusPerTask * params$slurm$doparFactor)

############################################################
# load snp data

tmpData = loadSnpGenoData(params$geno, params$plink$dataPathPrefix,
                          file.path(paramDir, params$snpSubsetFile))
snpData = tmpData$snpData
genoData = tmpData$genoData
rm(tmpData)

############################################################
# load grid data

gridTmp = loadGrid(procDir, params$plink$dataPathPrefix,
                   params$pheno$minRecLen, params$gwas, paramDir)
gridData = gridTmp[[1]]
covarColnames = gridTmp[[2]]
rm(gridTmp)

############################################################
# load phenotype data

phenoTmp = loadPheno(procDir, params$pheno, gridData,
                     file.path(paramDir, params$phecodeSubsetFile))
phenoData = phenoTmp[[1]]
phenoSummary = phenoTmp[[2]]
rm(phenoTmp)

############################################################

gwasMetadata = makeGwasMetadata(phecodeData, phenoData, phenoSummary, params$gwas)
gwasMetadata[, nControls := nrow(gridData) - nCases - nSinglets] # since minCases > 1

if (params$gwas$cox || params$gwas$logistic) {
  phenoList = prepPhenoDataForGwas(resultDir, gwasMetadata, phenoData,
                                   gridData, params$pheno$minEvents,
                                   params$pheno$ageBuffer)
  phenoFilenames = unlist(phenoList[, 2])
}

############################################################
# run cox regression

if (params$gwas$cox) {
  chunkIdxUnique = unique(snpData$chunkIdx)
  coxLog = createLogFile(resultDir, 'cox', nrow(gwasMetadata),
                         length(chunkIdxUnique))

  gwasChunkMetadata = foreach(chunkIdxNow = chunkIdxUnique, .combine = rbind) %dopar% {
    runGwasPhewasChunkCox(list(chunkIdx = chunkIdxNow),
                          snpData[chunkIdx == chunkIdxNow], genoData, gwasMetadata,
                          phenoFilenames, params, resultDir, coxLog)}

  gatherGwasChunks(gwasChunkMetadata, gwasMetadata, resultDir)
  compressFiles(file.path(resultDir, gwasMetadata$coxFilename))
  gwasMetadata[, coxFilename := paste0(coxFilename, '.gz')]

  unlink(file.path(resultDir, gwasChunkMetadata$filename))
  finishLogFile(coxLog)
}

unlink(file.path(resultDir, phenoFilenames))

############################################################
# prepare data for plink

if (params$gwas$logistic) {
  plinkTmp = prepForPlink(snpData, gridData, covarColnames, gwasMetadata,
                          phenoList[,1], resultDir)
  gwasMetadata = plinkTmp[[1]]
  plinkFilenames = plinkTmp[[2]]
  rm(plinkTmp)
}

############################################################
# run logistic regression in plink

if (params$gwas$logistic) {
  plinkFilepaths = lapply(plinkFilenames, function(f) file.path(resultDir, f))
  plinkArgs = makePlinkArgs(params$plink, plinkFilepaths)
  plinkLog = createLogFile(resultDir, 'logistic', nrow(gwasMetadata))

  done = foreach(phenoIdx = 1:nrow(gwasMetadata)) %dopar% {
    runGwasPlink(resultDir, gwasMetadata$phecodeStr[phenoIdx],
                 gwasMetadata$covarNum[phenoIdx], plinkArgs,
                 params$plink$execPath)

    cleanPlinkOutput(resultDir, gwasMetadata$phecodeStr[phenoIdx],
                     gwasMetadata$logisticFilename[phenoIdx])

    appendLogFile(plinkLog, gwasMetadata, phenoIdx)}

  gwasMetadata[, logisticFilename := paste0(logisticFilename, '.gz')]

  unlink(file.path(resultDir, plinkFilenames))
  finishLogFile(plinkLog)
}

############################################################

d = c('snpData', 'genoData', 'phenoList', 'done', 'd')
save(list = setdiff(ls(), d), file = file.path(resultDir, 'workspace.Rdata'))
write_tsv(gwasMetadata, file.path(resultDir, 'gwas_metadata.tsv'))
