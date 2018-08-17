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

registerDoParallel(cores = params$slurm$cpusPerTask * params$slurm$doParFactor)

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

gwasMetadata = makeGwasMetadata(phecodeData, phenoData, phenoSummary)

phenoList = foreach(phenoIdx = 1:nrow(gwasMetadata), .combine = rbind) %dopar% {
  whichSex = gwasMetadata$whichSex[phenoIdx]
  phenoDataNow = phenoData[phecode == gwasMetadata$phecode[phenoIdx], .(grid, age)]
  inputBase = makeInput(phenoDataNow, gridData, whichSex,
                        params$pheno$minEvents, params$pheno$ageBuffer)
  phenoFilename = tempfile('pheno_', tmpdir = '', fileext = '.rds')
  saveRDS(inputBase, file.path(resultDir, phenoFilename), compress = FALSE)
  phenoPlink = makePhenoPlink(inputBase, gwasMetadata$phecodeStr[phenoIdx])
  list(phenoPlink, phenoFilename)}

phenoFilenames = unlist(phenoList[,2])

############################################################
# run cox regression

# TODO: make log file into tsv
# datetime, phecode, phenoIdx, nPhecodes, chunkIdx, nChunks
# first row: datetime, 'starting', NA, nPhecodes, NA, nChunks

chunkIdxUnique = unique(snpData$chunkIdx)
coxLog = createLogFile(resultDir, 'cox', length(chunkIdxUnique))

gwasChunkMetadata = foreach(chunkIdxNow = chunkIdxUnique, .combine = rbind) %dopar% {
  runGwasPhewasChunkCox(list(chunkIdx = chunkIdxNow),
                        snpData[chunkIdx == chunkIdxNow], genoData,
                        gwasMetadata, phenoFilenames, params, resultDir,
                        coxLog)}

gatherGwasChunks(gwasChunkMetadata, gwasMetadata, resultDir)
compressFiles(file.path(resultDir, gwasMetadata$coxFilename))
gwasMetadata[, coxFilename := paste0(coxFilename, '.gz')]

unlink(file.path(resultDir, phenoFilenames))
unlink(file.path(resultDir, gwasChunkMetadata$filename))
finishLogFile(coxLog)

############################################################
# prepare data for plink

plinkTmp = prepForPlink(snpData, gridData, covarColnames,
                        gwasMetadata, phenoList[,1])
gwasMetadata = plinkTmp[[1]]
plinkPaths = plinkTmp[[2]]
rm(plinkTmp)

############################################################
# run logistic regression in plink

plinkArgs = makePlinkArgs(params$plink, plinkPaths)
plinkLog = createLogFile(resultDir, 'logistic')

done = foreach(phenoIdx = 1:nrow(gwasMetadata)) %dopar% {
  # run plink
  runGwasPlink(resultDir, gwasMetadata$phecodeStr[phenoIdx],
               gwasMetadata$covarNum[phenoIdx], plinkArgs, params$plink$execPath)

  # fix plink's stupid output spacing
  outputFile = cleanPlinkOutput(resultDir, gwasMetadata$phecodeStr[phenoIdx])

  # rename and compress
  filepath = file.path(resultDir, gwasMetadata$logisticFilename[phenoIdx])
  file.rename(outputFile, filepath)
  system2('gzip', paste('-f', filepath))

  appendLogFile(plinkLog, gwasMetadata, phenoIdx)}

gwasMetadata[, logisticFilename := paste0(logisticFilename, '.gz')]
finishLogFile(plinkLog)

############################################################

d = c('snpData', 'genoData', 'phenoList', 'done', 'd')
save(list = setdiff(ls(), d), file = file.path(resultDir, 'workspace.Rdata'))
write_tsv(gwasMetadata, file.path(resultDir, 'gwas_metadata.tsv'))
