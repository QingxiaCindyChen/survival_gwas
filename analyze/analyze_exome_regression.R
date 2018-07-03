source(file.path('analyze', 'analyze_setup.R'))

cmdArgs = commandArgs(trailingOnly = TRUE)

paramDir = 'params'

if (length(cmdArgs) == 0) {
  paramFile = 'exome_params.yaml'
} else {
  paramFile = cmdArgs[1]}
params = read_yaml(file.path(paramDir, paramFile))

registerDoParallel(cores = params$nCores)


resultDir = paste0(params$filePrefix, format(Sys.time(), '_%Y%m%d_%H%M%S'))
if (!file.exists(resultDir)) {
  dir.create(resultDir)}

############################################################
# load snp data

genoKeep = loadGeno(procDir, params$filePrefix, params$geno,
                    file.path(paramDir, params$snpSubsetFile))

############################################################
# load grid data

gridTmp = loadGrid(procDir, params$filePrefix, params$pheno$minRecLen,
                   params$gwas$nPC, params$gwas$splineDf, genoKeep$fam)
gridData = gridTmp[[1]]
covarColnames = gridTmp[[2]]
rm(gridTmp)

############################################################
# load phenotype data

phenoTmp = loadPheno(procDir, params$filePrefix, params$pheno, gridData,
                     file.path(paramDir, params$phecodeSubsetFile))
phenoData = phenoTmp[[1]]
phenoSummary = phenoTmp[[2]]
rm(phenoTmp)

############################################################

gwasMetadata = makeGwasMetadata(params$filePrefix, phecodeData, phenoData, phenoSummary)

# split analysis across servers
# gwasMetadata = gwasMetadata[1:round(nrow(gwasMetadata) / 2)]
# gwasMetadata = gwasMetadata[(1 + round(nrow(gwasMetadata) / 2)):nrow(gwasMetadata)]
# change names of progress files and workspace file

############################################################
# run cox regression

coxLog = createLogFile(resultDir, params$filePrefix, 'cox')

phenoPlinkList = foreach(ii = 1:nrow(gwasMetadata)) %dopar% {
  whichSex = gwasMetadata$whichSex[ii]
  phenoDataNow = phenoData[phecode == gwasMetadata$phecode[ii], .(grid, age)]

  inputBase = makeInput(phenoDataNow, gridData, whichSex, params$pheno$minEvents,
                        params$pheno$ageBuffer)
  coxStr = makeCoxStr(whichSex, params$gwas$nPC)

  gwasResult = runGwasCox(inputBase, genoKeep, coxStr)
  write_tsv(gwasResult, gzfile(file.path(resultDir, gwasMetadata$coxFilename[ii])))

  appendLogFile(coxLog, gwasMetadata, ii)
  makePhenoPlink(inputBase, gwasMetadata$phecodeStr[ii])}

finishLogFile(coxLog)

############################################################
# prepare data for plink

plinkTmp = prepForPlink(colnames(genoKeep$genotypes), gridData, covarColnames,
                        gwasMetadata, phenoPlinkList)
gwasMetadata = plinkTmp[[1]]
plinkPaths = plinkTmp[[2]]
rm(plinkTmp)

############################################################
# run logistic regression in plink

plinkArgs = makePlinkArgs(params$plink, plinkPaths)
plinkLog = createLogFile(resultDir, params$filePrefix, 'logistic')

done = foreach(ii = 1:nrow(gwasMetadata)) %dopar% {
  # run plink
  runGwasPlink(resultDir, params$filePrefix, gwasMetadata$phecodeStr[ii],
               gwasMetadata$covarNum[ii], plinkArgs, params$plink$execPath)

  # fix plink's stupid output spacing
  outputFilepath = cleanPlinkOutput(resultDir, params$filePrefix, gwasMetadata$phecodeStr[ii])

  # rename and compress
  file.rename(outputFilepath, file.path(resultDir, gwasMetadata$logisticFilename[ii]))
  system2('gzip', paste('-f', file.path(resultDir, gwasMetadata$logisticFilename[ii])))

  appendLogFile(plinkLog, gwasMetadata, ii)}

gwasMetadata[, logisticFilename := paste0(logisticFilename, '.gz')]
finishLogFile(plinkLog)

############################################################

d = c('genoAll', 'genoKeep', 'phenoPlinkList', 'done', 'd')
save(list = setdiff(ls(), d),
     file = file.path(resultDir, sprintf('%s_workspace.Rdata', params$filePrefix)))
