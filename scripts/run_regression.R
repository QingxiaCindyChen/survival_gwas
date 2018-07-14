source(file.path('scripts', 'setup.R'))

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) == 0) {
  paramDir = 'params'
  paramFile = 'exome_params_test1.yaml'
} else {
  paramDir = dirname(cmdArgs[1])
  paramFile = basename(cmdArgs[1])}

params = read_yaml(file.path(paramDir, paramFile))
procDir = file.path(procParent, params$datasetName)
resultDir = file.path(resultParent, params$datasetName,
                      format(Sys.time(), '%Y%m%d_%H%M%S'))

if (!file.exists(resultDir)) {
  dir.create(resultDir, recursive = TRUE)}
write_yaml(params, file.path(resultDir, 'params.yaml'))

if (Sys.getenv('SLURM_CPUS_PER_TASK') != '') {
  registerDoParallel(cores = Sys.getenv('SLURM_CPUS_PER_TASK'))
} else {
  registerDoParallel(cores = params$nCores)
}

############################################################
# load snp data

genoKeep = loadGeno(procDir, params$geno, file.path(paramDir, params$snpSubsetFile))

############################################################
# load grid data

gridTmp = loadGrid(procDir, params$pheno$minRecLen, params$gwas$nPC,
                   params$gwas$splineDf, genoKeep$fam)
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

# split analysis across servers
# gwasMetadata = gwasMetadata[1:round(nrow(gwasMetadata) / 2)]
# gwasMetadata = gwasMetadata[(1 + round(nrow(gwasMetadata) / 2)):nrow(gwasMetadata)]
# change names of progress files and workspace file

############################################################
# run cox regression

coxLog = createLogFile(resultDir, 'cox')

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
plinkLog = createLogFile(resultDir, 'logistic')

done = foreach(ii = 1:nrow(gwasMetadata)) %dopar% {
  # run plink
  runGwasPlink(resultDir, gwasMetadata$phecodeStr[ii],
               gwasMetadata$covarNum[ii], plinkArgs, params$plink$execPath)

  # fix plink's stupid output spacing
  outputFile = cleanPlinkOutput(resultDir, gwasMetadata$phecodeStr[ii])

  # rename and compress
  file.rename(outputFile, file.path(resultDir, gwasMetadata$logisticFilename[ii]))
  system2('gzip', paste('-f', file.path(resultDir, gwasMetadata$logisticFilename[ii])))

  appendLogFile(plinkLog, gwasMetadata, ii)}

gwasMetadata[, logisticFilename := paste0(logisticFilename, '.gz')]
finishLogFile(plinkLog)

############################################################

write_tsv(gwasMetadata, file.path(resultDir, 'gwas_metadata.tsv'))

d = c('genoKeep', 'phenoPlinkList', 'done', 'd')
save(list = setdiff(ls(), d), file = file.path(resultDir, 'workspace.Rdata'))
