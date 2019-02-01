source(file.path('scripts', 'setup_regression.R'))

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) == 0) {
  cmdArgs = 'params/mega/params_test1.yaml'}
paramDir = dirname(cmdArgs[1])
paramFile = basename(cmdArgs[1])

params = read_yaml(file.path(paramDir, paramFile))
procDir = file.path(procParent, params$datasetName)
resultDir = file.path(resultParent, params$datasetName,
                      format(Sys.time(), '%Y%m%d_%H%M%S'))

############################################################

if (is.null(params$slurm$nTasks)) {
  params$slurm$nTasks = 1}

writeResultDir(params, paramDir, resultDir)
writeSlurmRun(params$slurm, resultDir)

if (params$slurm$nTasks > 1) {
  writeSlurmGather(params$slurm, resultDir)

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
  writeTaskDirs(gwasMetadata, params, paramDir, resultDir)}
