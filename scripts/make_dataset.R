source(file.path('scripts', 'setup_make_dataset.R'))

con = odbcConnect('NZSQL', believeNRows = FALSE)

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) == 0) {
  cmdArgs = 'params/mega/params_test1.yaml'}
paramDir = dirname(cmdArgs[1])
paramFile = basename(cmdArgs[1])

params = read_yaml(file.path(paramDir, paramFile))
procParent = 'processed'
procDir = file.path(procParent, params$datasetName)

if (!file.exists(procDir)) {
  dir.create(procDir, recursive = TRUE)}
write_yaml(params, file.path(procDir, 'params.yaml'))

############################################################

phecodeIcdMapping = loadPhecodeIcdMapping(file.path(procParent, 'phecode_icd9_map_unrolled.csv'))

############################################################

phenoRaw = getPhenoRaw(con, unique(phecodeIcdMapping$icd),
                       params$gridTable, params$gridTableEuro)
phenoData = makePhenoData(phenoRaw, phecodeIcdMapping)
write_csv(phenoData, gzfile(file.path(procDir, 'phenotype_data.csv.gz')))

############################################################

gridData = makeGridData(con, params$gridTable, params$gridTableEuro)
write_csv(gridData, gzfile(file.path(procDir, 'grid_data.csv.gz')))

############################################################

genoData = makeGenoData(params$plink$dataPathPrefix)
saveRDS(genoData, file.path(procDir, 'genotype_data.rds'))

############################################################

if (!is.null(params$geno$aimsFile)) {
  snpgdsBED2GDS(paste0(params$plink$dataPathPrefix, '.bed'),
                paste0(params$plink$dataPathPrefix, '.fam'),
                paste0(params$plink$dataPathPrefix, '.bim'),
                paste0(params$plink$dataPathPrefix, '.gds'))
  gdsFile = snpgdsOpen(paste0(params$plink$dataPathPrefix, '.gds'))

  set.seed(3)
  pcData = makePcData(genoData$genoSummary, gdsFile, paramDir, params$geno)
  write_csv(pcData, gzfile(file.path(procDir, 'pc_data.csv.gz')))}
