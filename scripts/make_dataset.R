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

mapData = makeMapData(params$plink$dataPathPrefix)
write_csv(mapData, gzfile(file.path(procDir, 'map_data.csv.gz')))

############################################################

snpgdsBED2GDS(paste0(params$plink$dataPathPrefix, '.bed'),
              paste0(params$plink$dataPathPrefix, '.fam'),
              paste0(params$plink$dataPathPrefix, '.bim'),
              paste0(params$plink$dataPathPrefix, '.gds'))

set.seed(3)
pcDataTmp = makePcData(params$plink$dataPathPrefix, params$geno, paramDir,
                       params$slurm$cpusPerTask * params$slurm$doparFactor)
pcRaw = pcDataTmp[[1]]
pcData = pcDataTmp[[2]]

saveRDS(pcRaw, file.path(procDir, 'pc_data_raw.rds'))
write_csv(pcData, gzfile(file.path(procDir, 'pc_data.csv.gz')))
