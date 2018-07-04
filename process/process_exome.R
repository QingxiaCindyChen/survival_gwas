library('data.table')
library('readr')
library('RODBC')
library('snpStats')
library('SNPRelate')
library('yaml')

con = odbcConnect('NZSQL', believeNRows = FALSE)

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) == 0) {
  paramFile = 'exome_params.yaml'
} else {
  paramFile = cmdArgs[1]}

paramDir = 'params'
params = read_yaml(file.path(paramDir, paramFile))
procParent = 'processed'
procDir = file.path(procParent, params$datasetName)

if (!file.exists(procDir)) {
  dir.create(procDir, recursive = TRUE)}
write_yaml(params, file.path(resultDir, 'params.yaml'))

############################################################

loadPhecodeIcdMapping = function(filepath) {
  d = read_csv(filepath, col_types = 'cc')
  setDT(d)
  setnames(d, old = 'icd9', new = 'icd')
  return(d)}

phecodeIcdMapping = loadPhecodeIcdMapping(file.path(procParent, 'phecode_icd9_map_unrolled.csv'))

############################################################

getPhenoRaw = function(con, gridTable, icds) {
  queryStr = sprintf("select distinct a.GRID, a.CODE as icd, a.ENTRY_DATE
                     from ICD_CODES a
                     inner join %s b
                     on a.GRID = b.GRID
                     where a.CODE in ('%s')
                     and b.EURO = 1
                     order by a.GRID, a.ENTRY_DATE;",
                     gridTable, paste(icds, collapse="','"))
  queryStr = gsub(pattern = '\\s+', replacement = ' ', x = queryStr)
  d = setDT(sqlQuery(con, queryStr, stringsAsFactors = FALSE, as.is = TRUE))
  colnames(d) = tolower(colnames(d))
  d[, entry_date := as.Date(entry_date)]
  return(d)}

phenoRaw = getPhenoRaw(con, params$gridTable, unique(phecodeIcdMapping$icd))

makePhenoData = function(phenoRaw, phecodeIcdMapping) {
  d = merge(phenoRaw, phecodeIcdMapping, by = 'icd', allow.cartesian = TRUE)
  d = unique(d[order(grid, phecode, entry_date), .(grid, phecode, entry_date)])
  return(d)}

phenoData = makePhenoData(phenoRaw, phecodeIcdMapping)
write_csv(phenoData, gzfile(file.path(procDir, 'phenotype_data.csv.gz')))

############################################################

makeGridData = function(con, gridTable) {
  queryStr = sprintf('select a.GRID, c.GENDER_EPIC as gender, c.DOB,
                     min(a.ENTRY_DATE) as FIRST_ENTRY_DATE, max(a.ENTRY_DATE) as LAST_ENTRY_DATE
                     from ICD_CODES a
                     inner join %s b
                     on a.GRID = b.GRID
                     inner join SD_RECORD c
                     on a.GRID = c.GRID
                     where b.EURO = 1
                     group by a.GRID, c.GENDER_EPIC, c.DOB;',
                     gridTable)
  queryStr = gsub(pattern = '\\s+', replacement = ' ', x = queryStr)
  gridData = setDT(sqlQuery(con, queryStr, stringsAsFactors = FALSE, as.is = TRUE))
  colnames(gridData) = tolower(colnames(gridData))
  gridData[, first_entry_date := as.Date(first_entry_date)]
  gridData[, last_entry_date := as.Date(last_entry_date)]
  return(gridData[order(grid)])}

gridData = makeGridData(con, params$gridTable)
write_csv(gridData, gzfile(file.path(procDir, 'grid_data.csv.gz')))

############################################################

makeGenoData = function(dataPathPrefix) {
  genoFull = read.plink(dataPathPrefix)
  genoSummary = col.summary(genoFull$genotypes)
  return(list(genoFull = genoFull, genoSummary = genoSummary))}

genoData = makeGenoData(params$plink$dataPathPrefix)
saveRDS(genoData, file.path(procDir, 'genotype_data.rds'))

############################################################

if (!is.null(params$geno$aimsFile)) {
  makePcData = function(genoSummary, gdsFile, paramDir, p) {
    aims = read_csv(file.path(paramDir, p$aimsFile),
                    col_names = FALSE, col_types = 'c')$X1

    idx = (genoSummary$MAF >= p$minMaf) &
      (genoSummary$Call.rate >= p$minCallRate) &
      (2 * pnorm(-abs(genoSummary$z.HWE)) >= p$minHwePval) &
      (rownames(genoSummary) %in% aims)

    pcRaw = snpgdsPCA(gdsFile, snp.id = rownames(genoSummary)[idx],
                      eigen.cnt = 10, num.thread = 8, algorithm = 'randomized')
    pcData = data.table(pcRaw$sample.id, pcRaw$eigenvect)
    colnames(pcData) = c('grid', paste0('PC', 1:ncol(pcRaw$eigenvect)))
    return(pcData)}

  snpgdsBED2GDS(paste0(params$plink$dataPathPrefix, '.bed'),
                paste0(params$plink$dataPathPrefix, '.fam'),
                paste0(params$plink$dataPathPrefix, '.bim'),
                paste0(params$plink$dataPathPrefix, '.gds'))
  gdsFile = snpgdsOpen(paste0(params$plink$dataPathPrefix, '.gds'))

  set.seed(3)
  pcData = makePcData(genoData$genoSummary, gdsFile, paramDir, params$geno)
  write_csv(pcData, gzfile(file.path(procDir, 'pc_data.csv.gz')))
}
