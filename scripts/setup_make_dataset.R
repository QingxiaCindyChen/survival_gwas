library('data.table')
library('readr')
library('RODBC')
library('snpStats')
library('SNPRelate')
library('yaml')

loadPhecodeIcdMapping = function(filepath) {
  d = read_csv(filepath, col_types = 'cc')
  setDT(d)
  setnames(d, old = 'icd9', new = 'icd')
  return(d)}


getPhenoRaw = function(con, icds, gridTable, euro = TRUE) {
  queryStr = sprintf("select distinct a.GRID, a.CODE as icd, a.ENTRY_DATE
                     from ICD_CODES a
                     inner join %s b
                     on a.GRID = b.GRID
                     where a.CODE in ('%s')
                     %s
                     order by a.GRID, a.ENTRY_DATE;",
                     gridTable, paste(icds, collapse="','"),
                     ifelse(is.null(euro) || euro, 'and b.EURO = 1', ''))
  queryStr = gsub(pattern = '\\s+', replacement = ' ', x = queryStr)
  d = setDT(sqlQuery(con, queryStr, stringsAsFactors = FALSE, as.is = TRUE))
  colnames(d) = tolower(colnames(d))
  d[, entry_date := as.Date(entry_date)]
  return(d)}


makePhenoData = function(phenoRaw, phecodeIcdMapping) {
  d = merge(phenoRaw, phecodeIcdMapping, by = 'icd', allow.cartesian = TRUE)
  d = unique(d[order(grid, phecode, entry_date), .(grid, phecode, entry_date)])
  return(d)}


makeGridData = function(con, gridTable, euro = TRUE) {
  queryStr = sprintf('select a.GRID, c.GENDER_EPIC as gender, c.DOB,
                     min(a.ENTRY_DATE) as FIRST_ENTRY_DATE,
                     max(a.ENTRY_DATE) as LAST_ENTRY_DATE
                     from ICD_CODES a
                     inner join %s b
                     on a.GRID = b.GRID
                     inner join SD_RECORD c
                     on a.GRID = c.GRID
                     %s
                     group by a.GRID, c.GENDER_EPIC, c.DOB;',
                     gridTable,
                     ifelse(is.null(euro) || euro, 'where b.EURO = 1', ''))
  queryStr = gsub(pattern = '\\s+', replacement = ' ', x = queryStr)
  gridData = setDT(sqlQuery(con, queryStr, stringsAsFactors = FALSE, as.is = TRUE))
  colnames(gridData) = tolower(colnames(gridData))
  gridData[, first_entry_date := as.Date(first_entry_date)]
  gridData[, last_entry_date := as.Date(last_entry_date)]
  return(gridData[order(grid)])}


makeMapData = function(plinkDataPathPrefix) {
  mapData = read_table2(paste0(plinkDataPathPrefix, '.bim'), na = c('0', '-9'),
                        col_names = FALSE, col_types = 'iccicc')
  colnames(mapData) = c('chr', 'snp', 'cM', 'pos', 'allele1', 'allele2')
  return(mapData)}


makePcData = function(plinkDataPathPrefix, paramsGeno, paramDir, nThreads) {
  gdsFile = snpgdsOpen(paste0(params$plink$dataPathPrefix, '.gds'))

  # currently not using paramsGeno$qc
  if (is.null(paramsGeno$aimsFile)) {
    snps = unlist(snpgdsLDpruning(gdsFile, ld.threshold = 0.2))
  } else {
    # not finished for exome
    # aims = read_csv(file.path(paramDir, paramsGeno$aimsFile),
    #                 col_names = FALSE, col_types = 'c')$X1
    #
    # idx = (genoSummary$MAF >= paramsGeno$minMaf) &
    #   (genoSummary$Call.rate >= paramsGeno$minCallRate) &
    #   (2 * pnorm(-abs(genoSummary$z.HWE)) >= paramsGeno$minHwePval) &
    #   (rownames(genoSummary) %in% aims)
    # snps =
  }
  pcRaw = snpgdsPCA(gdsFile, snp.id = snps, eigen.cnt = 10,
                    num.thread = nThreads, algorithm = 'randomized')

  pcData = data.table(pcRaw$sample.id, pcRaw$eigenvect)
  colnames(pcData) = c('grid', paste0('PC', 1:ncol(pcRaw$eigenvect)))
  return(list(pcRaw, pcData))}
