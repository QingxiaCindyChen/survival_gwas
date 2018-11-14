library('curl')
library('data.table')
library('doParallel')
library('readr')

registerDoParallel(cores = 4)

procParent = 'processed'
procDir = file.path(procParent, 'gwas_catalog')
assocFilename = 'gwas_catalog_v1.0.2-associations_e93_r2018-08-28.tsv'

assocDataRaw = suppressWarnings(read_tsv(file.path(procDir, assocFilename)))
setDT(assocDataRaw)
colnames(assocDataRaw) = gsub('/|\\-|\\s', '_', tolower(colnames(assocDataRaw)))
colnames(assocDataRaw) = gsub('\\[|\\]|\\(|\\)', '', tolower(colnames(assocDataRaw)))

studyPhecodeMap = read_csv(file.path('params', 'catalog_study_phecode.csv'),
                           col_types = 'cc')
setDT(studyPhecodeMap)

mapData = read_csv(file.path(procParent, 'mega', 'map_data.csv.gz'),
                   col_types = 'iccicc')
setDT(mapData)

############################################################

getLdMat = function(snps, pop = 'EUR', r2 = TRUE) {
  urlBase = 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldmatrix?'
  snpsConcat = paste(snps, collapse = '%0A')
  popConcat = paste(pop, collapse = '%2B')
  urlFull = sprintf('%ssnps=%s&pop=%s&r2_d=%s', urlBase, snpsConcat, popConcat,
                    ifelse(r2, 'r2', 'd'))

  curlRaw = curl_fetch_memory(urlFull)
  ldDf = read_tsv(rawToChar(curlRaw$content),
                  col_types = cols(RS_number = 'c', .default = 'd'))
  ldMat = as.matrix(ldDf[, 2:ncol(ldDf)])
  ldMat[is.na(ldMat)] = 0 # monoallelic SNPs get returned as NA
  rownames(ldMat) = colnames(ldMat)

  snpsMissing = setdiff(snps, colnames(ldMat))
  for (snp in snpsMissing) {
    ldMat = cbind(ldMat, rep(0, nrow(ldMat)))
    ldMat = rbind(ldMat, matrix(c(rep(0, ncol(ldMat) - 1), 1),
                                nrow = 1, dimnames = list(snp)))}
  colnames(ldMat) = rownames(ldMat)
  ldMat = ldMat[snps, snps]
  return(ldMat)}


getLdBlocksFromLdMat = function(ldMat, minLd) {
  snps = colnames(ldMat)
  ldBlocks = snps
  diag(ldMat) = 0
  while (any(ldMat >= minLd, na.rm = TRUE)) {
    idx = arrayInd(which.max(ldMat), dim(ldMat))
    ldBlocks[idx[1]] = ldBlocks[idx[2]]
    ldMat[idx[1], idx[2]] = 0
    ldMat[idx[2], idx[1]] = 0}
  return(ldBlocks)}


getLdBlocks = function(snps, pop = 'EUR', r2 = TRUE, minLd = 0.8) {
  idx = !is.na(snps)
  if (length(snps) == 1 | sum(idx) == 1) {
    return(snps)}
  ldBlocks = rep(NA, length(snps))
  ldMat = getLdMat(snps[idx], pop, r2)
  ldBlocks[idx] = getLdBlocksFromLdMat(ldMat, minLd)
  return(ldBlocks)}


getSnpInfo = function(snps) {
  urlBase = 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldhap?'
  snpInfo = data.table(rsid = as.character(rep_len(NA, length(snps))),
                       chrPos = as.character(rep_len(NA, length(snps))))

  for (ii in 1:length(snps)) {
    urlFull = sprintf('%ssnps=%s&pop=EUR', urlBase, snps[ii])
    curlRaw = curl_fetch_memory(urlFull)
    if (startsWith(rawToChar(curlRaw$content), 'RS_Number')) {
      d = suppressWarnings(read_tsv(rawToChar(curlRaw$content), n_max = 1))
      if (!(startsWith(d[[1]][1], '#') || startsWith(d[[1]][1], '.'))) {
        set(snpInfo, i = ii, j = 1:2, as.list(d[, 1:2]))}}}
  return(snpInfo)}


getLdProxyRaw = function(snps, proxyDir, pop = 'EUR', r2 = TRUE,
                         logPath = NULL) {
  urlBase = 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?'
  d = foreach(snp = snps) %dopar% {
    urlFull = sprintf('%svar=%s&pop=%s&r2_d=%s',
                      urlBase, snp, pop, ifelse(r2, 'r2', 'd'))

    ldFilepath = file.path(proxyDir, sprintf('%s.tsv', snp))
    if (file.exists(ldFilepath)) {
      dNow = setDT(read_tsv(ldFilepath))
    } else {
      trying = TRUE
      while (trying) {
        tryCatch({
          curlRaw = curl_fetch_memory(urlFull)
          trying = FALSE
        }, error = function(e){}, finally = {})}
      dNow = setDT(read_tsv(rawToChar(curlRaw$content)))
      write_tsv(dNow, ldFilepath)}

    if (!is.null(logPath)) {
      write_tsv(data.table(datetime = Sys.time(), rsid = snp), logPath,
                append = TRUE)}
    dNow}
  names(d) = snps
  return(d)}


getLdProxy = function(dList, assocDataUnique, minLd = 0.8) {
  d = foreach(dRaw = dList, rsidRaw = names(dList), .combine = rbind) %dopar% {
    if (ncol(dRaw) == 10) {
      dNow = dRaw[R2 >= minLd, .(rsid = rsidRaw, rsid2 = RS_Number, chrPos = Coord, r2 = R2)]
      dNow[, c('chrTmp', 'posTmp') := tstrsplit(chrPos, ':', fixed = TRUE)]
      dNow[, chr := as.integer(substr(chrTmp, 4, nchar(chrTmp)))]
      dNow[, pos := as.integer(posTmp)]
      dNow[, c('chrTmp', 'posTmp', 'chrPos') := NULL]
      setcolorder(dNow, c('rsid', 'rsid2', 'chr', 'pos', 'r2'))
    } else {
      dNow = assocDataUnique[rsid == rsidRaw, .(rsid, rsid2 = rsidRaw, chr, pos, r2 = 1)]}
    dNow}
  return(d)}

############################################################
# get GRCh37 info for catalog snps

maxPval = 5e-8

assocData = assocDataRaw[(p_value <= maxPval) &
                           !(chr_id %in% c('X', 'Y')) &
                           !is.na(chr_pos)]
setnames(assocData, 'snps', 'snp')

assocData = merge(assocData, studyPhecodeMap, by = 'study_accession')
assocData = unique(assocData[, .(snp, phecode)])
setorderv(assocData, c('phecode', 'snp'))
setcolorder(assocData, c('phecode', 'snp'))

assocDataUnique = unique(assocData[, .(snp)])
assocDataUnique[, c('rsid', 'chrPos') := getSnpInfo(snp)]
assocDataUnique[, c('chrTmp', 'posTmp') := tstrsplit(chrPos, ':', fixed = TRUE)]
assocDataUnique[, chr := as.integer(substr(chrTmp, 4, nchar(chrTmp)))]
assocDataUnique[, pos := as.integer(posTmp)]
assocDataUnique[, c('chrPos', 'chrTmp', 'posTmp') := NULL]

############################################################
# consolidate snps for each phecode into blocks based on LD

assocData = merge(assocData[, .(snp, phecode)], assocDataUnique[!is.na(rsid)],
                  by = 'snp')
assocData[, ldBlock := getLdBlocks(rsid), by = .(phecode, chr)]

assocData[, .(nAssoc = uniqueN(rsid)), by = phecode]
assocData[, .(nAssoc = uniqueN(ldBlock)), by = phecode]

############################################################
# get the proxy snps

minLd = 0.8

proxyDir = file.path(procDir, 'ldproxy')
dir.create(proxyDir, recursive = TRUE, showWarnings = FALSE)
logPath = file.path(procDir, 'ldproxy_progress.tsv')
writeLines(paste(c('datetime', 'rsid'), collapse = '\t'), con = logPath)

ldProxyList = getLdProxyRaw(unique(assocData$rsid), proxyDir = proxyDir,
                            logPath = logPath)

ldProxyData = getLdProxy(ldProxyList, assocDataUnique, minLd)

# maybe the lines below should be in a separate script,
# since they're specific to a platform (i.e., mega)
ldProxyMega = merge(ldProxyData, mapData[, .(snp, chr, pos)], by = c('chr', 'pos'))

expandedAssocData = merge(assocData[, .(phecode, rsid, ldBlock)],
                          ldProxyMega[, .(rsid, rsid2, snp)], by = 'rsid')
expandedAssocData = unique(expandedAssocData[, .(phecode, ldBlock, rsid2, snp)])

write_tsv(expandedAssocData, file.path(procParent, 'mega', 'catalog_assoc_data.tsv'))
