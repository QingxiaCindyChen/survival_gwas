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

studyPhecodeMap = read_csv(file.path(procDir, 'catalog_study_phecode.csv'),
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


# getLdProxy = function(snps, pop = 'EUR', r2 = TRUE, proxyDir = NULL, logPath = NULL) {
#   urlBase = 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?'
#   d = foreach(snp = snps, .combine = rbind) %dopar% {
#     urlFull = sprintf('%svar=%s&pop=%s&r2_d=%s',
#                       urlBase, snp, pop, ifelse(r2, 'r2', 'd'))
#     curlRaw = curl_fetch_memory(urlFull)
#     dNow = setDT(read_tsv(rawToChar(curlRaw$content)))
#     if (!is.null(proxyDir)) {
#       write_tsv(dNow, file.path(proxyDir, sprintf('%s.tsv', snp)))}
#     if (!is.null(logPath)) {
#       write_tsv(data.table(datetime = Sys.time(), rsid = snp), logPath,
#                 append = TRUE)}
#     dNow[, rsid := snp]
#     dNow}
#   setcolorder(d, 'rsid')
#   return(d)}


getLdProxy = function(snps, pop = 'EUR', r2 = TRUE, proxyDir = NULL,
                      logPath = NULL) {
  urlBase = 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?'
  d = foreach(snp = snps, .combine = rbind) %dopar% {
    urlFull = sprintf('%svar=%s&pop=%s&r2_d=%s',
                      urlBase, snp, pop, ifelse(r2, 'r2', 'd'))
    trying = TRUE
    while (trying) {
      tryCatch({
        curlRaw = curl_fetch_memory(urlFull)
        trying = FALSE
      }, error = function(e){}, finally = {})}
    # curlRaw = curl_fetch_memory(urlFull)

    dNow = setDT(read_tsv(rawToChar(curlRaw$content)))
    if (!is.null(proxyDir)) {
      write_tsv(dNow, file.path(proxyDir, sprintf('%s.tsv', snp)))}
    if (!is.null(logPath)) {
      write_tsv(data.table(datetime = Sys.time(), rsid = snp), logPath,
                append = TRUE)}
    dNow[, rsid := snp]
    dNow}
  setcolorder(d, 'rsid')
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

proxyDir = file.path(procDir, 'ldproxy')
dir.create(proxyDir, recursive = TRUE, showWarnings = FALSE)
logPath = file.path(procDir, 'ldproxy_progress.tsv')
writeLines(paste(c('datetime', 'rsid'), collapse = '\t'), con = logPath)

d = getLdProxy(unique(assocData$rsid), proxyDir = proxyDir, logPath = logPath)
# d = getLdProxy(unique(assocData$rsid)[1:20], proxyDir = proxyDir, logPath = logPath)

minLd = 0.8
d1 = d[R2 >= minLd, .(rsid, rsid2 = RS_Number, chrPos = Coord, r2 = R2)]
d1[, c('chrTmp', 'posTmp') := tstrsplit(chrPos, ':', fixed = TRUE)]
d1[, chr := as.integer(substr(chrTmp, 4, nchar(chrTmp)))]
d1[, pos := as.integer(posTmp)]
d1[, c('chrTmp', 'posTmp', 'chrPos') := NULL]

# maybe this part and below should be in a separate script,
# since it's specific to a platform (i.e., mega)
d2 = merge(d1, mapData[, .(snp, chr, pos)], by = c('chr', 'pos'))

d3 = merge(assocData[, .(phecode, rsid, ldBlock)],
           d2[, .(rsid, rsid2, snp)], by = 'rsid')
d3 = unique(d3[, .(phecode, ldBlock, rsid2, snp)])

write_tsv(d3, file.path(procParent, 'mega', 'catalog_assoc_data.tsv'))
