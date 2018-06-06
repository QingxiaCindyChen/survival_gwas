############################################################
# where necessary, run logistic regression with speedglm

timeStarted = Sys.time()
cat(sprintf('%s started analysis (glm)\n', timeStarted), file = logFilepath, append = TRUE)

glmFilenames = foreach(ii = 1:nrow(phecodeDataKeep), .combine = c) %dopar% {
  logisticFilepath = file.path(resultDir, phecodeDataKeep$logisticFilename[ii])
  resultPlink = setDT(read_tsv(logisticFilepath, col_types = 'dcdccddddddd'))
  snpsNow = resultPlink[is.na(P), SNP]
  # return(length(snpsNow))

  if (length(snpsNow) == 0) {
    return(NA)
  } else {
    whichSex = phecodeDataKeep$whichSex[ii]
    glmStr = makeGlmStr(whichSex, nPC, splineDf)

    phenoNow = phenoGlm[, c('grid', phecodeDataKeep$phecodeStr[ii]), with = FALSE]
    colnames(phenoNow)[2] = 'status'
    inputBase = merge(gridData, phenoNow[!is.na(status)], by = 'grid')

    resultNow = foreach(snp = snpsNow, .combine = rbind) %do% {
      inputNow = addSnpToInput(inputBase, genoData$genoFull, snp)
      glmFit = runGlm(glmStr, inputNow)
      data.table(coef(summary(glmFit))[2, 1:3, drop = FALSE])}

    colnames(resultNow) = c('beta', 'se', 'z')
    resultNow[, pval := 2 * pnorm(-abs(z))] # speedglm encodes p-values as factors
    resultNow[, snp := snpsNow]

    glmFilename = sprintf('%s_%s_logistic_2.tsv.gz', filePrefix, phecodeDataKeep$phecodeStr[ii])
    write_tsv(resultNow, gzfile(file.path(resultDir, glmFilename)))

    cat(sprintf('%s completed phecode %s\n', Sys.time(), phecodeDataKeep$phecode[ii]),
        file = logFilepath, append = TRUE)
    glmFilename}}

phecodeDataKeep[, logisticFilename2 := glmFilenames]

timeElapsed = Sys.time() - timeStarted
cat(sprintf('Time elapsed of %.2f %s (glm)\n', timeElapsed, attr(timeElapsed, 'units')),
    file = logFilepath, append = TRUE)

