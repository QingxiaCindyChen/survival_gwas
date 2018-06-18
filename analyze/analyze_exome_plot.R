source(file.path('analyze', 'analyze_setup.R'))

analysisDir = 'exome_full'
filePrefix = 'exome'

############################################################
# ensures that plotting uses exactly the same data used in the analysis,
# including gridData, phenoData, snpData, and all parameters and functions

load(file.path(resultDir, analysisDir, sprintf('%s_workspace.Rdata', filePrefix)))
# resultDir = 'results'

genoData = readRDS(file.path(procDir, sprintf('%s_genotype_data.rds', filePrefix)))
snpData = merge(snpData, data.table(snp = rownames(genoData$genoSummary), genoData$genoSummary), by = 'snp')

# phecodeDataKeep = phecodeDataKeepOrig[phecode %in% c('274.1', '714.1', '335', '185', '427.21', '290.11')]

############################################################

registerDoParallel(cores = 2)

maxPvalLoad = 1e-3

resultList = foreach(ii = 1:nrow(phecodeDataKeep)) %dopar% {
  coxFilepath = file.path(resultDir, analysisDir, phecodeDataKeep$coxFilename[ii])
  dCox = setDT(read_tsv(coxFilepath, col_types = 'ddddc')) # col_type 'n' chokes on exponential notation
  dCox[, beta := -beta] # coefficients from coxph are for the major allele
  dCox[, method := 'cox']

  logisticFilepath = file.path(resultDir, analysisDir, phecodeDataKeep$logisticFilename[ii])
  dLogistic = setDT(read_tsv(logisticFilepath, col_types = 'dcdccddddddd'))
  colnames(dLogistic) = tolower(colnames(dLogistic))
  dLogistic = dLogistic[, .(beta, se, z = stat, pval = p, snp = snp, method = 'logistic')]

  # if (!is.na(phecodeDataKeep$logisticFilename2[ii])) {
  #   glmFilepath = file.path(resultDir, analysisDir, phecodeDataKeep$logisticFilename2[ii])
  #   dGlm = setDT(read_tsv(glmFilepath, col_types = 'ddddc'))
  #   dGlm[, beta := -beta]
  #   dGlm[, method:= 'logistic']
  #   dLogistic = rbind(dLogistic[!(snp %in% dGlm$snp)], dGlm)}

  d = rbind(dCox, dLogistic)
  d[, phecode := phecodeDataKeep$phecode[ii]]
  dLambda = d[, .(lambdaMed = median(z^2, na.rm = TRUE) / qchisq(0.5, 1)), by = .(phecode, method)]

  d = d[, if (mean(log(pval), na.rm = TRUE) <= log(maxPvalLoad)) .SD, by = snp]
  list(d, dLambda)}

result = rbindlist(lapply(resultList, function(d) d[[1]]), use.names = TRUE)
result = merge(result, snpData, by = 'snp', sort = FALSE)
result = merge(result, phecodeData[, .(phecode, phenotype, group, groupnum, color)], by = 'phecode')
result = merge(result, phecodeDataKeep[, .(phecode, nCases)], by = 'phecode')

resultLambda = rbindlist(lapply(resultList, function(d) d[[2]]), use.names = TRUE)
resultLambda = dcast(resultLambda, phecode ~ method, value.var = 'lambdaMed')

############################################################

# r = unique(result[, .(phecode, phenotype, method)])
#
# done = foreach(ii = 1:nrow(r)) %do% {
#   resultNow = merge(result, r[ii,], by = colnames(r))
#
#   plotTitle = sprintf('%s (%s), %s regression', r$phenotype[ii], r$phecode[ii], r$method[ii])
#   pheFileText = gsub('.', 'p', r$phecode[ii], fixed = TRUE)
#
#   filename = sprintf('%s_phe%s_%s_man.pdf', filePrefix, pheFileText, r$method[ii])
#   pdf(file.path(resultDir, analysisDir, filename), width = 6, height = 4)
#   manhattan(resultNow, p = 'pval', snp = 'snp', chr = 'chr', bp = 'pos', main = plotTitle)
#   dev.off()
#
#   filename = sprintf('%s_phe%s_%s_qq.pdf', filePrefix, pheFileText, r$method[ii])
#   pdf(file.path(resultDir, analysisDir, filename), width = 6, height = 4)
#   qq(resultNow$pval, main = plotTitle)
#   dev.off()
# }

############################################################

ptShp = 16
ptSz = 0.5
# ptShp = 21
# ptSz = 1
ptAlph = 0.2

lnSz = 0.5
lnCol = 'gray'

# hazard ratios are extremely similar to odds ratios
maxPval = 1e-5
resultSig = result[, if (mean(-log(pval), na.rm = TRUE) >= -log(maxPval)) .SD, by = .(phecode, snp)]
resultSig[, c('logRatio', 'negLogPval') := .(log2(exp(beta)), -log10(pval))]
resultSigEffect = dcast(resultSig, phecode + snp ~ method, value.var = 'logRatio')

p = ggplot(resultSigEffect) +
  geom_abline(slope = 1, intercept = 0, color = lnCol, size = lnSz) +
  geom_point(aes(x = logistic, y = cox), shape = 16, size = 0.5, alpha = ptAlph) +
  labs(title = 'Effect size', x = 'log2(hazard ratio)', y = 'log2(odds ratio)')

paramList = list(col = 'white', fill = 'darkgray', size = 0.25)
p1 = ggExtra::ggMarginal(p, type = 'histogram', binwidth = 0.15, boundary = 0,
                         xparams = paramList, yparams = paramList)

cor(resultSigEffect$logistic, resultSigEffect$cox, use = 'na.or.complete')

resultSigPval = dcast(resultSig, phecode + snp ~ method, value.var = 'negLogPval')
p2 = ggplot(resultSigPval) +
  geom_hline(yintercept = 0, color = lnCol, size = lnSz) +
  geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
  geom_smooth(aes(x = (logistic + cox) / 2, y = cox - logistic), size = 0.5, method = 'loess', span = 0.5) +
  labs(title = '-log10(p)')

# standard errors are slightly lower for cox
resultSigSe = dcast(resultSig, phecode + snp ~ method, value.var = 'se')
p3 = ggplot(resultSigSe) +
  geom_histogram(aes(x = cox - logistic), binwidth = 0.0005, boundary = 0,
                 size = 0.25, fill = 'darkgray', color = 'white') +
  labs(title = 'Standard error') +
  scale_x_continuous(limits = c(-0.015, 0.005))

# genomic inflation factors are similar or slightly higher
p4 = ggplot(resultLambda) +
  geom_hline(yintercept = 0, color = lnCol, size = lnSz) +
  geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
  labs(title = 'Lambda median') +
  scale_x_continuous(limits = c(0.8, 1.2)) +
  scale_y_continuous(limits = c(-0.5, 0.5))

p = plot_grid(p1, p2, p3, p4, align = 'hv', axis = 'tb', nrow = 2)
ggsave(file.path(resultDir, 'exome_full.pdf'), plot = p, width = 6, height = 5.75)

############################################################

resultSigPval = dcast(resultSig, phecode + snp + MAF + nCases + group + color ~ method,
                      value.var = 'negLogPval')

p = ggplot(resultSigPval) +
  geom_point(aes(x = log10(MAF), y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph)
print(p)

p = ggplot(resultSigPval) +
  geom_point(aes(x = log10(nCases), y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
  geom_smooth(aes(x = log10(nCases), y = cox - logistic), size = 0.5, method = 'loess', span = 0.5)
print(p)

p = ggplot(resultSigPval) +
  geom_point(aes(x = log10(nCases), y = cox - logistic, color = group),
             shape = ptShp, size = ptSz, alpha = ptAlph) +
  geom_smooth(aes(x = log10(nCases), y = cox - logistic), size = 0.5, method = 'loess', span = 0.5) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
print(p)



a = dcast(resultSig, phecode + snp ~ method, value.var = 'negLogPval')
a = merge(a, resultSig[method == 'cox', .(phecode, snp, logRatio)], by = c('phecode', 'snp'))

p = ggplot(a) +
  geom_point(aes(x = logRatio, y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
  geom_smooth(aes(x = logRatio, y = cox - logistic), size = 0.5, method = 'loess', span = 0.5)
print(p)


result1 = dcast(result, phecode + snp ~ method, value.var = 'pval')
result2 = result1[, .(r = cor(-log10(logistic), -log10(cox))), by = phecode]
p = ggplot(result2) +
  stat_ecdf(aes(x = r), pad = FALSE)
print(p)
