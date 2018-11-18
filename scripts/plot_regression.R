source(file.path('scripts', 'setup_regression.R'))

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) == 0) {
  resultDir = 'results/mega/20181030_104520'
} else {
  resultDir = cmdArgs[1]}

params = read_yaml(file.path(resultDir, 'params.yaml'))
procDir = file.path(procParent, params$datasetName)
plotDir = file.path(resultDir, 'plots')
dir.create(plotDir, recursive = TRUE)

registerDoParallel()
maxPvalLoad = 1e-2

############################################################

gwasMetadata = read_tsv(file.path(resultDir, 'gwas_metadata.tsv'),
                        col_types = 'cccccdc')
setDT(gwasMetadata)

mapData = read_csv(file.path(procDir, 'map_data.csv.gz'), col_types = 'iccicc')
setDT(mapData)

############################################################

gwasDataTmp = loadGwas(resultDir, gwasMetadata, maxPvalLoad)
gwasData = gwasDataTmp[[1]]
gwasNaData = gwasDataTmp[[2]]
gwasCorData = gwasDataTmp[[3]]
gwasLambdaData = gwasDataTmp[[4]]

gwasData[, pval := ifelse(is.na(pval), 1, pval)]
gwasData = mergeAll(gwasData, phecodeData, gwasMetadata, mapData)

rm(gwasDataTmp)

############################################################

maxPval = 1e-5
gwasDataSig = filterForSignificance(gwasData, maxPval)

############################################################

ptShp = 16
ptSz = 0.5
ptAlph = 0.2
lnSz = 0.5
lnCol = 'gray'

# correlations of p-values
pCor = ggplot(gwasCorData) +
  stat_ecdf(aes(x = r), pad = FALSE) +
  labs(x = expression(cor(p[Cox]*', '*p[logistic])),
       y = 'Cumulative fraction\nof phenotypes')

# genomic inflation factors
pLambda = plotLambda(gwasLambdaData, lnCol, lnSz, ptShp, 1, 0.5)

# p-values
resultTmp = plotPval(gwasDataSig, lnCol, lnSz, ptShp, ptSz, ptAlph)
gwasDataPval = resultTmp[[1]]
pPval = resultTmp[[2]]

pPvalZoom = pPval +
  scale_x_continuous(limits = c(5, 20)) +
  scale_y_continuous(limits = c(-5, 6))

# standard errors
resultTmp = plotSe(gwasDataSig)
gwasDataSe = resultTmp[[1]]
pSe = resultTmp[[2]]

# effect sizes
resultTmp = plotEffectSize(gwasDataSig, lnCol, lnSz, ptShp, ptSz, ptAlph)
gwasDataEffect = resultTmp[[1]]
pEffect = resultTmp[[2]]

cor(gwasDataEffect$logistic, gwasDataEffect$cox, use = 'na.or.complete')

p = plot_grid(pCor, pLambda, pEffect, pPval, pPvalZoom, pSe,
              labels = 'AUTO', ncol = 2, align = 'hv', axis = 'tb')
ggsave(file.path(plotDir, 'summary_stats.pdf'),
       plot = p, width = 6, height = 8)

############################################################

plotManhattan(gwasData, mapData, plotDir, width = 4, height = 4)

pList = plotManhattan(gwasData[phecode %in% c('185', '274.1')],
                      mapData, plotDir, save = FALSE)

p = plot_grid(plotlist = pList, nrow = 1, labels = 'AUTO')
ggsave(file.path(plotDir, 'example_manhattan.pdf'),
       plot = p, width = 8, height = 4.25)

############################################################

# a = merge(gwasDataPval, gwasDataSig, by = c('phecode', 'snp'))
#
# p = ggplot(a) +
#   geom_point(aes(x = log10(maf), y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph)
# print(p)
#
# p = ggplot(a) +
#   geom_point(aes(x = log10(nCases), y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
#   geom_smooth(aes(x = log10(nCases), y = cox - logistic), size = 0.5, method = 'loess', span = 0.5)
# print(p)
#
# p = ggplot(a) +
#   geom_point(aes(x = logRatio, y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
#   geom_smooth(aes(x = logRatio, y = cox - logistic), size = 0.5, method = 'loess', span = 0.5)
# print(p)
