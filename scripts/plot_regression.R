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

# gwasData[, plotManhattan(.BY, .SD, plotDir, cex = 0.5),
#          by = .(phecode, phecodeStr, phenotype, method)]

# gwasData[, plotQq(.BY, .SD, plotDir),
#          by = .(phecode, phecodeStr, phenotype, method)]

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

sz = 0.25

mapDataPlot = mapData[order(chr, pos)]
mapDataPlot[, posIdx := .I]
mapDataPlot[, chrMod := (chr + 1) %% 2]
mapDataTicks = mapDataPlot[(chr %% 2) == 1, .SD[round(.N / 2)], by = chr]

phecodeNow = '185'
phenoLabel = phecodeData[phecode == phecodeNow,
                         .(label = sprintf('%s (%s)', phenotype, phecode))][[1]]
gdNow = gwasData[phecode == phecodeNow]
gdNow[, method := ifelse(method == 'cox', 'Cox', method)]
gdNow = merge(gdNow, mapDataPlot, by = c('chr', 'pos'))

p1 = ggplot(gdNow) +
  facet_grid(method ~ .) +
  geom_point(aes(x = posIdx, y = -log10(pval), color = factor(chrMod)),
             size = sz) +
  geom_hline(yintercept = -log10(5e-8), color = '#33a02c') +
  geom_hline(yintercept = -log10(1e-5), color = '#b2df8a') +
  labs(title = phenoLabel, x = 'Chromosome', y = expression(-log[10](p))) +
  scale_x_continuous(breaks = mapDataTicks$posIdx, labels = mapDataTicks$chr) +
  scale_color_brewer(palette = 'Paired') +
  guides(color = FALSE) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(size = 7))


phecodeNow = '274.1'
phenoLabel = phecodeData[phecode == phecodeNow,
                         .(label = sprintf('%s (%s)', phenotype, phecode))][[1]]
gdNow = gwasData[phecode == phecodeNow]
gdNow[, method := ifelse(method == 'cox', 'Cox', method)]
gdNow = merge(gdNow, mapDataPlot, by = c('chr', 'pos'))

p2 = ggplot(gdNow) +
  facet_grid(method ~ .) +
  geom_point(aes(x = posIdx, y = -log10(pval), color = factor(chrMod)),
             size = sz) +
  geom_hline(yintercept = -log10(5e-8), color = '#33a02c') +
  geom_hline(yintercept = -log10(1e-5), color = '#b2df8a') +
  labs(title = phenoLabel, x = 'Chromosome', y = expression(-log[10](p))) +
  scale_x_continuous(breaks = mapDataTicks$posIdx, labels = mapDataTicks$chr) +
  scale_color_brewer(palette = 'Paired') +
  guides(color = FALSE) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(size = 7))


p = plot_grid(p1, p2, nrow = 1, labels = 'AUTO')
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
