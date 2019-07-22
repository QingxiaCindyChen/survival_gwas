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
                        col_types = 'cccccdddc')
setDT(gwasMetadata)

mapData = read_csv(file.path(procDir, 'map_data.csv.gz'), col_types = 'iccicc')
setDT(mapData)

############################################################

gwasMetadataNice = merge(gwasMetadata[, .(phecode, whichSex, nCases, nControls)],
                         phecodeData[, .(phecode, phenotype)],
                         by = 'phecode')
setcolorder(gwasMetadataNice, 'phenotype')
write_csv(gwasMetadataNice, file.path(resultDir, 'gwas_metadata_nice.csv'))

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

gwasData[pval <= 5e-8, .N, by = method]

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
resultTmp = plotEffectSize(gwasDataSig, lnCol, lnSz, ptShp, ptSz, ptAlph,
                           md = FALSE)
gwasDataEffect = resultTmp[[1]]
pEffect = resultTmp[[2]]

cor(gwasDataEffect$logistic, gwasDataEffect$cox, use = 'na.or.complete')

p = plot_grid(pCor, pLambda, pEffect, pPval, pPvalZoom, pSe,
              labels = 'AUTO', ncol = 2, align = 'hv', axis = 'tb')
ggsave(file.path(plotDir, 'summary_stats.pdf'),
       plot = p, width = 6, height = 8)

############################################################

plotManhattan(gwasData, mapData, plotDir, width = 4, height = 4)

pList = plotManhattan(gwasData[phecode %in% c('165.1', '185', '250.2', '411.2')],
                      mapData, plotDir, save = FALSE)

p = plot_grid(plotlist = pList, nrow = 2, labels = 'AUTO')
ggsave(file.path(plotDir, 'example_manhattan.pdf'),
       plot = p, width = 8, height = 8.25)

############################################################

gdSelect = loadGwas(resultDir, gwasMetadata[phecode %in% c('165.1', '185', '250.2', '411.2')], 1)[[1]]
gdSelect[method == 'cox', method := 'Cox']
gdSelect = merge(gdSelect, phecodeData[, .(phecode, phenotype)], by = 'phecode')

transNegLog10 = scales::trans_new('neglog10', function(x) -log10(x), function(x) 10^(-x))
scaleBreaks = 10^(seq(-39, 0, 3))
scaleLabels = as.character(-log10(scaleBreaks))

phecodes = sort(unique(gdSelect$phecode))
pList = foreach(phecodeNow = phecodes) %dopar% {
  gdNow = gdSelect[phecode == phecodeNow]
  phenoLabel = sprintf('%s (%s)', gdNow$phenotype[1], gdNow$phecode[1])
  p = ggplot(gdNow, aes(sample = pval)) +
    facet_grid(. ~ method) +
    geom_abline(slope = 1, intercept = 0, color = 'gray') +
    stat_qq(size = 0.25, distribution = stats::qunif) +
    labs(title = phenoLabel,
         x = expression(-log[10](p)~expected),
         y = expression(-log[10](p)~observed)) +
    scale_x_continuous(trans = transNegLog10, breaks = scaleBreaks, labels = scaleLabels) +
    scale_y_continuous(trans = transNegLog10, breaks = scaleBreaks, labels = scaleLabels)
}
names(pList) = phecodes

p = plot_grid(plotlist = pList, align = 'hv', axis = 'tblr', ncol = 2, labels = 'AUTO')
ggsave(file.path(plotDir, 'example_qq.png'), plot = p, width = 7.5, height = 5)
