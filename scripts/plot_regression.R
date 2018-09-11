source(file.path('scripts', 'setup_regression.R'))

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) == 0) {
  resultDir = 'results/mega/20180905_093452'
} else {
  resultDir = cmdArgs[1]}

params = read_yaml(file.path(resultDir, 'params.yaml'))
procDir = file.path(procParent, params$datasetName)
plotDir = file.path(resultDir, 'plots')
dir.create(plotDir, recursive = TRUE)

registerDoParallel()
# maxPvalLoad = 1
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
gwasLambdaData = gwasDataTmp[[2]]

gwasData[, pval := ifelse(is.na(pval), 1, pval)]
gwasData = mergeAll(gwasData, phecodeData, gwasMetadata, mapData)

rm(gwasDataTmp)

############################################################

gwasData[, plotManhattanAndQq(.BY, .SD, plotDir, cex = 0.5),
         by = .(phecode, phecodeStr, phenotype, method)]

############################################################

maxPval = 1e-5
# maxPval = 1
gwasDataSig = filterForSignificance(gwasData, maxPval)

############################################################

ptShp = 16
ptSz = 0.5
ptAlph = 0.2
lnSz = 0.5
lnCol = 'gray'

# hazard ratios are extremely similar to odds ratios
resultTmp = plotEffectSize(gwasDataSig, lnCol, lnSz, ptShp, ptSz, ptAlph)
gwasDataEffect = resultTmp[[1]]
pEffect = resultTmp[[2]]

cor(gwasDataEffect$logistic, gwasDataEffect$cox, use = 'na.or.complete')

# p-values at the low end are smaller
resultTmp = plotPval(gwasDataSig, lnCol, lnSz, ptShp, ptSz, ptAlph)
gwasDataPval = resultTmp[[1]]
pPval = resultTmp[[2]]

# standard errors are slightly lower for cox
resultTmp = plotSe(gwasDataSig)
gwasDataSe = resultTmp[[1]]
pSe = resultTmp[[2]]

# genomic inflation factors are similar or slightly higher
pLambda = plotLambda(gwasLambdaData, lnCol, lnSz, ptShp, ptSz, ptAlph)

# p = plot_grid(p1, p2, p3, p4, align = 'hv', axis = 'tb', nrow = 2)
# ggsave(file.path(resultDir, 'exome_full.pdf'), plot = p, width = 6, height = 5.75)

############################################################

a = merge(gwasDataPval, gwasDataSig, by = c('phecode', 'snp'))

# p = ggplot(a) +
#   geom_point(aes(x = log10(maf), y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph)
# print(p)

p = ggplot(a) +
  geom_point(aes(x = log10(nCases), y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
  geom_smooth(aes(x = log10(nCases), y = cox - logistic), size = 0.5, method = 'loess', span = 0.5)
print(p)

p = ggplot(a) +
  geom_point(aes(x = logRatio, y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
  geom_smooth(aes(x = logRatio, y = cox - logistic), size = 0.5, method = 'loess', span = 0.5)
print(p)


a1 = dcast(gwasData, phecode + snp ~ method, value.var = 'pval')
a2 = a1[, .(r = cor(-log10(logistic), -log10(cox))), by = phecode]
p = ggplot(a2) +
  stat_ecdf(aes(x = r), pad = FALSE)
print(p)
