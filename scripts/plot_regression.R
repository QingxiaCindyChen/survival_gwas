source(file.path('scripts', 'setup_regression.R'))

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) == 0) {
  resultDir = 'results/exome/20180720_110952'
} else {
  resultDir = cmdArgs[1]}

params = read_yaml(file.path(resultDir, 'params.yaml'))
procDir = file.path(procParent, params$datasetName)
plotDir = file.path(resultDir, 'plots')
dir.create(plotDir, recursive = TRUE)
# load(file.path(resultDir, 'workspace.Rdata'))

registerDoParallel(cores = 2)

testing = TRUE
if (testing) {
  maxPvalLoad = 1
} else {
  maxPvalLoad = 1e-3}

############################################################

genoData = loadGeno(procDir, params$geno, file.path(resultDir, params$snpSubsetFile))
gwasMetadata = read_tsv(file.path(resultDir, 'gwas_metadata.tsv'), col_types = 'cccccdc')
setDT(gwasMetadata)

############################################################

gwasDataTmp = loadGwas(resultDir, gwasMetadata, maxPvalLoad)
gwasData = gwasDataTmp[[1]]
gwasLambdaData = gwasDataTmp[[2]]

gData = mergeAll(gwasData, phecodeData, gwasMetadata, genoData)

############################################################

gData[, plotManhattanAndQq(.BY, .SD, plotDir),
      by = .(phecode, phecodeStr, phenotype, method)]

############################################################

# maxPval = 1e-5
maxPval = 1
gDataSig = filterForSignificance(gData, maxPval)

############################################################

ptShp = 16
ptSz = 0.5
ptAlph = 0.2
lnSz = 0.5
lnCol = 'gray'

# hazard ratios are extremely similar to odds ratios
resultTmp = plotEffectSize(gDataSig, lnCol, lnSz, ptShp, ptSz, ptAlph)
gDataEffect = resultTmp[[1]]
pEffect = resultTmp[[2]]

cor(gDataEffect$logistic, gDataEffect$cox, use = 'na.or.complete')

# p-values at the low end are smaller
resultTmp = plotPval(gDataSig, lnCol, lnSz, ptShp, ptSz, ptAlph)
gDataPval = resultTmp[[1]]
pPval = resultTmp[[2]]

# standard errors are slightly lower for cox
resultTmp = plotSe(gDataSig)
gDataSe = resultTmp[[1]]
pSe = resultTmp[[2]]

# genomic inflation factors are similar or slightly higher
pLambda = plotLambda(gwasLambdaData, lnCol, lnSz, ptShp, ptSz, ptAlph)

# p = plot_grid(p1, p2, p3, p4, align = 'hv', axis = 'tb', nrow = 2)
# ggsave(file.path(resultDir, 'exome_full.pdf'), plot = p, width = 6, height = 5.75)

############################################################

a = merge(gDataPval, gDataSig, by = c('phecode', 'snp'))

p = ggplot(a) +
  geom_point(aes(x = log10(maf), y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph)
print(p)

p = ggplot(a) +
  geom_point(aes(x = log10(nCases), y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
  geom_smooth(aes(x = log10(nCases), y = cox - logistic), size = 0.5, method = 'loess', span = 0.5)
print(p)

p = ggplot(a) +
  geom_point(aes(x = logRatio, y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
  geom_smooth(aes(x = logRatio, y = cox - logistic), size = 0.5, method = 'loess', span = 0.5)
print(p)


a1 = dcast(gData, phecode + snp ~ method, value.var = 'pval')
a2 = a1[, .(r = cor(-log10(logistic), -log10(cox))), by = phecode]
p = ggplot(a2) +
  stat_ecdf(aes(x = r), pad = FALSE)
print(p)
