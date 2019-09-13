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
maxPvalLoad = 1e-4

############################################################

gwasMetadata = read_tsv(file.path(resultDir, 'gwas_metadata.tsv'),
                        col_types = 'cccccdddc')
setDT(gwasMetadata)

mapData = read_csv(file.path(procDir, 'map_data.csv.gz'), col_types = 'iccicc')
setDT(mapData)

############################################################

gwasDataTmp = loadGwas(resultDir, gwasMetadata, maxPvalLoad)
gwasData = gwasDataTmp[[1]]

gwasData[, pval := ifelse(is.na(pval), 1, pval)]
gwasData = mergeAll(gwasData, phecodeData, gwasMetadata, mapData)

rm(gwasDataTmp)

############################################################

catalogAssocData = read_tsv(file.path(procDir, 'catalog_assoc_data.tsv'),
                            col_types = 'cccc')
setDT(catalogAssocData)

nAssoc = nrow(unique(catalogAssocData[, .(phecode, ldBlock)]))

a = unique(merge(gwasData[, .(phecode, snp, method, pval)],
                 catalogAssocData[, .(phecode, snp, ldBlock)],
                 by = c('phecode', 'snp')))

############################################################

pvalSeqCutoff = 1e-4 # must be <= maxPvalLoad
aCast = dcast(a, phecode + ldBlock + snp ~ method, value.var = 'pval')
aCast[, sequential := ifelse(logistic <= pvalSeqCutoff, cox, logistic)]
a = melt(aCast, measure.vars = c('logistic', 'cox', 'sequential'),
         variable.name = 'method', value.name = 'pval')

############################################################

negLogPvalCutoffs = seq(5, 30, 0.02)

a1 = foreach(negLogPvalNow = negLogPvalCutoffs, .combine = rbind) %dopar% {
  aNow = a[pval <= 10^(-negLogPvalNow),
           .(nDetected = uniqueN(.SD, by = c('phecode', 'ldBlock'))),
           by = method]
  aNow[, negLogPvalCutoff := negLogPvalNow][]}

a1$method = factor(a1$method, levels(a1$method), c('logistic', 'Cox', 'sequential'))
a1[, fracDetected := nDetected / nAssoc]

p1 = ggplot(a1) +
  geom_line(aes(x = negLogPvalCutoff, y = fracDetected, color = method,
                linetype = method), size = 0.8) +
  labs(x = expression(-log[10](p)~cutoff),
       y = 'Sensitivity for\nknown associations', color = 'Method',
       linetype = 'Method') +
  scale_x_continuous(limits = c(5, 9)) +
  scale_y_continuous(limits = c(0.02, 0.055)) +
  scale_color_manual(values = c('#33a02c', '#1f78b4', '#a6cee3')) +
  scale_linetype_manual(values = c('solid', 'solid', 'dashed')) +
  guides(color = guide_legend(keywidth = 2, keyheight = 1),
         linetype = guide_legend(keywidth = 2, keyheight = 1))


a2 = a1[, .(relDiffDetected = (nDetected[method == 'Cox'] -
                                 nDetected[method == 'logistic']) /
              nDetected[method == 'logistic']),
        by = negLogPvalCutoff]

p2 = ggplot(a2) +
  geom_line(aes(x = negLogPvalCutoff, y = relDiffDetected, color = 'raw'),
            size = 0.8) +
  geom_smooth(aes(x = negLogPvalCutoff, y = relDiffDetected, color = 'smoothed'),
              span = 0.5, se = FALSE, size = 0.75,
              show.legend = TRUE) +
  labs(x = expression(-log[10](p)~cutoff),
       y = 'Rel. change in sensitivity\nfrom logistic to Cox') +
  scale_x_continuous(limits = c(5, 9)) +
  scale_y_continuous(limits = c(0, 0.2)) +
  scale_color_manual(name = 'Type', values = c('darkgray', 'black')) +
  guides(color = guide_legend(keywidth = 2, keyheight = 1))

p = plot_grid(p1, p2, ncol = 1, labels = 'AUTO', align = 'v', axis = 'lr')
ggsave(file.path(plotDir, 'summary_catalog_sensitivity.pdf'),
       plot = p, width = 4.5, height = 5.1)
