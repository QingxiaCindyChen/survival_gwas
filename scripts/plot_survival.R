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
maxPvalLoad = 1e-5

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
# load snp data

tmpData = loadSnpGenoData(params$geno, params$plink$dataPathPrefix,
                          file.path(resultDir, params$snpSubsetFile))
snpData = tmpData$snpData
genoData = tmpData$genoData
rm(tmpData)

############################################################
# load grid data

gridTmp = loadGrid(procDir, params$plink$dataPathPrefix,
                   params$pheno$minRecLen, params$gwas)
gridData = gridTmp[[1]]
covarColnames = gridTmp[[2]]
rm(gridTmp)

############################################################
# load phenotype data

phenoTmp = loadPheno(procDir, params$pheno, gridData,
                     file.path(resultDir, params$phecodeSubsetFile))
phenoData = phenoTmp[[1]]
phenoSummary = phenoTmp[[2]]
rm(phenoTmp)

############################################################

survPhecodeSnp = data.table(phecode = c('185', '290.11', '335'),
                            snp = c('rs7931342', 'rs157582', 'rs3129889'))

sfdPre = foreach(ii = 1:nrow(survPhecodeSnp), .combine = rbind) %dopar% {
  phecodeNow = survPhecodeSnp$phecode[ii]
  snpNow = survPhecodeSnp$snp[ii]
  sfNow = getSurvfit(phecodeNow, snpNow, gwasMetadata, phenoData, gridData,
                     whichSex, params$pheno, genoData)
  sfdNow = getSurvfitDt(sfNow)
  sfdNow[, phecode := phecodeNow]
  sfdNow[, snp := snpNow]
  sfdNow}

############################################################

sfd = merge(sfdPre, phecodeData[, .(phecode, phenotype)], by = 'phecode')
sfd = merge(sfd, mapData[, .(snp, allele1, allele2)], by = 'snp')

sfd[, phenoLabel := sprintf('%s (%s)', phenotype, phecode)]
sfd[, phenoLabel := factor(phenoLabel, sort(unique(sfd$phenoLabel),
                                            decreasing = TRUE))]
sfd[, snpLabel := sprintf('%s-%s', snp, allele1)]

p = ggplot(sfd) +
  facet_wrap(~ phenoLabel + snpLabel, scales = 'free_y', nrow = 1) +
  geom_step(aes(x = age, y = surv, color = factor(genotype)), size = 1) +
  labs(x = 'Age (y)', y = 'Survival\n(frac. undiagnosed)',
       color = 'Allele\ncount') +
  scale_x_continuous(limits = c(20, 89)) +
  scale_color_viridis(direction = 1, discrete = TRUE)

ggsave(file.path(plotDir, 'example_survival.pdf'),
       plot = p, width = 7, height = 2.5)
