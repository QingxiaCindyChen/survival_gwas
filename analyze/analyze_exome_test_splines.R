source(file.path('analyze', 'analyze_setup.R'))

# testing = FALSE
# if (testing) {
#   resultDir = 'results_test'
# }

registerDoParallel(cores = 24)
# registerDoParallel(cores = 2)
filePrefix = 'exome'

# plinkFilepath = '../plink-1.90b6/plink'
# plinkMemSize = 100 # MB
# plinkVif = 1e6 # reduce plink's tendency to output NA due to multicollinearity of the spline bases

# genoDir = '../genotype_data/exome'
# genoPrefix = 'Exome_GRID_Euro'

minRecLen = 0 # years
minEvents = 2 # to be defined as a case; controls have zero events
maxAgeAtEvent = 90 # years; dates in the SD after this age are untrustworthy
minCases = 50 # to analyze the phecode

nPC = 2 # for both cox and logistic regression
buffer = 1 # years; for calculating age1 in cox regression

minMaf = 0.01
minCallRate = 0.95
minHwePval = 0.001

splineDf = 3:4 # for last_age in logistic regression (ns::splines)

############################################################
# load snp data

snpData = setDT(read_csv(file.path(procDir, 'exome_map.csv.gz'), col_types = cols()))
snpData = snpData[, .(snp = snp.name, chr = chromosome, pos = position)]

genoData = readRDS(file.path(procDir, sprintf('%s_genotype_data.rds', filePrefix)))

idx = (genoData$genoSummary$MAF >= minMaf) &
  (genoData$genoSummary$Call.rate >= minCallRate) &
  (2 * pnorm(-abs(genoData$genoSummary$z.HWE)) >= minHwePval)
snps = intersect(colnames(genoData$genoFull$genotypes)[idx], snpData[chr <= 22, snp])

snps = unique(read_csv('results/exome_pilot_top20.csv', col_types = cols())$snp)

# if (testing) {
  # snps = unique(read_csv('results/exome_pilot_top20.csv', col_types = cols())$snp)
  # set.seed(17)
  # snps = snps[sample.int(length(snps), 1000)]
# }

# snpFilepath = tempfile('snp_', fileext = '.tsv')
# write_tsv(data.table(snps), snpFilepath, col_names = FALSE)

############################################################
# load grid data

gridData = setDT(read_csv(file.path(procDir, sprintf('%s_grid_data.csv.gz', filePrefix)),
                          col_types = 'ccDDD'))
gridData[, first_age := time_length(first_entry_date - dob, 'years')]
gridData[, last_age := time_length(last_entry_date - dob, 'years')]
gridData[, rec_len := last_age - first_age]
gridData = gridData[first_age >= 0 & rec_len >= minRecLen]

# gridData = cbind(gridData, data.table(splines::ns(gridData$last_age, df = splineDf)))
# colnames(gridData)[(ncol(gridData) - splineDf + 1):ncol(gridData)] = paste0('last_age', 1:splineDf)

# covarColnames = c('rec_len', paste0('last_age', 1:splineDf))

# load PC data
if (nPC > 0) {
  pcData = setDT(read_csv(file.path(procDir, sprintf('%s_pc_data.csv.gz', filePrefix)), col_types = cols()))
  gridData = merge(gridData, pcData[, 1:(1 + nPC)], by = 'grid')
  # covarColnames = c(covarColnames, colnames(pcData)[2:(1 + nPC)])
}

genoTmp = data.table(grid = genoData$genoFull$fam$pedigree, sex = genoData$genoFull$fam$sex)
gridData = merge(gridData, genoTmp, by = 'grid')

# covarData = gridData[, .(FID = grid, IID = grid)]
# covarData = cbind(covarData, gridData[, c(covarColnames, 'sex'), with = FALSE])
#
# covarFilepath = tempfile('covar_', fileext = '.tsv')
# write_tsv(covarData, covarFilepath)

ageData = foreach(dfNow = splineDf) %dopar% {
  d = data.table(splines::ns(gridData$last_age, df = dfNow))
  colnames(d) = paste0('last_age', 1:dfNow)
  cbind(gridData[, .(grid)], d)
}

############################################################
# load phenotype data

phenoData = setDT(read_csv(file.path(procDir, sprintf('%s_phenotype_data.csv.gz', filePrefix)),
                           col_types = 'ccD'))

phenoData = merge(phenoData, gridData[, .(grid, dob, sex)], by = 'grid')
phenoData[, age := time_length(entry_date - dob, 'years')]
phenoData = phenoData[age <= maxAgeAtEvent]

phenoData = merge(phenoData, phecodeData[, .(phecode, whichSex)], by = 'phecode')
phenoData = phenoData[(whichSex == 'both') | (whichSex == 'male' & sex == 1) |
                        (whichSex == 'female' & sex == 2)]

phenoTmp = phenoData[, .N, by = .(grid, phecode)]
phenoTmp = phenoTmp[, .(nCases = sum(N >= minEvents)), by = phecode]
phenoData = merge(phenoData, phenoTmp[nCases >= minCases], by = 'phecode')
phenoData = phenoData[, .(grid, phecode, age)]

# if (testing) {
#   phenoData = phenoData[phecode %in% c('274.1', '714.1', '335', '185', '427.21', '290.11')]
# }

############################################################

phecodeDataKeep = phecodeData[phecode %in% unique(phenoData$phecode), .(phecode, whichSex)]
phecodeDataKeep[, phecodeStr := paste0('phe', gsub('.', 'p', phecode, fixed = TRUE))]

############################################################

d = foreach(ii = 1:nrow(phecodeDataKeep)) %dopar% {
  whichSex = phecodeDataKeep$whichSex[ii]
  phenoDataNow = phenoData[phecode == phecodeDataKeep$phecode[ii], .(grid, age)]
  inputBase1 = makeInput(phenoDataNow, gridData, whichSex, minEvents, buffer)

  d = foreach(jj = 1:length(splineDf)) %do% {
    glmStr = makeGlmStr(whichSex, nPC, splineDf[jj])
    inputBase2 = merge(inputBase1, ageData[[jj]], by = 'grid')

    d = foreach(snp = snps) %do% {
      inputNow = addSnpToInput(inputBase2, genoData$genoFull, snp)
      glmFit = glm(formula(glmStr), family = 'binomial', data = inputNow)
      d1 = coef(summary(glmFit))[2,]
      d2 = car::vif(glmFit)
      list(d1, d2)
    }

    d1 = data.table(do.call(rbind, lapply(d, function(a) a[[1]])))
    colnames(d1) = c('beta', 'se', 'z', 'pval')
    d1[, snp := snps]
    d1[, ageDf := splineDf[jj]]

    d2 = data.table(do.call(rbind, lapply(d, function(a) a[[2]])))
    d2[, snp := snps]
    d2[, ageDf := splineDf[jj]]
    d2 = melt(d2, id.vars = c('snp', 'ageDf'), value.name = 'vif', variable.factor = FALSE)

    list(d1, d2)
  }

  d1 = rbindlist(lapply(d, function(a) a[[1]]))
  d1[, phecode := phecodeDataKeep$phecode[ii]]
  d2 = rbindlist(lapply(d, function(a) a[[2]]))
  d2[, phecode := phecodeDataKeep$phecode[ii]]
  list(d1, d2)
}

d1 = rbindlist(lapply(d, function(a) a[[1]]))
d2 = rbindlist(lapply(d, function(a) a[[2]]))

save(d1, d2, file = file.path(resultDir, sprintf('%s_test_splines.Rdata', filePrefix)))
