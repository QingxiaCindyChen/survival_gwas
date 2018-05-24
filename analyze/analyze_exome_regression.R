source(file.path('analyze', 'analyze_setup.R'))

filePrefix = 'exome'

minRecLen = 0 # years
minEvents = 2 # to be defined as a case; controls have zero events
maxAgeAtEvent = 90 # years; dates in the SD after this age are untrustworthy
minCases = 50 # to analyze the phecode

nPC = 2 # for both cox and logistic regression
splineDf = 4 # for last_age in logistic regression (ns::splines)
buffer = 1 # years; for calculating age1 in cox regression

minMaf = 0.01
minCallRate = 0.95
minHwePval = 0.001

############################################################
# load snp data

snpData = setDT(read_csv(file.path(procDir, 'exome_map.csv.gz'), col_types = cols()))
snpData = snpData[,.(snp = snp.name, chr = chromosome, pos = position)]

genoData = readRDS(file.path(procDir, sprintf('%s_genotype_data.rds', filePrefix)))

idx = (genoData$genoSummary$MAF >= minMaf) &
	(genoData$genoSummary$Call.rate >= minCallRate) &
	(2 * pnorm(-abs(genoData$genoSummary$z.HWE)) >= minHwePval)
snps = intersect(colnames(genoData$genoFull$genotypes)[idx], snpData[chr <= 22, snp])

# testing
# snps = unique(read_csv(file.path(resultDir, 'exome_pilot_top20.csv'), col_types = cols())$snp)
# set.seed(17)
# snps = snps[sample.int(length(snps), 5000)]

############################################################
# load grid data

gridData = setDT(read_csv(file.path(procDir, sprintf('%s_grid_data.csv.gz', filePrefix)),
								  col_types = 'ccDDD'))
gridData[, first_age := time_length(first_entry_date - dob, 'years')]
gridData[, last_age := time_length(last_entry_date - dob, 'years')]
gridData[, rec_len := last_age - first_age]
gridData = gridData[first_age >= 0 & rec_len >= minRecLen,]

gridData = cbind(gridData, data.table(splines::ns(gridData$last_age, df = splineDf)))
colnames(gridData)[(ncol(gridData) - splineDf + 1):ncol(gridData)] = paste0('last_age', 1:splineDf)

# load PC data
pcData = setDT(read_csv(file.path(procDir, sprintf('%s_pc_data.csv.gz', filePrefix)), col_types = cols()))
if (nPC > 0) {
	gridData = merge(gridData, pcData[, 1:(1 + nPC)], by = 'grid')}

genoTmp = data.table(grid = genoData$genoFull$fam$pedigree, sex = genoData$genoFull$fam$sex)
gridData = merge(gridData, genoTmp, by = 'grid')

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

# testing
# phenoData = phenoData[phecode %in% c('274.1', '714.1', '335', '185', '427.21', '290.11')]

############################################################

phecodeDataKeep = phecodeData[phecode %in% unique(phenoData$phecode), .(phecode, whichSex)]
phecodeDataKeep[, filename := sprintf('%s_phe%s.csv.gz', filePrefix, gsub('.', 'p', phecode, fixed = TRUE))]

# genoFlat = data.table(grid = gridData$grid, as(genoData$genoFull$genotypes[gridData$grid, snps], 'numeric'))
# genoFlat = melt(genoFlat, id.vars = 'grid', variable.name = 'snpId', value.name = 'snp')
# genoFlat = genoFlat[!is.na(snp)]
# genoFlat[, snp := as.integer(snp)]

methodVec = rep.int(c('logistic', 'cox'), times = length(snps))
snpVec = rep(snps, each = 2)

############################################################

registerDoParallel(cores = 20)

logFilepath = file.path(resultDir, sprintf('%s_progress.txt', filePrefix))
timeStarted = Sys.time()
cat(sprintf('%s started analysis\n', timeStarted), file = logFilepath)

# microbenchmark::microbenchmark({

done = foreach(ii = 1:nrow(phecodeDataKeep)) %dopar% {
	whichSex = phecodeDataKeep$whichSex[ii]
	phenoDataNow = phenoData[phecode == phecodeDataKeep$phecode[ii], .(grid, age)]

	inputBase = makeInput(phenoDataNow, gridData, whichSex, minEvents, buffer)
	glmStr = getGlmStr(whichSex, nPC, splineDf)
	coxStr = getCoxStr(whichSex, nPC)

	resultNow = foreach(snp = snps, .combine = rbind) %do% {
		inputNow = addSnpToInput(inputBase, genoData$genoFull, snp)
		# inputNow = addSnpToInput1(inputBase, genoFlat, snp)
		glmFit = runGlm(glmStr, inputNow)
		coxFit = runCox(coxStr, inputNow)
		data.table(rbind(coef(summary(glmFit))[2, 1:3], coef(summary(coxFit))[1, c(1, 3:4)]))}

	colnames(resultNow) = c('coef', 'se', 'z')

	# speedglm encodes pvals as factors, and coxph sets small pvals to zero.
	resultNow[, pval := 2 * pnorm(-abs(z))]
	resultNow[, method := methodVec]
	resultNow[, snp := snpVec]
	resultNow[, phecode := phecodeDataKeep$phecode[ii]]

	write_csv(resultNow, gzfile(file.path(resultDir, phecodeDataKeep$filename[ii])))
	cat(sprintf('%s completed phecode %s (%d of %d)\n', Sys.time(), phecodeDataKeep$phecode[ii],
					ii, nrow(phecodeDataKeep)), file = logFilepath, append = TRUE)
}

# }, times = 1)

timeElapsed = Sys.time() - timeStarted
cat(sprintf('Time elapsed of %.2f %s\n', timeElapsed, attr(timeElapsed, 'units')),
	 file = logFilepath, append = TRUE)

save(list = setdiff(ls(), 'genoData'),
	  file = file.path(resultDir, sprintf('%s_workspace.Rdata', filePrefix)))
