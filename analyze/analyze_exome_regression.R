source(file.path('analyze', 'analyze_setup.R'))

filePrefix = 'exome'

minEvents = 2
maxAgeAtEvent = 90 # years; dates in the SD after this age are untrustworthy
buffer = 1 # years; used for cox regression

# load grid data
minRecLen = 0 # years
gridData = setDT(read_csv(file.path(procDir, sprintf('%s_grid_data.csv.gz', filePrefix)),
								  col_types = 'ccDDD'))
gridData[, first_age := time_length(first_entry_date - dob, 'years')]
gridData[, last_age := time_length(last_entry_date - dob, 'years')]
gridData[, rec_len := last_age - first_age]
gridData = gridData[first_age >= 0 & rec_len >= minRecLen,]

# load snp data
minMaf = 0.01
minCallRate = 0.95
minHwePval = 0.001
snpData = setDT(read_csv(file.path(procDir, 'exome_map.csv.gz'), col_types = cols()))
snpData = snpData[,.(snp = snp.name, chr = chromosome, pos = position)]

genoData = readRDS(file.path(procDir, sprintf('%s_genotype_data.rds', filePrefix)))
idx = (genoData$genoSummary$MAF >= minMaf) &
	(genoData$genoSummary$Call.rate >= minCallRate) &
	(2 * pnorm(-abs(genoData$genoSummary$z.HWE)) >= minHwePval)
snps = intersect(colnames(genoData$genoFull$genotypes)[idx], snpData[chr <= 22, snp])

# load phenotype data
minGrids = 50
phenoData = setDT(read_csv(file.path(procDir, sprintf('%s_phenotype_data.csv.gz', filePrefix)),
									col_types = 'ccD'))
phenoData = phenoData[, if (length(unique(grid)) >= minGrids) .SD, by = phecode]
phenoData = phenoData[order(grid, phecode, entry_date), .(grid, phecode, entry_date)]

############################################################
# testing

# phenoData = phenoData[phecode %in% c('274.1', '714.1', '335', '185', '427.21', '290.11'),]
# snps = unique(read_csv(file.path(resultDir, 'exome_pilot_top20.csv'))$snp)

############################################################

phenoData = merge(phenoData, gridData, by = 'grid')
phenoData[, age := time_length(entry_date - dob, 'years')]
phenoData = phenoData[age <= maxAgeAtEvent,]

############################################################

registerDoParallel(cores = 20)
# registerDoParallel(cores = 2)

phenoLoop = phecodeData[phecode %in% unique(phenoData$phecode),]

logFilepath = file.path(resultDir, sprintf('%s_progress.txt', filePrefix))
timeStarted = Sys.time()
cat(sprintf('%s started analysis\n', timeStarted), file = logFilepath)

gwasFilenames = foreach(ii = 1:nrow(phenoLoop), .combine = c) %dopar% {
	phecode = phenoLoop$phecode[ii]
	whichSex = phenoLoop$whichSex[ii]
	phenoDataNow = phenoData[phecode == phenoLoop$phecode[ii], .(grid, age)]

	if (whichSex == 'both') {
		runGlmNow = runGlmSexBoth
		runCoxphNow = runCoxphSexBoth
	} else {
		runGlmNow = runGlmSexOne
		runCoxphNow = runCoxphSexOne}

	inputBase = makeInput(phenoDataNow, gridData, genoData$genoFull, whichSex, minEvents, buffer)

	resultNow = foreach(snp = snps, .combine = rbind) %do% {
		inputNow = addSnpToInput(inputBase, genoData$genoFull, snp)
		glmFit = runGlmNow(inputNow)
		coxFit = runCoxphNow(inputNow)
		rbind(coef(summary(glmFit))[2, 1:3], coef(summary(coxFit))[1, c(1, 3:4)])}

	setDT(resultNow)
	colnames(resultNow) = c('coef', 'se', 'z')

	# speedglm encodes pvals as factors, and coxph sets small pvals to zero.
	resultNow[, pval := 2 * pnorm(-abs(z))]
	resultNow[, method := rep.int(c('logistic', 'cox'), times = length(snps))]
	resultNow[, snp := rep(snps, each = 2)]
	resultNow[, phecode := phenoLoop$phecode[ii]]

	filename = sprintf('%s_phe%s.csv.gz', filePrefix, gsub('.', 'p', phecode, fixed = TRUE))
	write_csv(resultNow, gzfile(file.path(resultDir, filename)))

	cat(sprintf('%s completed phecode %s (%d of %d)\n', Sys.time(), phecode, ii, nrow(phenoLoop)),
		 file = logFilepath, append = TRUE)
	filename}

timeElapsed = Sys.time() - timeStarted
cat(sprintf('Time elapsed of %.2f %s\n', timeElapsed, attr(timeElapsed, 'units')),
	 file = logFilepath, append = TRUE)

save(list = setdiff(ls(), 'genoData'),
	  file = file.path(resultDir, sprintf('%s_workspace.Rdata', filePrefix)))
