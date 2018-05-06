source(file.path('analyze', 'analyze_setup.R'))

procDir = 'processed'
resultDir = 'results'
filePrefix = 'exome'

minEvents = 2
buffer = 1 # years, used for cox regression

# load phecode data
phecodeData = read_csv(file.path(procDir, 'phecode_definitions1.2.csv'), col_types = 'ccccii') %>%
	rename(phecode = jd_code,
			 phenotype = jd_string,
			 controlExcludeRange = jd_control_exclude_range,
			 whichSex = sex) %>%
	mutate(whichSex = tolower(ifelse(is.na(whichSex), 'both', whichSex)))

# load grid data
minRecLen = 0 # years
gridData = read_csv(file.path(procDir, 'grid_data.csv.gz'), col_types = 'ccDDD')
gridData = gridData %>%
	mutate(first_age = time_length(first_entry_date - dob, 'years'),
			 last_age = time_length(last_entry_date - dob, 'years'),
			 rec_len = last_age - first_age) %>%
	filter(first_age >= 0,
			 rec_len >= minRecLen)

# load snp data
minMaf = 0.01
minCallRate = 0.95
genoData = readRDS(file.path(procDir, 'genotype_data_exome.rds'))
idx = (genoData$genoSummary$MAF >= minMaf) & (genoData$genoSummary$Call.rate >= minCallRate)
snpData = read_csv(file.path(procDir, 'exome_map.csv.gz'), col_types = cols())
snps = intersect(colnames(genoData$genoFull$genotypes)[idx], snpData$snp.name[snpData$chromosome <= 22])

# load phenotype data
minGrids = 50
phenoData = read_csv(file.path(procDir, 'phenotype_data.csv.gz'), col_types = 'cDc') %>%
	group_by(phecode) %>%
	filter(n_distinct(grid) >= minGrids) %>%
	ungroup()

########################################
# testing

pilotPhenotypes = read_tsv(file.path(procDir, 'pilot_phenotypes.tsv'), col_types = 'cc')
phecodeData = semi_join(phecodeData, pilotPhenotypes, by = 'phecode')

# set.seed(4)
# snps = snps[sample.int(length(snps), 10)]

########################################

registerDoParallel(cores = 16)

logFilepath = file.path(resultDir, 'progress.txt')
timeStarted = Sys.time()
cat(sprintf('%s started analysis\n', timeStarted), file = logFilepath)


done = foreach(ii = 1:nrow(phecodeData)) %dopar% {
	whichSex = phecodeData$whichSex[ii]
	phenoDataNow = semi_join(phenoData, phecodeData[ii,], by = 'phecode')


	glmInput = makeGlmInput(phenoDataNow, gridData, genoData$genoFull, minEvents, whichSex)

	resultGlm = foreach(snp = snps, .combine = rbind) %dopar% {
		glmFit = runGlm(glmInput, genoData$genoFull, snp, whichSex)
		bind_cols(tibble(snp = snp), as_tibble(t(coef(summary(glmFit))[2,])))}

	resultGlm = resultGlm %>%
		transmute(method = 'logistic', snp = snp, coef = Estimate, se = `Std. Error`, z = `z value`,
					 pval = `Pr(>|z|)`)


	coxInput = makeCoxInput(phenoDataNow, gridData, genoData$genoFull, minEvents, buffer, whichSex)

	resultCox = foreach(snp = snps, .combine = rbind) %do% {
		coxFit = runCoxph(coxInput, genoData$genoFull, snp, whichSex)
		bind_cols(tibble(snp = snp), as_tibble(t(coef(summary(coxFit))[1,])))}

	resultCox = resultCox %>%
		transmute(method = 'cox', snp = snp, coef = coef, se = `se(coef)`, z = z, pval = `Pr(>|z|)`)


	resultNow = bind_rows(resultGlm, resultCox) %>%
		mutate(phecode = phecodeData$phecode[ii]) %>%
		arrange(method, pval)

	filename = sprintf('%s_minEvents%d_phe%s.csv.gz', filePrefix, minEvents,
							 gsub('.', 'p', phecodeData$phecode[ii], fixed = TRUE))
	write_csv(resultNow, gzfile(file.path(resultDir, filename)))

	cat(sprintf('%s completed phecode %d of %d\n', Sys.time(), ii, nrow(phecodeData)),
		 file = logFilepath, append = TRUE)
}

timeElapsed = Sys.time() - timeStarted
cat(sprintf('Time elapsed of %.2f %s\n', timeElapsed, attr(timeElapsed, 'units')),
	 file = logFilepath, append = TRUE)
