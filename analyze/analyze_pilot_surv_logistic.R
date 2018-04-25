library('readr')
library('dplyr')
library('lubridate')
library('doParallel')
library('snpStats')


# expected colnames in phenoRaw: grid, entry_date
# expected colnames in gridInfo: grid, last_age, rec_len
makeGlmInput = function(phenoRaw, gridInfo, exome, minEvents) {
	pheno = phenoRaw %>%
		distinct(grid, entry_date) %>%
		count(grid) %>%
		right_join(gridInfo, by = 'grid') %>%
		mutate(n = ifelse(is.na(n), 0, n)) %>%
		filter(n==0 | n >= minEvents) %>%
		mutate(status = as.integer(n >= minEvents))

	glmInput = tibble(grid = exome$fam$pedigree, sex = exome$fam$sex) %>%
		inner_join(pheno, by = 'grid') %>%
		select(grid, sex, last_age, rec_len, status)}


# expected colnames in coxInput: grid, sex, last_age, rec_len, status
runGlm = function(glmInput, exome, snp, splineDf = 4) {
	glmInput$snp = as(exome$genotypes[glmInput$grid, snp], 'numeric')[,1]
	glmFit = glm(status ~ snp + sex + rec_len + splines::ns(last_age, df = splineDf),
					 data = glmInput, family = binomial)}


# expected colnames in phenoRaw: grid, entry_date
# expected colnames in gridInfo: grid, dob, first_age, last_age
makeCoxInput = function(phenoRaw, gridInfo, exome, minEvents, buffer) {
	phenoCase = phenoRaw %>%
		semi_join(gridInfo, by = 'grid') %>%
		group_by(grid) %>%
		distinct(entry_date) %>%
		arrange(entry_date) %>%
		filter(row_number() == minEvents) %>%
		ungroup()

	phenoControl = gridInfo %>%
		select(grid) %>%
		anti_join(phenoRaw, by = 'grid')

	pheno = bind_rows(phenoCase, phenoControl) %>%
		right_join(gridInfo, by = 'grid') %>%
		mutate(age = time_length(entry_date - dob, 'years')) %>%
		select(grid, first_age, last_age, age)

	coxInput = tibble(grid = exome$fam$pedigree, sex = exome$fam$sex) %>%
		inner_join(pheno, by = 'grid') %>%
		group_by(grid) %>%
		mutate(status = as.integer(!is.na(age)),
				 age2 = ifelse(status, age, last_age),
				 age1 = min(first_age, age2 - buffer)) %>%
		ungroup() %>%
		select(grid, sex, age1, age2, status)}


# expected colnames in coxInput: grid, sex, age1, age2, status
runCoxph = function(coxInput, exome, snp) {
	coxInput$snp = as(exome$genotypes[coxInput$grid, snp], 'numeric')[,1]
	coxFit = coxph(Surv(age1, age2, status) ~ snp + sex, data = coxInput)}

########################################

registerDoParallel(cores = 16)

minMaf = 0.01
minCallRate = 0.95

phenotypes = c('alzheimers', 'atrial_fibrillation', 'gout',
					'multiple_sclerosis', 'prostate_cancer', 'rheumatoid_arthritis')
minEvents = 2
minRecLen = 0 # years
# minLastAge = 18 # years
# minAgeAtEvent = 18 # better?
buffer = 1 # years, used for cox regression

procDir = 'processed'
resultDir = 'results'

# load snp data
genotypeDir = '../genotype_data/exome'
exome = read.plink(file.path(genotypeDir, 'Exome_GRID_Euro'))
exomeSummary = col.summary(exome$genotypes)
idx = (exomeSummary$MAF >= minMaf) & (exomeSummary$Call.rate >= minCallRate)
snpInfo = read_csv(file.path(procDir, 'exome_map.csv.gz'), col_types = cols())
snps = intersect(colnames(exome$genotypes)[idx], snpInfo$snp.name[snpInfo$chromosome <= 22])

# set.seed(4)
# snps = snps[sample.int(length(snps), 10)]

# load grid data
gridInfo = read_csv(file.path(procDir, 'grid_info.csv.gz'), col_types = 'ccDTT')
colnames(gridInfo) = tolower(colnames(gridInfo))
gridInfo = gridInfo %>%
	rename(gender = gender_epic) %>%
	mutate(first_entry_date = as.Date(first_entry_date),
			 last_entry_date = as.Date(last_entry_date),
			 first_age = time_length(first_entry_date - dob, 'years'),
			 last_age = time_length(last_entry_date - dob, 'years'),
			 rec_len = last_age - first_age) %>%
	filter(first_age >= 0,
			 rec_len >= minRecLen)


for (phenotype in phenotypes) {
	phenoRaw = read_csv(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype)), col_types = 'ccT')
	colnames(phenoRaw) = tolower(colnames(phenoRaw))
	phenoRaw$entry_date = as.Date(phenoRaw$entry_date)


	coxInput = makeCoxInput(phenoRaw, gridInfo, exome, minEvents, buffer)

	coxDf = foreach(snp = snps, .combine = rbind) %dopar% {
		coxFit = runCoxph(coxInput, exome, snp)
		bind_cols(tibble(snp = snp), as_tibble(t(coef(summary(coxFit))[1,])))}

	coxDf = coxDf %>%
		rename(expCoef = `exp(coef)`, seCoef = `se(coef)`, pval = `Pr(>|z|)`) %>%
		arrange(pval) %>%
		write_csv(gzfile(file.path(resultDir, sprintf('surv%d_%s.csv.gz', minEvents, phenotype))))


	glmInput = makeGlmInput(phenoRaw, gridInfo, exome, minEvents)

	glmDf = foreach(snp = snps, .combine = rbind) %dopar% {
		glmFit = runGlm(glmInput, exome, snp)
		bind_cols(tibble(snp = snp), as_tibble(t(coef(summary(glmFit))[2,])))}

	glmDf = glmDf %>%
		rename(coef = Estimate, seCoef = `Std. Error`, z = `z value`, pval = `Pr(>|z|)`) %>%
		arrange(pval) %>%
		write_csv(gzfile(file.path(resultDir, sprintf('logistic%d_%s.csv.gz', minEvents, phenotype))))
}
