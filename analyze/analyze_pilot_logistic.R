library('readr')
library('dplyr')
library('lubridate')
library('doParallel')
library('snpStats')

registerDoParallel(cores=16)

phenotypeNames = c('alzheimers', 'atrial_fibrillation', 'gout',
						 'multiple_sclerosis', 'prostate_cancer', 'rheumatoid_arthritis')
minMaf = 0.01
minEvents = 1
minRecLen = 1 # years

# load snp data
exome = read.plink(file.path('..', 'exome', 'Exome_GRID_Euro'))
exomeSummary = col.summary(exome$genotypes)
idx = exomeSummary$MAF >= minMaf
snpInfo = read_csv(file.path('processed', 'exome_map.csv'), col_types=cols())
snpNames = intersect(colnames(exome$genotypes)[idx], snpInfo$snp.name[snpInfo$chromosome <= 22])

# load grid data
gridInfo = read_csv(file.path('processed', 'grid_info.csv'), col_types='ccDTT')
colnames(gridInfo) = tolower(colnames(gridInfo))
gridInfo = gridInfo %>%
	rename(gender = gender_epic) %>%
	mutate(first_entry_date = as.Date(first_entry_date),
			 last_entry_date = as.Date(last_entry_date),
			 first_age = time_length(first_entry_date - dob, 'years'),
			 last_age = time_length(last_entry_date - dob, 'years'))

gridInfo = gridInfo %>%
	filter(last_age >= 18,
			 last_age - first_age >= minRecLen)

for (phenotypeName in phenotypeNames) {
	# load phenotype data
	phenoRaw = read_csv(file.path('processed', sprintf('pheno_%s.csv', phenotypeName)), col_types='ccT')
	colnames(phenoRaw) = tolower(colnames(phenoRaw))
	phenoRaw$entry_date = as.Date(phenoRaw$entry_date)

	# prepare phenotype data for analysis
	pheno = phenoRaw %>%
		count(grid) %>%
		right_join(gridInfo, by='grid') %>%
		mutate(n = ifelse(is.na(n), 0, n)) %>%
		filter(n==0 | n >= minEvents) %>%
		mutate(status = as.integer(n >= minEvents))

	# run glm analysis
	glmBase = tibble(grid = exome$fam$pedigree, sex = exome$fam$sex) %>%
		inner_join(pheno, by='grid') %>%
		select(grid, sex, last_age, status)

	glmDf = foreach(snpName=snpNames, .combine = rbind) %dopar% {
		glmInput = glmBase
		glmInput$snp = as(exome$genotypes[glmBase$grid, snpName], 'numeric')[,1]
		glmFit = glm(status ~ snp + sex + splines::ns(last_age, df=4), data=glmInput, family=binomial)
		tibble(snpName = snpName) %>%
			bind_cols(as_tibble(t(coef(summary(glmFit))[2,])))}

	glmDf = glmDf %>%
		rename(coef = Estimate, seCoef = `Std. Error`, z = `z value`, pval = `Pr(>|z|)`) %>%
		arrange(pval) %>%
		write_csv(file.path('results', sprintf('logistic%d_%s.csv', minEvents, phenotypeName)))
}
