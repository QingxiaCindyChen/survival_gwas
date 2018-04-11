library('readr')
library('dplyr')
library('lubridate')
library('doParallel')
library('snpStats')

registerDoParallel(cores=16)

phenotypeNames = c('alzheimers', 'atrial_fibrillation', 'gout',
						 'multiple_sclerosis', 'prostate_cancer', 'rheumatoid_arthritis')
buffer = 1 # years
minMaf = 0.01
nEvents = 1

# load snp data
genotypeDir = '../genotype_data/exome'
exome = read.plink(file.path(genotypeDir, 'Exome_GRID_Euro'))
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

for (phenotypeName in phenotypeNames) {
	# load phenotype data
	phenoRaw = read_csv(file.path('processed', sprintf('pheno_%s.csv', phenotypeName)), col_types='ccT')
	colnames(phenoRaw) = tolower(colnames(phenoRaw))
	phenoRaw$entry_date = as.Date(phenoRaw$entry_date)

	# prepare phenotype data for analysis
	pheno = phenoRaw %>%
		right_join(gridInfo, by='grid') %>%
		mutate(age = time_length(entry_date - dob, 'years')) %>%
		group_by(grid, gender, first_age, last_age) %>%
		summarize(age = min(age))

	# run survival analysis
	coxBase = tibble(grid = exome$fam$pedigree, sex = exome$fam$sex) %>%
		inner_join(pheno, by='grid') %>%
		group_by(grid) %>%
		mutate(status = as.integer(!is.na(age)),
				 age2 = ifelse(status, age, last_age),
				 age1 = min(first_age, age2 - buffer)) %>%
		ungroup() %>%
		select(grid, sex, age1, age2, status)

	coxDf = foreach(snpName=snpNames, .combine = rbind) %dopar% {
		coxInput = coxBase
		coxInput$snp = as(exome$genotypes[coxBase$grid, snpName], 'numeric')[,1]
		coxFit = coxph(Surv(age1, age2, status) ~ snp + sex + cluster(grid), data=coxInput)
		tibble(snpName = snpName) %>%
			bind_cols(as_tibble(t(coef(summary(coxFit))[1,])))}

	coxDf = coxDf %>%
		rename(expCoef = `exp(coef)`, seCoef = `se(coef)`, pval = `Pr(>|z|)`) %>%
		arrange(pval) %>%
		write_csv(file.path('results', sprintf('surv%d_%s.csv', nEvents, phenotypeName)))
}
