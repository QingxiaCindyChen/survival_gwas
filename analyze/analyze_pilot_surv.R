library('readr')
library('dplyr')
library('lubridate')
library('doParallel')
library('snpStats')

registerDoParallel(cores = 16)

minMaf = 0.01
minCallRate = 0.95

phenotypes = c('alzheimers', 'atrial_fibrillation', 'gout',
					'multiple_sclerosis', 'prostate_cancer', 'rheumatoid_arthritis')
minEvents = 2
buffer = 1 # years

procDir = 'processed'
resultDir = 'results'

# load snp data
genotypeDir = '../genotype_data/exome'
exome = read.plink(file.path(genotypeDir, 'Exome_GRID_Euro'))
exomeSummary = col.summary(exome$genotypes)
idx = (exomeSummary$MAF >= minMaf) & (exomeSummary$Call.rate >= minCallRate)
snpInfo = read_csv(file.path(procDir, 'exome_map.csv.gz'), col_types=cols())
snpNames = intersect(colnames(exome$genotypes)[idx], snpInfo$snp.name[snpInfo$chromosome <= 22])

# load grid data
gridInfo = read_csv(file.path(procDir, 'grid_info.csv.gz'), col_types='ccDTT')
colnames(gridInfo) = tolower(colnames(gridInfo))
gridInfo = gridInfo %>%
	rename(gender = gender_epic) %>%
	mutate(first_entry_date = as.Date(first_entry_date),
			 last_entry_date = as.Date(last_entry_date),
			 first_age = time_length(first_entry_date - dob, 'years'),
			 last_age = time_length(last_entry_date - dob, 'years'))

# functions:
# make coxDf for a phenotype: phenoRaw, exome, gridInfo, snpNames, buffer, minEvents
# make coxBase for a phenotype: phenoRaw, gridInfo, exome, buffer, minEvents
# make coxFit for a phenotype and snp: coxBase, exome, snpName

for (phenotype in phenotypes) {
	# load phenotype data
	phenoRaw = read_csv(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype)), col_types='ccT')
	colnames(phenoRaw) = tolower(colnames(phenoRaw))
	phenoRaw$entry_date = as.Date(phenoRaw$entry_date)

	# prepare phenotype data for analysis
	phenoCase = phenoRaw %>%
		group_by(grid) %>%
		distinct(entry_date) %>%
		arrange(entry_date) %>%
		filter(row_number() == minEvents) %>%
		ungroup()

	phenoControl = gridInfo %>%
		select(grid) %>%
		anti_join(phenoRaw, by='grid')

	pheno = bind_rows(phenoCase, phenoControl) %>%
		right_join(gridInfo, by='grid') %>%
		mutate(age = time_length(entry_date - dob, 'years')) %>%
		select(grid, gender, first_age, last_age, age)

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
		coxFit = coxph(Surv(age1, age2, status) ~ snp + sex, data=coxInput)
		tibble(snpName = snpName) %>%
			bind_cols(as_tibble(t(coef(summary(coxFit))[1,])))}

	coxDf = coxDf %>%
		rename(expCoef = `exp(coef)`, seCoef = `se(coef)`, pval = `Pr(>|z|)`) %>%
		arrange(pval) %>%
		write_csv(gzfile(file.path(resultDir, sprintf('surv%d_%s.csv.gz', minEvents, phenotype))))
}
