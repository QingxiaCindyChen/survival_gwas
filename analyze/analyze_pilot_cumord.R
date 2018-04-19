library('readr')
library('dplyr')
library('lubridate')
library('doParallel')
library('ordinal')
library('snpStats')

registerDoParallel(cores = 16)

minMaf = 0.01
minCallRate = 0.95

phenotypes = c('alzheimers', 'atrial_fibrillation', 'gout',
					'multiple_sclerosis', 'prostate_cancer', 'rheumatoid_arthritis')
maxEvents = 3
minRecLen = 1 # years

procDir = 'processed'
resultDir = 'results'

# load snp data
exome = read.plink(file.path('..', 'exome', 'Exome_GRID_Euro'))
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

gridInfo = gridInfo %>%
	filter(last_age >= 18,
			 last_age - first_age >= minRecLen)

for (phenotype in phenotypes) {
	# load phenotype data
	phenoRaw = read_csv(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype)), col_types='ccT')
	colnames(phenoRaw) = tolower(colnames(phenoRaw))
	phenoRaw$entry_date = as.Date(phenoRaw$entry_date)

	# prepare phenotype data for analysis
	pheno = phenoRaw %>%
		count(grid) %>%
		right_join(gridInfo, by='grid') %>%
		mutate(n = ifelse(is.na(n), 0, n),
				 status = factor(ifelse(n >= maxEvents, maxEvents, n), ordered=TRUE))

	# run clm analysis
	clmBase = tibble(grid = exome$fam$pedigree, sex = exome$fam$sex) %>%
		inner_join(pheno, by='grid') %>%
		select(grid, sex, last_age, status)

	clmDf = foreach(snpName=snpNames, .combine = rbind) %dopar% {
		clmInput = clmBase
		clmInput$snp = as(exome$genotypes[clmBase$grid, snpName], 'numeric')[,1]
		clmFit = clm(status ~ snp + sex + splines::ns(last_age, df=4), link='logit', data=clmInput)
		tibble(snpName = snpName) %>%
			bind_cols(as_tibble(t(coef(summary(clmFit))[maxEvents + 1,])))}

	clmDf = clmDf %>%
		rename(coef = Estimate, seCoef = `Std. Error`, z = `z value`, pval = `Pr(>|z|)`) %>%
		arrange(pval) %>%
		write_csv(gzfile(file.path(resultDir, sprintf('cumord%d_%s.csv.gz', maxEvents, phenotype))))
}
