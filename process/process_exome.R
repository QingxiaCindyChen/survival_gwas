library(data.table)
library(readr)
library(RODBC)
library(snpStats)
library(SNPRelate)

gridSet = 'LB_EXOME'
procDir = 'processed'
filePrefix = 'exome'

genotypeDir = '../genotype_data/exome'
genoPrefix = 'Exome_GRID_Euro'

phecodeIcdMapping = fread(file.path(procDir, 'phecode_icd9_map_unrolled.csv'))
setnames(phecodeIcdMapping, old = 'icd9', new = 'icd')

phecodeData = fread(file.path(procDir, 'phecode_definitions1.2.csv'))
phecodeData = phecodeData[,.(phecode = jd_code,
									  phenotype = jd_string,
									  controlExcludeRange = jd_control_exclude_range,
									  whichSex = tolower(ifelse(sex == '', 'both', sex)),
									  rollup, leaf)]

con = odbcConnect('NZSQL', believeNRows = FALSE)

############################################################

getPhenoRaw = function(con, icds, gridSet) {
	queryStr = sprintf("select distinct a.GRID, a.CODE as icd, a.ENTRY_DATE
							 from ICD_CODES a
							 inner join %s b
							 on a.GRID = b.GRID
							 where a.CODE in ('%s')
							 and b.EURO = 1
							 order by a.GRID, a.ENTRY_DATE;",
							 gridSet, paste(icds, collapse="','"))
	queryStr = gsub(pattern = '\\s+', replacement = ' ', x = queryStr)
	d = setDT(sqlQuery(con, queryStr, stringsAsFactors = FALSE, as.is = TRUE))
	colnames(d) = tolower(colnames(d))
	d[, entry_date := as.Date(entry_date)]
	d}

icds = unique(phecodeIcdMapping[phecode %in% phecodeData$phecode, icd]) # in case phecodeData is a subset
phenoRaw = getPhenoRaw(con, icds, gridSet)

phenoMapped = merge(phenoRaw, phecodeIcdMapping, by = 'icd', allow.cartesian = TRUE)
phenoMapped = merge(phenoMapped, phecodeData[, .(phecode)], by = 'phecode') # in case phecodeData is a subset
phenoMapped = unique(phenoMapped[, .(grid, phecode, entry_date)])

write_csv(phenoMapped, gzfile(file.path(procDir, sprintf('%s_phenotype_data.csv.gz', filePrefix))))

############################################################

queryStr = sprintf('select a.GRID, c.GENDER_EPIC as gender, c.DOB,
						 min(a.ENTRY_DATE) as FIRST_ENTRY_DATE, max(a.ENTRY_DATE) as LAST_ENTRY_DATE
						 from ICD_CODES a
						 inner join %s b
						 on a.GRID = b.GRID
						 inner join SD_RECORD c
						 on a.GRID = c.GRID
						 where b.EURO = 1
						 group by a.GRID, c.GENDER_EPIC, c.DOB;',
						 gridSet)
queryStr = gsub(pattern = '\\s+', replacement = ' ', x = queryStr)

gridData = setDT(sqlQuery(con, queryStr, stringsAsFactors = FALSE, as.is = TRUE))
colnames(gridData) = tolower(colnames(gridData))

gridData[, first_entry_date := as.Date(first_entry_date)]
gridData[, last_entry_date := as.Date(last_entry_date)]

write_csv(gridData, gzfile(file.path(procDir, sprintf('%s_grid_data.csv.gz', filePrefix))))

############################################################

genoFull = read.plink(file.path(genotypeDir, genoPrefix))
genoSummary = col.summary(genoFull$genotypes)

saveRDS(list(genoFull = genoFull, genoSummary = genoSummary),
		  file.path(procDir, sprintf('%s_genotype_data.rds', filePrefix)))

snpgdsBED2GDS(file.path(genotypeDir, paste0(genoPrefix, '.bed')),
				  file.path(genotypeDir, paste0(genoPrefix, '.fam')),
				  file.path(genotypeDir, paste0(genoPrefix, '.bim')),
				  file.path(genotypeDir, paste0(genoPrefix, '.gds')))
genoFile = snpgdsOpen(file.path(genotypeDir, paste0(genoPrefix, '.gds')))

minMaf = 0.01
minCallRate = 0.95
minHwePval = 0.001
idx = (genoSummary$MAF >= minMaf) &
	(genoSummary$Call.rate >= minCallRate) &
	(2 * pnorm(-abs(genoSummary$z.HWE)) >= minHwePval)

aims = read_csv(file.path(procDir, 'exome_aims.txt'), col_names = FALSE, col_types = cols())$X1
snps = intersect(colnames(genoFull$genotypes)[idx], aims)

pcaResult = snpgdsPCA(genoFile, snp.id = snps, eigen.cnt = 10, num.thread = 8, algorithm = 'randomized')

pcData = data.table(pcaResult$sample.id, pcaResult$eigenvect)
colnames(pcData) = c('grid', paste0('PC', 1:ncol(pcaResult$eigenvect)))

write_csv(pcData, gzfile(file.path(procDir, sprintf('%s_pc_data.csv.gz', filePrefix))))
