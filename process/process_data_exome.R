library('dplyr')
library('readr')
library('RODBC')
library('snpStats')

con = odbcConnect('NZSQL', believeNRows = FALSE)
gridSet = 'LB_EXOME'
procDir = 'processed'

phecodeIcdMapping = read_csv(file.path(procDir, 'phecode_icd9_map_unrolled.csv'), col_types = 'cc') %>%
	rename(icd = icd9)

phecodeData = read_csv(file.path(procDir, 'phecode_definitions1.2.csv'), col_types = 'ccccii') %>%
	rename(phecode = jd_code,
			 phenotype = jd_string,
			 controlExcludeRange = jd_control_exclude_range,
			 whichSex = sex) %>%
	mutate(whichSex = tolower(ifelse(is.na(whichSex), 'both', whichSex)))

########################################
# testing

pilotPhenotypes = read_tsv(file.path(procDir, 'pilot_phenotypes.tsv'), col_types = 'cc')
phecodeData = semi_join(phecodeData, pilotPhenotypes, by = 'phecode')

########################################

getPhenoRaw = function(con, icds, gridSet) {
	queryStr = sprintf("select a.GRID, a.CODE as icd, a.ENTRY_DATE
							 from ICD_CODES a
							 inner join %s b
							 on a.GRID = b.GRID
							 where a.CODE in ('%s')
							 and b.EURO = 1
							 order by a.GRID, a.ENTRY_DATE;",
							 gridSet, paste(icds, collapse="','"))
	queryStr = gsub(pattern = '\\s+', replacement = ' ', x = queryStr)
	d = as_tibble(sqlQuery(con, queryStr, stringsAsFactors = FALSE))
	colnames(d) = tolower(colnames(d))
	d = d %>%
		mutate(icd = as.character(icd),
				 entry_date = as.Date(entry_date))}

icds = semi_join(phecodeIcdMapping, phecodeData, by = 'phecode')$icd
phenoRaw = getPhenoRaw(con, icds, gridSet)
phenoMapped = phenoRaw %>%
	inner_join(phecodeIcdMapping, by = 'icd') %>%
	semi_join(phecodeData1, by = 'phecode') %>%
	distinct(grid, entry_date, phecode)
write_csv(phenoMapped, gzfile(file.path(procDir, 'phenotype_data.csv.gz')))

########################################

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

gridData = as_tibble(sqlQuery(con, queryStr, stringsAsFactors = FALSE))
colnames(gridData) = tolower(colnames(gridData))
gridInfo = gridData %>%
	mutate(first_entry_date = as.Date(first_entry_date),
			 last_entry_date = as.Date(last_entry_date))
write_csv(gridData, gzfile(file.path(procDir, 'grid_data.csv.gz')))

########################################

genotypeDir = '../genotype_data/exome'
genoFull = read.plink(file.path(genotypeDir, 'Exome_GRID_Euro'))
genoSummary = col.summary(genoFull$genotypes)
saveRDS(list(genoFull = genoFull, genoSummary = genoSummary), file.path(procDir, 'genotype_data_exome.rds'))
