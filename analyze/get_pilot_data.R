# library('tidyverse')
library('readr')
library('RODBC')

# need to clean this up, include mapping and maybe save everything as one big dataframe

# mapUnrolled = read_csv(file.path(procDir, 'phecode_icd9_map_unrolled.csv'), col_types='cc')
# phecodeDef = read_csv(file.path(procDir, 'phecode_definitions1.2.csv'), col_types='ccccii')

con = odbcConnect('NZSQL', believeNRows = FALSE)
gridSet = 'LB_EXOME'
procDir = 'processed'

getPhenoRaw = function(con, phecodes, gridSet) {
	queryStr = sprintf("select a.GRID, a.CODE, a.ENTRY_DATE
							 from ICD_CODES a
							 inner join %s b
							 on a.GRID = b.GRID
							 where a.CODE in ('%s')
							 and b.EURO = 1
							 order by a.GRID, a.ENTRY_DATE;",
							 gridSet, paste(phecodes, collapse="','"))
	queryStr = gsub(pattern = '\\s+', replacement = ' ', x = queryStr)
	df = sqlQuery(con, queryStr)}

########################################

queryStr = sprintf('select a.GRID, c.GENDER_EPIC, c.DOB,
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

gridInfo = sqlQuery(con, queryStr)
write_csv(gridInfo, gzfile(file.path(procDir, 'grid_info.csv.gz')))

########################################

phenotype = 'gout'
phecodes = c('274', '274.0', '274.00', '274.01', '274.02', '274.03', '274.1', '274.10',
				 '274.11', '274.19', '274.8', '274.81', '274.82', '274.89', '274.9')
phenoRaw = getPhenoRaw(con, phecodes, gridSet)
write_csv(phenoRaw, gzfile(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype))))

########################################

phenotype = 'rheumatoid_arthritis'
phecodes = c('714.0', '714.1', '714.2', '714.81')
phenoRaw = getPhenoRaw(con, phecodes, gridSet)
write_csv(phenoRaw, gzfile(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype))))

########################################

phenotype = 'multiple_sclerosis'
phecodes = c('340')
phenoRaw = getPhenoRaw(con, phecodes, gridSet)
write_csv(phenoRaw, gzfile(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype))))

########################################

phenotype = 'prostate_cancer'
phecodes = c('185', '233.4', 'V10.46')
phenoRaw = getPhenoRaw(con, phecodes, gridSet)
write_csv(phenoRaw, gzfile(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype))))

########################################

phenotype = 'atrial_fibrillation'
phecodes = c('427.31')
phenoRaw = getPhenoRaw(con, phecodes, gridSet)
write_csv(phenoRaw, gzfile(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype))))

########################################

phenotype = 'alzheimers'
phecodes = c('331.0', '331.00')
phenoRaw = getPhenoRaw(con, phecodes, gridSet)
write_csv(phenoRaw, gzfile(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype))))

########################################

close(con)

########################################

# library('dplyr')
# library('lubridate')

# gridInfo1 = read_csv(file.path(procDir, 'old/grid_info.csv.gz'), col_types='ccDTT')
# colnames(gridInfo1) = tolower(colnames(gridInfo1))
# gridInfo1 = gridInfo1 %>%
# 	rename(gender = gender_epic) %>%
# 	mutate(first_entry_date = as.Date(first_entry_date),
# 			 last_entry_date = as.Date(last_entry_date),
# 			 first_age = time_length(first_entry_date - dob, 'years'),
# 			 last_age = time_length(last_entry_date - dob, 'years')) %>%
# 	arrange(grid)
#
# gridInfo2 = read_csv(file.path(procDir, 'grid_info.csv.gz'), col_types='ccDTT')
# colnames(gridInfo2) = tolower(colnames(gridInfo2))
# gridInfo2 = gridInfo2 %>%
# 	rename(gender = gender_epic) %>%
# 	mutate(first_entry_date = as.Date(first_entry_date),
# 			 last_entry_date = as.Date(last_entry_date),
# 			 first_age = time_length(first_entry_date - dob, 'years'),
# 			 last_age = time_length(last_entry_date - dob, 'years')) %>%
# 	arrange(grid)
#
# print(identical(gridInfo1$grid, gridInfo2$grid))
# print(identical(gridInfo1$first_age, gridInfo2$first_age))
# which(gridInfo1$first_age != gridInfo2$first_age)


# phenotypes = c('alzheimers', 'atrial_fibrillation', 'gout',
# 					'multiple_sclerosis', 'prostate_cancer', 'rheumatoid_arthritis')
#
# for (phenotype in phenotypes) {
# 	phenoRaw1 = read_csv(file.path(procDir, sprintf('old/pheno_%s.csv.gz', phenotype)), col_types='ccT')
# 	colnames(phenoRaw1) = tolower(colnames(phenoRaw1))
# 	phenoRaw1$entry_date = as.Date(phenoRaw1$entry_date)
# 	phenoRaw2 = read_csv(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype)), col_types='ccT')
# 	colnames(phenoRaw2) = tolower(colnames(phenoRaw2))
# 	phenoRaw2$entry_date = as.Date(phenoRaw2$entry_date)
# 	print(identical(phenoRaw1, phenoRaw2))
# }


# for (phenotype in phenotypes) {
# 	phenoRaw = readRDS(file.path(procDir, sprintf('pheno_%s.rds', phenotype)))
# 	write_csv(phenoRaw, gzfile(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype))))
# }
