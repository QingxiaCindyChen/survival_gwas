library('tidyverse')
library('cowplot')

procDir = 'processed'

mapUnrolled = read_csv(file.path(procDir, 'phecode_icd9_map_unrolled.csv'), col_types='cc')
mapRolled = read_csv(file.path(procDir, 'phecode_icd9_rolled.csv'), col_types='cccccciii')
phecodeDef = read_csv(file.path(procDir, 'phecode_definitions1.2.csv'), col_types='ccccii')

pilotPhenotypes = read_tsv(file.path(procDir, 'pilot_phenotypes.tsv'), col_types='cc')

a = pilotPhenotypes %>%
	inner_join(mapUnrolled, by='phecode') %>%
	write_tsv('pilot_icd9.tsv')

# a = phecodeDef %>%
# 	filter(str_detect(tolower(jd_string), 'atrial fibrillation'))

# a1 = mapRolled %>%
# 	filter(PheCode=='427.21')
	# filter(str_detect(tolower(Phenotype), 'gout'))

# a2 = mapUnrolled %>%
# 	filter(phecode=='427.21')
