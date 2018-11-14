library('readr')
library('data.table')

efoPhecodeMapping = read_csv('GWASCatalogMapping/EFOMapping/MapTables/AllEFO_toPhecodes.csv',
                             col_types = 'cccccccc')
colnames(efoPhecodeMapping) = tolower(colnames(efoPhecodeMapping))
setDT(efoPhecodeMapping)

####################

sdn = unique(studyData[(association_count > 0) &
                         grepl('\\d{1}\\,\\d{3}', initial_sample_size) &
                         grepl('european', initial_sample_size, ignore.case = TRUE) &
                         !grepl(', ', mapped_trait_uri),
                       .(study_accession, pubmedid, date, initial_sample_size,
                         study, disease_trait, mapped_trait_uri)])

sdn[, efo := tstrsplit(mapped_trait_uri, split = 'http://www.ebi.ac.uk/efo/',
                       fixed = TRUE, keep = 2)]

####################

a = merge(phecodeDataKeep, efoPhecodeMapping[, .(efo, phecode)], by = 'phecode')
a = merge(a, sdn[, .(efo)], by = 'efo')
a = unique(a[, .(phecode)])
phecodeDataKeep1 = merge(phecodeDataKeep, a, by = 'phecode')

write_tsv(phecodeDataKeep1[, .(phecode, phenotype, nCases)], 'phecode_summary_in_catalog.tsv')
