library('data.table')
library('readr')

procDir = 'processed'

studyFilename = 'gwas_catalog_v1.0.2-studies_r2018-08-28.tsv'
studyData = read_tsv(file.path(procDir, 'gwas_catalog', studyFilename))
setDT(studyData)
colnames(studyData) = gsub('/|\\[|\\]|\\s', '_', tolower(colnames(studyData)))

phecodeData = read_csv(file.path(procDir, 'phecode_data.csv.gz'),
                       col_types = 'ccc??????')
setDT(phecodeData)

studyPhecodeMapping = read_csv(file.path(procDir, 'gwas_catalog',
                                         'catalog_study_phecode.csv'), col_types = 'cc')
setDT(studyPhecodeMapping)

d = merge(studyPhecodeMapping, phecodeData[, .(phecode, phenotype)],
          by = 'phecode')
d = merge(d, studyData[, .(study_accession, disease_trait, initial_sample_size, study)],
          by = 'study_accession')

dSummary = d[, .(nStudies = .N), by = .(phecode, phenotype, disease_trait)]
setnames(dSummary, 'phenotype', 'phecode_description')

dPhecode = unique(d[order(phecode), .(phecode)])
write_tsv(dPhecode, file.path(procDir, 'gwas_catalog', 'phecodes_catalog50.tsv'),
          col_names = FALSE)
