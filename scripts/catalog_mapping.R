library('data.table')
library('readr')

# criteria:
# disease_trait corresponds to phecode (e.g., not severity, survival, or interaction)
# >=1000 European ancestry cases
# >=1000 European ancestry controls

procDir = 'processed'
studyFilename = 'gwas_catalog_v1.0.2-studies_r2018-08-28.tsv'

studyData = read_tsv(file.path(procDir, 'gwas_catalog', studyFilename))
setDT(studyData)
colnames(studyData) = gsub('/|\\[|\\]|\\s', '_', tolower(colnames(studyData)))

phecodeData = read_csv(file.path(procDir, 'phecode_data.csv.gz'),
                       col_types = 'ccc??????')
setDT(phecodeData)

############################################################

ignoredWords = c('and', 'or', 'in', 'of', 'to', 'the', 'on', 'with', 'as',
                 'but', 'not', 'origin', 'disease', 'diseases', 'disorder',
                 'disorders', 'syndrome','signs', 'symptoms', 'nos', 'other',
                 'due', 'use')
gsubPattern = sprintf('(\\b%s\\b)', paste(ignoredWords, collapse = '\\b|\\b'))

# leave as is: ' -
phecodeData[, pheno_clean := gsub(gsubPattern, '', tolower(phenotype))]
phecodeData[, pheno_clean := gsub('\\,|\\;|\\(|\\)|\\[|\\]|\\.|\\&', '', pheno_clean)]
phecodeData[, pheno_clean := gsub('(^| )[a-z]( |$)', '', pheno_clean)]

# make this "and" instead of "or"?
phecodeData[, pattern := paste(strsplit(pheno_clean, '\\s+')[[1]],
                               collapse = '\\b|\\b'),
            by = 1:nrow(phecodeData)]
phecodeData[, pattern := sprintf('(?:\\b%s\\b)', pattern)]
setcolorder(phecodeData, c('pattern', 'pheno_clean', 'phenotype'))

############################################################

# phecodesTest = read_tsv(file.path('params/mega/phecodes_test.tsv'),
#                         col_names = FALSE, col_types = 'c')$X1
# phecodeDataKeep = phecodeData[phecode %in% phecodesTest]

# phecodeSummary = read_tsv(file.path('processed/mega/phecode_summary_in_catalog.tsv'))
phecodeSummary = read_tsv(file.path('processed/mega/phecode_summary.tsv'))
setDT(phecodeSummary)
phecodeDataKeep = merge(phecodeData, phecodeSummary, by = 'phecode')
phecodeDataKeep = phecodeDataKeep[nCases >= 100]

############################################################

getCartesianProd = function(d1, d2) {
  d1New = d1
  d1New[, dummy := 1]
  d2New = d2
  d2New[, dummy := 1]
  d = merge(d1New, d2New, by = 'dummy', allow.cartesian = TRUE)
  d[, dummy := NULL]
  return(d)}

sdn = unique(studyData[(association_count > 0) &
                         grepl('\\d{1}\\,\\d{3}', initial_sample_size) &
                         grepl('european', initial_sample_size, ignore.case = TRUE) &
                         !grepl(', ', mapped_trait_uri),
                       .(study_accession, pubmedid, date, initial_sample_size,
                         study, disease_trait)])

# a1 = getCartesianProd(sdn, phecodeDataKeep[, .(pattern, phenotype, phecode)])
#
# a2 = a1[, if (grepl(pattern, disease_trait, ignore.case = TRUE)) .SD,
#         by = 1:nrow(a1)]
#
# a3 = a2[phecode == '250.2']
# a3 = a2[phecode == phecodesTest[1]]

# a3 = a2[grepl('gout|uric|urate', disease_trait, ignore.case = TRUE)]
# a3 = a2[grepl('multiple sclerosis', disease_trait, ignore.case = TRUE)]
# a3 = a2[grepl('prostate cancer', disease_trait, ignore.case = TRUE)]

#################################

phecodesDone = read_csv('processed/gwas_catalog/catalog_phecodes_done.csv', col_types = 'c')$phecode
pheDone = phecodeDataKeep[phecode %in% phecodesDone]
pheRemaining = phecodeDataKeep[!(phecode %in% phecodesDone)]

a = sdn[grepl('urate', disease_trait, ignore.case = TRUE)]
          # !grepl('paediatric', disease_trait, ignore.case = TRUE)]

a = sdn[tolower(disease_trait) %in% c('urate levels', 'gout')]

a1 = a[!(study_accession %in% c('GCST000130')),
       .(phecode = '274.1', study_accession)]
write_csv(a1, 'catalog_study_now.csv')

# a2 = phecodeData[phecode %in% pheDone]
# a3 = phecodeDataKeep[grepl('breast', phenotype, ignore.case = TRUE)]
