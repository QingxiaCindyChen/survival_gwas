library('data.table')
library('readr')

procDir = 'processed'

############################################################

catalogAssocRaw = read_tsv(file.path(procDir, 'gwas_catalog_v1.0.2-associations_e92_r2018-05-29.tsv'),
                           guess_max = 1e5)
setDT(catalogAssocRaw)
colnames(catalogAssocRaw) = tolower(colnames(catalogAssocRaw))

############################################################

minEvents = 2
phenoData = read_csv(file.path(procDir, 'exome_phenotype_data.csv.gz'), col_types = 'ccD')
setDT(phenoData)
phenoTmp = phenoData[, .N, by = .(grid, phecode)
                     ][, .(nCasesApprox = sum(N >= minEvents)), by = phecode]

phecodeData = setDT(read_csv(file.path(procDir, 'phecode_data.csv.gz'), col_types = 'ccc??????'))
phecodeData = merge(phecodeData, phenoTmp, by = 'phecode')[order(nCasesApprox, decreasing = TRUE)]
setcolorder(phecodeData, c('phecode', 'phenotype', 'nCasesApprox'))

############################################################

traitData = unique(catalogAssocRaw[, .(mapped_trait)])

phecodeTraitMapping = data.table(phecode = c('185', '274.1', '290.11', '335', '427.21', '714.1'),
                                 mapped_trait = c('prostate carcinoma', 'gout', 'Alzheimers disease',
                                                  'multiple sclerosis', 'atrial fibrillation',
                                                  'rheumatoid arthritis'))

############################################################

pvalCutoff = 5e-8

# https://www.ebi.ac.uk/gwas/docs/fileheaders#_file_headers_for_catalog_version_1_0_2
catalogAssoc = merge(catalogAssocRaw, phecodeTraitMapping, by = 'mapped_trait')
catalogAssoc[, assoc_id := 1:nrow(catalogAssoc)]
catalogAssoc = catalogAssoc[pvalue_mlog >= -log10(pvalCutoff)]
# possibly filter for more stuff, e.g., ethnicity of cohort

catalogAssocTmp = catalogAssoc[, strsplit(`strongest snp-risk allele`, '; ', fixed = TRUE), by = assoc_id
                               ][, c('snp', 'risk_allele') := tstrsplit(V1, '-', fixed = TRUE)]
catalogAssoc = merge(catalogAssoc, catalogAssocTmp[, .(assoc_id, snp, risk_allele)], by = 'assoc_id')

catalogAssocUnique = unique(catalogAssoc[, .(phecode, snp)])
write_csv(catalogAssocUnique, gzfile(file.path(procDir, 'catalog_associations_unique.csv.gz')))
