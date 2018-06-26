#Updated GWAS catalog v1.0.2 - includes a snomed icd match too
#Inputs - GWASCatalog_ICD9_OMIM_SNOMEDCT_Map.csv, phecode_icd9_rolled.csv, ExomeIDs.csv, snomed_ICD9, gwas_catalog-ancestry.tsv
#Outputs - MultiPheROC_ICD_GWASv102.pdf, MultiPheROC_SNO_GWASv102.pdf, MultiPheROC_ICD_SNO_GWASv102.pdf
library('readr')
library('cowplot')
library('data.table')
library('survival')
library('snpStats')
library('survminer')
library('dplyr')
library('plotROC')
options(warn = -1)

filepath = file.path('/Users', 'srhoades', 'Dropbox (VUMC)/', 'JakeColab_CountingGWAS', 'counting_gwas', 'analyze')
#ROC Curve setup/functions
source(file.path(filepath, 'GWASCatalog', 'ROC_setup.R'))

#File names
Coxfiles = list.files(file.path(filepath, 'GWASCatalog', 'exome_full'))[c(grep('_cox.tsv.gz', list.files(file.path(filepath, 'GWASCatalog', 'exome_full'))))]
Logfiles = list.files(file.path(filepath, 'GWASCatalog', 'exome_full'))[c(grep('_logistic.tsv.gz', list.files(file.path(filepath, 'GWASCatalog', 'exome_full'))))]

#Get file for phecode of interest, and make sure its matched between cox and logistic
FilePhecodes = gsub('exome_|_cox.tsv.gz', '', Coxfiles)

#New with omim and snomedct
GWASCat = as.data.frame(read_csv(file = file.path(filepath, 'GWASCatalog', 'GWASCatalog_ICD9_OMIM_SNOMEDCT_Map.csv')))
GWASCatInternal = setDT(read_delim(file = file.path(filepath, 'GWASCatalog', 'catalog_phecode_map.txt'), delim = '\t', col_types = cols()))

#Exome
SnpIDs = setDT(read_csv(file = file.path(filepath, 'GWASCatalog', 'ExomeIDs.csv'), col_types = cols()))

GWASCatInternal = GWASCatInternal[c(GWASCatInternal$SNP %in% SnpIDs$rsID), ]

#Phecodes from the internal (Lisa's) list
MatchPhes = gsub('\\.', 'p', GWASCatInternal$jd_code)
MatchPhes = paste0('phe', MatchPhes)
InternalFilePhecodes = FilePhecodes[FilePhecodes %in% MatchPhes] #82

#v1.0.2
#Don't need OMIM here - match p-value
GWASCat = GWASCat %>% select(-OMIM) %>% unique()
GWASCat = GWASCat[GWASCat$`P-VALUE` < 5e-08, ]

#SNOMED
snoIDs = setDT(read_csv(file = file.path(filepath, 'GWASCatalog', 'snomed_ICD9'), col_types = 'cc'))

#Ancestry filter
Anc = as.data.frame(read_tsv(file = file.path(filepath, 'GWASCatalog', 'gwas_catalog-ancestry.tsv')))
#Messy names, just anything with European in it? many are mixed
AncList = as.character(Anc$PUBMEDID[c(grep('Euro', Anc$`BROAD ANCESTRAL CATEGORY`))])

phemap = read_csv(file = file.path(filepath, 'GWASCatalog', 'phecode_icd9_rolled.csv'))

#Trim rows without an ICD9
GWASCatICD = GWASCat[!(is.na(GWASCat$ICD)), ] #roughly half of the rows removed, but only 180 unique ICD9s here

#Keep NAs - since we want the rows which have SNOMED but not ICD (if we think of ICD as gold standard)
GWASCatSNO = GWASCat[(is.na(GWASCat$ICD)), ] #roughly half of the rows removed, but only 180 unique ICD9s here

#European participants
GWASCatICD = GWASCatICD[c(GWASCatICD$PUBMEDID %in% AncList), ] 
GWASCatSNO = GWASCatSNO[c(GWASCatSNO$PUBMEDID %in% AncList), ]

#Snp format
GWASCatICD = data.table(GWASCatICD)[, list(SNPS = unlist(strsplit(SNPS, '; '))), by=c('P-VALUE', 'ICD', 'PUBMEDID')] #Remove duplicates from the same study
GWASCatICD = GWASCatICD %>% unique()# %>% select(-PUBMEDID) - keep for later to make sure you aren't duplicating the ICD-SNOMED DFs after phecode mapping
GWASCatICD = GWASCatICD[c(GWASCatICD$SNP %in% SnpIDs$rsID), ] #Exome
GWASCatICD$ICD = as.character(GWASCatICD$ICD)

GWASCatSNO = data.table(GWASCatSNO)[, list(SNPS = unlist(strsplit(SNPS, '; '))), by=c('P-VALUE', 'SNOMEDCT', 'PUBMEDID')]
GWASCatSNO = GWASCatSNO %>% unique()# %>% select(-PUBMEDID)
GWASCatSNO = GWASCatSNO[c(GWASCatSNO$SNP %in% SnpIDs$rsID), ] #Exome
GWASCatSNO$SNOMEDCT = as.character(GWASCatSNO$SNOMEDCT)

#Catch formatting
#if length of string before period is 2, add a 0 in front, if the length is 1, add two zeros in front
for(phe in 1:nrow(GWASCatICD)){
  if(nchar(gsub('\\.\\d?\\d?', '', GWASCatICD$ICD[phe]))==2){
    GWASCatICD$ICD[phe] = paste0('0', GWASCatICD$ICD[phe])
  }
  if(nchar(gsub('\\.\\d?\\d?', '', GWASCatICD$ICD[phe]))==1){
    GWASCatICD$ICD[phe] = paste0('00', GWASCatICD$ICD[phe])
  }
} 

GWASCatICD = inner_join(GWASCatICD, phemap, by = c('ICD' = 'ICD9')) #089.52 ICD omitted
GWASCatSNO = inner_join(GWASCatSNO, snoIDs, by = c('SNOMEDCT' = 'SNOMED_CID'))
GWASCatSNO = inner_join(GWASCatSNO, phemap, by = c('ICD_CODE' = 'ICD9')) #Ranges and "V" phecodes omitted

GWASCatICD = GWASCatICD %>% rename(SNP = SNPS, jd_code = PheCode) %>% dplyr::select(SNP, jd_code, PUBMEDID)
GWASCatSNO = GWASCatSNO %>% rename(SNP = SNPS, jd_code = PheCode) %>% dplyr::select(SNP, jd_code, PUBMEDID)

GWASCatSNO = GWASCatSNO[!(is.na(GWASCatSNO$jd_code)), ]
GWASCatICD = GWASCatICD[!(is.na(GWASCatICD$jd_code)), ]
GWASBoth = rbind(GWASCatICD, GWASCatSNO) %>% unique() %>% select(-PUBMEDID)
GWASCatSNO = select(GWASCatSNO, -PUBMEDID)
GWASCatICD = select(GWASCatICD, -PUBMEDID)

MatchPhestestICD = unique(gsub('\\.', 'p', GWASCatICD$jd_code))
for(item in 1:length(MatchPhestestICD)){
  #Matches to phecode filenames require removing starting zeros
  if(substr(MatchPhestestICD[item], start = 1, stop = 1)=='0'){
    MatchPhestestICD[item] = substr(MatchPhestestICD[item], start = 2, stop = nchar(MatchPhestestICD[item]))
  }
}

MatchPhestestSNO = unique(gsub('\\.', 'p', GWASCatSNO$jd_code))
for(item in 1:length(MatchPhestestSNO)){
  if(substr(MatchPhestestSNO[item], start = 1, stop = 1)=='0'){
    MatchPhestestSNO[item] = substr(MatchPhestestSNO[item], start = 2, stop = nchar(MatchPhestestSNO[item]))
  }
}

MatchPhestestJoin = unique(gsub('\\.', 'p', GWASBoth$jd_code))
for(item in 1:length(MatchPhestestJoin)){
  if(substr(MatchPhestestJoin[item], start = 1, stop = 1)=='0'){
    MatchPhestestJoin[item] = substr(MatchPhestestJoin[item], start = 2, stop = nchar(MatchPhestestJoin[item]))
  }
}

MatchPhestestICD = paste0('phe', MatchPhestestICD)
MatchPhestestSNO = paste0('phe', MatchPhestestSNO)
MatchPhestestJoin = paste0('phe', MatchPhestestJoin)

#Lists of phecodes to compile/plot
FilePhecodesICD = FilePhecodes[FilePhecodes %in% MatchPhestestICD] 
FilePhecodesSNO = FilePhecodes[FilePhecodes %in% MatchPhestestSNO] 
FilePhecodesBoth = FilePhecodes[FilePhecodes %in% MatchPhestestJoin] 
FilePhecodesBoth = unique(FilePhecodesBoth)

#Build ROC DFs
InternalROC = PhePhetchROC(InternalFilePhecodes, GWASCatInternal)
ICDROC = PhePhetchROC(FilePhecodesICD, GWASCatICD)
SNOROC = PhePhetchROC(FilePhecodesSNO, GWASCatSNO)
BothROC = PhePhetchROC(FilePhecodesBoth, GWASBoth)

#GGplots
OutputInternal = BuildROCCurve(InternalROC, 'GWAS Catalog - Internal Mapping')
OutputICD = BuildROCCurve(ICDROC, 'GWAS Catalog v1.0.2 - ICD Map')
OutputSNO = BuildROCCurve(SNOROC, 'GWAS Catalog v1.0.2 - SNOMEDCT Map')
OutputBoth = BuildROCCurve(BothROC, 'GWAS Catalog v1.0.2 - ICD and SNOMEDCT Map')

ggsave(filename = 'MultiPheROC_ICD_InternalMap.pdf', path = file.path(filepath, 'GWASCatalog'), plot = OutputInternal, width = 5, height = 3.5, units = c('in'))
ggsave(filename = 'MultiPheROC_ICD_InternalMap.png', path = file.path(filepath, 'GWASCatalog'), plot = OutputInternal, width = 5, height = 3.5, units = c('in'), dpi = 300)
ggsave(filename = 'MultiPheROC_ICD_GWASv102.pdf', path = file.path(filepath, 'GWASCatalog'), plot = OutputICD, width = 5, height = 3.5, units = c('in'))
ggsave(filename = 'MultiPheROC_ICD_GWASv102.png', path = file.path(filepath, 'GWASCatalog'), plot = OutputICD, width = 5, height = 3.5, units = c('in'), dpi = 300)
ggsave(filename = 'MultiPheROC_SNO_GWASv102.pdf', path = file.path(filepath, 'GWASCatalog'), plot = OutputSNO, width = 5, height = 3.5, units = c('in'))
ggsave(filename = 'MultiPheROC_SNO_GWASv102.png', path = file.path(filepath, 'GWASCatalog'), plot = OutputSNO, width = 5, height = 3.5, units = c('in'), dpi = 300)
ggsave(filename = 'MultiPheROC_ICD_SNO_GWASv102.pdf', path = file.path(filepath, 'GWASCatalog'), plot = OutputBoth, width = 5, height = 3.5, units = c('in'))
ggsave(filename = 'MultiPheROC_ICD_SNO_GWASv102.png', path = file.path(filepath, 'GWASCatalog'), plot = OutputBoth, width = 5, height = 3.5, units = c('in'), dpi = 300)
