library('data.table')
library('readr')

procDir = 'processed'

filePath = file.path(procDir, 'phecode_info_1.2.rda')
download.file('https://github.com/PheWAS/PheWAS/raw/068aa4164deeefe3cd7c19da297c869e2bc197fa/data/pheinfo.rda',
              filePath)
load(filePath) # dataframe called pheinfo
setDT(pheinfo)
pheinfo[, description := NULL]

phecodeData = setDT(read_csv(file.path(procDir, 'phecode_definitions1.2.csv'), col_types = 'ccc???'))
phecodeData = phecodeData[, .(phecode = jd_code,
                              phenotype = jd_string,
                              controlExcludeRange = jd_control_exclude_range,
                              whichSex = tolower(ifelse(sex == '' | is.na(sex), 'both', sex)),
                              rollup, leaf)]

phecodeData = merge(phecodeData, pheinfo, by = 'phecode', all.x = TRUE)
write_csv(phecodeData, gzfile(file.path(procDir, 'phecode_data.csv.gz')))
