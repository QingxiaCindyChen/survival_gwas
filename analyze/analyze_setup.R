library('data.table')
library('speedglm')
library('readr')
library('lubridate')
library('doParallel')
library('snpStats')
library('cowplot')
library('qqman')

procDir = 'processed'
resultDir = 'results'

phecodeData = setDT(read_csv(file.path(procDir, 'phecode_definitions1.2.csv'), col_types = 'ccc???'))
phecodeData = phecodeData[, .(phecode = jd_code,
                              phenotype = jd_string,
                              controlExcludeRange = jd_control_exclude_range,
                              whichSex = tolower(ifelse(sex == '' | is.na(sex), 'both', sex)),
                              rollup, leaf)]

############################################################

eb = element_blank()
theme_set(theme_light() +
            theme(axis.text = element_text(color = 'black'), strip.text = element_text(color = 'black'),
                  panel.grid.minor = eb, legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))


# expected colnames in phenoData: grid, age
# expected colnames in gridData: grid, first_age, last_age
makeInput = function(phenoData, gridData, whichSex, minEvents, buffer) {
  phenoCase = phenoData
  setkeyv(phenoCase, c('grid', 'age'))
  phenoCase = phenoCase[, if (.N >= minEvents) .SD[minEvents,], by = grid]
  phenoControl = fsetdiff(gridData[, .(grid)], phenoData[, .(grid)])

  input = rbind(phenoCase, phenoControl, fill = TRUE)
  input = merge(input, gridData, by = 'grid')

  if (whichSex == 'male') {
    input = input[sex == 1]
  } else if (whichSex == 'female') {
    input = input[sex == 2]}

  input[, status := ifelse(is.na(age), 0, 1)]
  input[, age2 := ifelse(status, age, last_age)]
  input[, age1 := min(first_age, max(0, age2 - buffer)), by = grid]
  input}


addSnpToInput = function(input, genoFull, snp) {
  input[, genotype := as(genoFull$genotypes[grid, snp], 'numeric')[,1]]
  input[!is.na(genotype),]}


makeCoxStr = function(whichSex, nPC) {
  formStr = 'Surv(age1, age2, status) ~ genotype'
  if (whichSex == 'both') {
    formStr = paste(formStr, '+ sex')}
  if (nPC > 0) {
    formStr = paste(formStr, '+', paste0('PC', 1:nPC, collapse = ' + '))}
  formStr}


runCox = function(formulaStr, input) {
  coxph(formula(formulaStr), data = input)}


makeGlmStr = function(whichSex, nPC, splineDf) {
  formStr = sprintf('status ~ genotype + rec_len + %s', paste0('last_age', 1:splineDf, collapse = ' + '))
  if (whichSex == 'both') {
    formStr = paste(formStr, '+ sex')}
  if (nPC > 0) {
    formStr = paste(formStr, '+', paste0('PC', 1:nPC, collapse = ' + '))}
  formStr}


runGlm = function(formulaStr, input) {
  speedglm(formula(formulaStr), family = binomial(), data = input)}
