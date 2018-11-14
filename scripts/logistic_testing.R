library('speedglm')

plinkCovarData = setDT(read_tsv(plinkPaths$covar))
setnames(plinkCovarData, 'FID', 'grid')
plinkCovarData[, IID := NULL]

plinkPhenoData = setDT(read_tsv(plinkPaths$pheno))
setnames(plinkPhenoData, 'FID', 'grid')
plinkPhenoData[, IID := NULL]

plinkGenoData = as.data.table(genoData[, snpData$snpIdx], keep.rownames = TRUE)
setnames(plinkGenoData, 'rn', 'grid')

phecodeStrNow = 'phe290p11'
snpNow = colnames(plinkGenoData)[5]

input = merge(plinkCovarData, plinkPhenoData[, c('grid', phecodeStrNow),
                                             with = FALSE], by = 'grid')
input = merge(input, plinkGenoData[, c('grid', snpNow), with = FALSE],
              by = 'grid')

setnames(input, c(phecodeStrNow, snpNow), c('pheno', 'snp'))
input[, pheno := ifelse(pheno == -9, NA, pheno)]
table(input$pheno, useNA = 'always')

formulaStr = 'pheno ~ snp + rec_len + last_age1 + PC1 + PC2 + sex'
speedglmFit = speedglm(formula(formulaStr), family = binomial(), data = input)
glmFit = glm(formula(formulaStr), family = 'binomial', data = input)

car::vif(glmFit)


