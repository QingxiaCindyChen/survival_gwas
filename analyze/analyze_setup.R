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
phecodeData = phecodeData[,.(phecode = jd_code,
									  phenotype = jd_string,
									  controlExcludeRange = jd_control_exclude_range,
									  whichSex = tolower(ifelse(sex == '' | is.na(sex), 'both', sex)),
									  rollup, leaf)]

############################################################

eb = element_blank()
theme_set(theme_light() +
			 	theme(axis.text = element_text(color = 'black'), strip.text = element_text(color = 'black'),
			 			panel.grid.minor = eb, legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))


filterInputBySex = function(input, whichSex) {
	if (whichSex == 'male') {
		input[sex == 1,]
	} else if (whichSex == 'female') {
		input[sex == 2,]
	} else {
		input}}


# expected colnames in phenoData: grid, age
# expected colnames in gridData: grid, first_age, last_age
makeInput = function(phenoData, gridData, genoFull, whichSex, minEvents, buffer) {
	phenoCase = merge(phenoData, gridData[, .(grid)], by = 'grid')
	setkeyv(phenoCase, c('grid', 'age'))
	phenoCase = phenoCase[, if (.N >= minEvents) .SD[minEvents,], by = grid]
	phenoControl = fsetdiff(gridData[, .(grid)], phenoData[, .(grid)])

	pheno = rbind(phenoCase, phenoControl, fill = TRUE)
	pheno = merge(pheno, gridData, by = 'grid')

	input = data.table(grid = genoFull$fam$pedigree, sex = genoFull$fam$sex)
	input = merge(input, pheno, by = 'grid')
	input[, status := ifelse(is.na(age), 0, 1)]
	input[, age2 := ifelse(status, age, last_age)]
	input[, age1 := min(first_age, max(0, age2 - buffer)), by = grid]
	input}


addSnpToInput = function(input, genoFull, snpName) {
	input[, snp := as(genoFull$genotypes[grid, snpName], 'numeric')[,1]]
	input[!is.na(snp),]}


getGlmStr = function(whichSex = 'both', nPC = 3, splineDf = 4) {
	formStr = sprintf('status ~ snp + rec_len + splines::ns(last_age, df = %d)', splineDf)
	if (whichSex == 'both') {
		formStr = paste(formStr, '+ sex')}
	if (nPC > 0) {
		formStr = paste(formStr, '+', paste0('PC', 1:nPC, collapse = ' + '))}
	formStr}


getCoxStr = function(whichSex = 'both', nPC = 3) {
	formStr = 'Surv(age1, age2, status) ~ snp'
	if (whichSex == 'both') {
		formStr = paste(formStr, '+ sex')}
	if (nPC > 0) {
		formStr = paste(formStr, '+', paste0('PC', 1:nPC, collapse = ' + '))}
	formStr}


runGlm = function(formulaStr, input) {
	speedglm(formula(formulaStr), family = binomial(), data = input)}


runCox = function(formulaStr, input) {
	coxph(formula(formulaStr), data = input)}


makeSurvPlotDataframe = function(fit, data, variable = NULL) {
	# adapted from survminer
	if (is.null(variable)) {
		pred = survexp(~ 1, data = data, ratetable = fit)
		d = data.frame(time = c(0, pred$time), surv = c(1, pred$surv))

	} else {
		lev = sort(unique(data[[variable]]))
		df0 = data
		df0[[variable]] = '_reference_' # collision unlikely, but not impossible
		form = update(formula(fit), as.formula(sprintf('%s ~ . - %s', variable, variable)))
		environment(form) = parent.frame()

		rwt = numeric(nrow(data))
		for (level in lev) {
			idx = which(data[[variable]] == level)
			if (length(idx) > 0) {
				df1 = data[idx,, drop = FALSE]
				ndf = rbind(df0, df1)
				ndf[[variable]] = factor(ndf[[variable]])
				model = glm(formula = form, data = ndf, family = binomial)
				allRes = predict(model, newdata = data, type = 'response')
				rwt[idx] = 1/allRes[idx]}}

		nform = as.formula(sprintf('%s ~ %s', as.character(formula(fit))[2], variable))
		nfit = coxph(formula = nform, data = data, weights = rwt)
		pred = survexp(as.formula(sprintf('~ %s', variable)), data = data, ratetable = nfit)

		if (length(dim(pred$surv)) == 2) {
			for (ii in 1:ncol(pred$surv)) {
				for (jj in nrow(pred$surv):2) {
					if (pred$surv[jj, ii] > pred$surv[jj - 1, ii]) {
						pred$surv[jj - 1, ii] = 1 }}}}

		d = data.frame(time = rep(c(0, pred$time), length(lev)), surv = c(rbind(1, pred$surv)),
							value = rep(lev, each = 1 + length(pred$time)))}}
