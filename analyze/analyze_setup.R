library('readr')
library('dplyr')
library('tidyr')
library('lubridate')
library('doParallel')
library('snpStats')
library('cowplot')
library('qqman')

eb = element_blank()
theme_set(theme_light() +
			 	theme(axis.text = element_text(color = 'black'), strip.text = element_text(color = 'black'),
			 			panel.grid.minor = eb, legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))


filterInputSex = function(inputDf, whichSex) {
	if (whichSex == 'male') {
		filter(inputDf, sex == 1)
	} else if (whichSex == 'female') {
		filter(inputDf, sex == 2)
	} else {
		inputDf}}


# expected colnames in phenoRaw: grid, entry_date
# expected colnames in gridData: grid, last_age, rec_len
makeGlmInput = function(phenoData, gridData, genoFull, minEvents, whichSex) {
	pheno = phenoData %>%
		count(grid) %>%
		right_join(gridData, by = 'grid') %>%
		mutate(n = ifelse(is.na(n), 0, n)) %>%
		filter(n==0 | n >= minEvents) %>%
		mutate(status = as.integer(n >= minEvents))

	glmInput = tibble(grid = genoFull$fam$pedigree, sex = genoFull$fam$sex) %>%
		inner_join(pheno, by = 'grid') %>%
		select(grid, sex, last_age, rec_len, status)

	glmInput = filterInputSex(glmInput, whichSex)}


# expected colnames in coxInput: grid, sex, last_age, rec_len, status
runGlm = function(glmInput, genoFull, snpName, whichSex, splineDf = 4) {
	glmInput$snp = as(genoFull$genotypes[glmInput$grid, snpName], 'numeric')[,1]
	glmInput = filter(glmInput, !is.na(snp))
	if (whichSex == 'both') {
		glmFit = glm(status ~ snp + sex + rec_len + splines::ns(last_age, df = splineDf),
						 data = glmInput, family = binomial)
	} else {
		glmFit = glm(status ~ snp + rec_len + splines::ns(last_age, df = splineDf),
						 data = glmInput, family = binomial)}}


# expected colnames in phenoRaw: grid, entry_date
# expected colnames in gridData: grid, dob, first_age, last_age
makeCoxInput = function(phenoData, gridData, genoFull, minEvents, buffer, whichSex) {
	phenoCase = phenoData %>%
		semi_join(gridData, by = 'grid') %>%
		group_by(grid) %>%
		arrange(entry_date) %>%
		filter(row_number() == minEvents) %>%
		ungroup()

	phenoControl = gridData %>%
		select(grid) %>%
		anti_join(phenoData, by = 'grid')

	pheno = bind_rows(phenoCase, phenoControl) %>%
		inner_join(gridData, by = 'grid') %>%
		mutate(age = time_length(entry_date - dob, 'years')) %>%
		select(grid, first_age, last_age, age)

	coxInput = tibble(grid = genoFull$fam$pedigree, sex = genoFull$fam$sex) %>%
		inner_join(pheno, by = 'grid') %>%
		group_by(grid) %>%
		mutate(status = as.integer(!is.na(age)),
				 age2 = ifelse(status, age, last_age),
				 age1 = min(first_age, age2 - buffer)) %>%
		ungroup() %>%
		select(grid, sex, age1, age2, status)

	coxInput = filterInputSex(coxInput, whichSex)}


# expected colnames in coxInput: grid, sex, age1, age2, status
runCoxph = function(coxInput, genoFull, snpName, whichSex) {
	coxInput$snp = as(genoFull$genotypes[coxInput$grid, snpName], 'numeric')[,1]
	coxInput = filter(coxInput, !is.na(snp))
	if (whichSex == 'both') {
		coxFit = coxph(Surv(age1, age2, status) ~ snp + sex, data = coxInput)
	} else {
		coxFit = coxph(Surv(age1, age2, status) ~ snp, data = coxInput)}}


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

