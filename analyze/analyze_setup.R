library('data.table')
library('speedglm') # load before dplyr
library('readr')
# library('dplyr')
# library('tidyr')
library('lubridate')
library('doParallel')
library('snpStats')
library('cowplot')
library('qqman')

procDir = 'processed'
resultDir = 'results'

phecodeData = fread(file.path(procDir, 'phecode_definitions1.2.csv'))
phecodeData = phecodeData[,.(phecode = jd_code,
									  phenotype = jd_string,
									  controlExcludeRange = jd_control_exclude_range,
									  whichSex = tolower(ifelse(sex == '', 'both', sex)),
									  rollup, leaf)]

############################################################

eb = element_blank()
theme_set(theme_light() +
			 	theme(axis.text = element_text(color = 'black'), strip.text = element_text(color = 'black'),
			 			panel.grid.minor = eb, legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))


filterInputBySex = function(input, whichSex) {
	if (whichSex == 'male') {
		input[sex == 1,]
		# filter(input, sex == 1)
	} else if (whichSex == 'female') {
		input[sex == 2,]
		# filter(input, sex == 2)
	} else {
		input}}


addSnpToInput = function(input, genoFull, snpName) {
	input[, snp := as(genoFull$genotypes[grid, snpName], 'numeric')[,1]]
	input[!is.na(snp),]}
	# input$snp = as(genoFull$genotypes[input$grid, snpName], 'numeric')[,1]
	# input[!is.na(input$snp),]}


# expected colnames in phenoData: grid, entry_date
# expected colnames in gridData: grid, last_age, rec_len
makeInputGlm = function(phenoData, gridData, genoFull, minEvents, whichSex) {
	# pheno = phenoData %>%
	# 	count(grid) %>%
	# 	right_join(gridData, by = 'grid') %>%
	# 	mutate(n = ifelse(is.na(n), 0, n)) %>%
	# 	filter(n==0 | n >= minEvents) %>%
	# 	mutate(status = as.integer(n >= minEvents))

	pheno = phenoData[, .N, by = grid]
	pheno = merge(pheno, gridData, by = 'grid', all.y = TRUE)
	pheno = pheno[is.na(N) | (N >= minEvents),]
	pheno[, status := as.integer(!is.na(N))]

	input = data.table(grid = genoFull$fam$pedigree, sex = genoFull$fam$sex)
	input = merge(input, pheno, by = 'grid')[, .(grid, sex, last_age, rec_len, status)]
	filterInputBySex(input, whichSex)}
	# glmInput = tibble(grid = genoFull$fam$pedigree, sex = genoFull$fam$sex) %>%
	# 	inner_join(pheno, by = 'grid') %>%
	# 	select(grid, sex, last_age, rec_len, status)
	# filterInputBySex(glmInput, whichSex)}


# expected colnames in input: grid, sex, last_age, rec_len, status
runGlmSexBoth = function(input, splineDf = 4) {
	speedglm(status ~ snp + sex + rec_len + splines::ns(last_age, df = splineDf),
				family = binomial(), data = input)}

runGlmSexOne = function(input, splineDf = 4) {
	speedglm(status ~ snp + rec_len + splines::ns(last_age, df = splineDf),
				family = binomial(), data = input)}


# expected colnames in phenoData: grid, entry_date
# expected colnames in gridData: grid, dob, first_age, last_age
makeInputCox = function(phenoData, gridData, genoFull, minEvents, buffer, whichSex) {
	# phenoCase1 = phenoData %>%
	# 	semi_join(gridData, by = 'grid') %>%
	# 	group_by(grid) %>%
	# 	arrange(entry_date) %>%
	# 	filter(row_number() == minEvents) %>%
	# 	ungroup()
	phenoCase = merge(phenoData, gridData[, .(grid)])
	setkeyv(phenoCase, c('grid', 'entry_date'))
	phenoCase = phenoCase[, if (.N >= minEvents) .SD[minEvents,], by = grid]

	# phenoControl1 = gridData %>%
	# 	select(grid) %>%
	# 	anti_join(phenoData, by = 'grid')
	phenoControl = fsetdiff(gridData[, .(grid)], phenoData[, .(grid)])

	# pheno = bind_rows(phenoCase, phenoControl) %>%
	# 	inner_join(gridData, by = 'grid') %>%
	# 	mutate(age = time_length(entry_date - dob, 'years')) %>%
	# 	select(grid, first_age, last_age, age)
	pheno = rbind(phenoCase, phenoControl, fill = TRUE)
	pheno = merge(pheno, gridData, by = 'grid')
	pheno[, c('status', 'age') := .(as.integer(!is.na(entry_date)),
											  time_length(entry_date - dob, 'years'))]

	input = data.table(grid = genoFull$fam$pedigree, sex = genoFull$fam$sex)
	input = merge(input, pheno[, .(grid, first_age, last_age, age, status)], by = 'grid')
	input[, age2 := ifelse(status, age, last_age)]
	input[, age1 := min(first_age, max(0, age2 - buffer)), by = .(grid, sex, age2, status)]
	input[, .(grid, sex, age1, age2, status)]}
	# coxInput = tibble(grid = genoFull$fam$pedigree, sex = genoFull$fam$sex) %>%
	# 	inner_join(pheno, by = 'grid') %>%
	# 	group_by(grid) %>%
	# 	mutate(status = as.integer(!is.na(age)),
	# 			 age2 = ifelse(status, age, last_age),
	# 			 age1 = min(first_age, max(0, age2 - buffer))) %>%
	# 	ungroup() %>%
	# 	select(grid, sex, age1, age2, status)
	# filterInputBySex(coxInput, whichSex)}


# expected colnames in input: grid, sex, age1, age2, status
runCoxphSexBoth = function(input) {
	coxph(Surv(age1, age2, status) ~ snp + sex, data = input)}

runCoxphSexOne = function(input) {
	coxph(Surv(age1, age2, status) ~ snp, data = input)}


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
