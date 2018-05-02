library('readr')
library('dplyr')
library('lubridate')
library('cowplot')
library('doParallel')
library('snpStats')

eb = element_blank()
theme_set(theme_light() +
			 	theme(axis.text = element_text(color = 'black'), strip.text = element_text(color = 'black'),
			 			panel.grid.minor = eb, legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))


makeSurvPlotDataframe = function(fit, data, variable = NULL) {
	if (is.null(variable)) {
		pred = survexp(~ 1, data = data, ratetable = fit)
		curve = data.frame(time = c(0, pred$time), surv = c(1, pred$surv))

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

		# remove leading zeros, while survexp returns non monotonic results
		if (length(dim(pred$surv)) == 2) {
			for (ii in 1:ncol(pred$surv)) {
				for (jj in nrow(pred$surv):2) {
					if (pred$surv[jj, ii] > pred$surv[jj - 1, ii]) {
						pred$surv[jj - 1, ii] = 1 }}}}

		curve = data.frame(time = rep(c(0, pred$time), length(lev)), surv = c(rbind(1, pred$surv)),
								 value = rep(lev, each = 1 + length(pred$time)))}}

########################################

load('pilot_workspace.Rdata')

exampleDf = tibble(phenotype = c('alzheimers', 'atrial_fibrillation', 'multiple_sclerosis', 'prostate_cancer'),
						 snp = c('rs769449', 'rs6843082', 'rs9271366', 'rs10993994')) %>%
	inner_join(phenoMetadata, by = 'phenotype')

registerDoParallel(cores = 1)

d = foreach(ii = 1:nrow(exampleDf), .combine = rbind) %dopar% {
	phenotype = exampleDf$phenotype[ii]
	snp = exampleDf$snp[ii]
	whichSex = exampleDf$whichSex[ii]

	phenoRaw = read_csv(file.path(procDir, sprintf('pheno_%s.csv.gz', phenotype)), col_types = 'ccT')
	colnames(phenoRaw) = tolower(colnames(phenoRaw))
	phenoRaw$entry_date = as.Date(phenoRaw$entry_date)

	coxInput = makeCoxInput(phenoRaw, gridInfo, exome, minEvents, buffer, whichSex)
	coxInput$snp = as(exome$genotypes[coxInput$grid, snp], 'numeric')[,1]
	coxInput = filter(coxInput, !is.na(snp))
	coxFit = runCoxph(coxInput, exome, snp, whichSex)

	d = makeSurvPlotDataframe(coxFit, coxInput, 'snp')
	d = d %>%
		mutate(phenotype = phenotype,
				 snp = snp)}

########################################

d = mutate(d, phenotypeSnp = sprintf('%s, %s', phenotype, snp))

p = ggplot(d) +
	facet_wrap(~ phenotypeSnp, nrow = 2, scales = 'free') +
	geom_step(aes(x = time, y = 1 - surv, color = factor(2 - value, levels = 2:0)), size = 1) +
	scale_color_manual(values = c('#225ea8', '#41b6c4', '#a1dab4')) +
	labs(x = 'Age (y)', y = 'Adjusted cumulative hazard', color = 'Minor allele\ncount')

ggsave(file.path(resultDir, 'pilot_examples.pdf'), plot = p, width = 6.5, height = 5)
