source(file.path('analyze', 'analyze_setup.R'))

analysisDir = 'exome_test'
filePrefix = 'exome'

############################################################
# ensures that plotting uses exactly the same data used in the analysis,
# including gridData, phenoData, snpData, and all parameters and functions

load(file.path(resultDir, analysisDir, sprintf('%s_workspace.Rdata', filePrefix)))
resultDir = 'results'

genoData = readRDS(file.path(procDir, sprintf('%s_genotype_data.rds', filePrefix)))

snpData = merge(snpData, data.table(snp = rownames(genoData$genoSummary), genoData$genoSummary), by = 'snp')

# phecodeDataKeep = phecodeDataKeep[filename %in% list.files(file.path(resultDir, analysisDir))]
phecodeDataKeep = phecodeDataKeep[c(1:2, 4:6)]

############################################################

registerDoParallel(cores = 1)

maxPvalLoad = 1

resultList = foreach(ii = 1:nrow(phecodeDataKeep)) %dopar% {
	coxFilepath = file.path(resultDir, analysisDir, phecodeDataKeep$coxFilename[ii])
	dCox = setDT(read_tsv(coxFilepath, col_types = '????c'))
	dCox[, beta := -beta] # coefficients from coxph are for the major allel
	dCox[, method := 'cox']

	logisticFilepath = file.path(resultDir, analysisDir, phecodeDataKeep$logisticFilename[ii])
	dLogistic = setDT(read_tsv(logisticFilepath, col_types = '?c?cc???????'))
	colnames(dLogistic) = tolower(colnames(dLogistic))
	dLogistic = dLogistic[, .(beta, se, z = stat, pval = p, snp = snp, method = 'logistic')]
	# TODO: verify that the stat value from plink's logistic regression is a z-value

	d = rbind(dCox, dLogistic)
	d[, phecode := phecodeDataKeep$phecode[ii]]
	dLambda = d[, .(lambdaMed = median(z^2) / qchisq(0.5, 1)), by = .(phecode, method)]

	d = d[, if (mean(log(pval), na.rm = TRUE) <= log(maxPvalLoad)) .SD, by = snp]
	list(d, dLambda)
}

result = rbindlist(lapply(resultList, function(d) d[[1]]), use.names = TRUE)
result = merge(result, snpData, by = 'snp', sort = FALSE)
result = merge(result, phecodeData[, .(phecode, phenotype)], by = 'phecode')

resultLambda = rbindlist(lapply(resultList, function(d) d[[2]]), use.names = TRUE)

############################################################

# phecodeNow = '250.1'
# snpNow = 'rs9275495'

# phecodeNow = '250.2'
# snpNow = 'rs7903146'

phecodeNow = '335'
snpNow = 'rs3130683'

# phecodeNow = '714.1'
# snpNow = 'rs7774434'

# phecodeNow = '185'
# snpNow = 'rs10505477'

whichSexNow = phecodeData[phecode == phecodeNow, whichSex]
phenoDataNow = phenoData[phecode == phecodeNow, .(grid, age)]

inputBase = makeInput(phenoDataNow, gridData, whichSexNow, minEvents, buffer)
inputNow = addSnpToInput(inputBase, genoData$genoFull, snpNow)

coxStr = makeCoxStr(whichSexNow, 0) # can't use PCs to make empirical survival curves
a = survfit(as.formula(coxStr), data = inputNow)

############################################################

aDt = data.table(time = a$time, surv = a$surv, strata = names(a$strata)[rep(1:length(a$strata), a$strata)])
if (whichSexNow == 'both') {
	aStrata = data.table(strata = names(a$strata), snp = rep(0:2, each = 2), sex = rep(c('male', 'female'), 3))
} else {
	aStrata = data.table(strata = names(a$strata), snp = 0:2)
}
aDt = merge(aDt, aStrata, by = 'strata')

if (whichSexNow == 'both') {
	p = ggplot(aDt) +
		geom_step(aes(x = time, y = surv, color = factor(2 - snp), size = sex)) +
		labs(x = 'Age (y)', y = 'Survival', color = 'Minor allele\ncount', size = 'Sex') +
		# scale_color_brewer(type = 'seq', palette = 'Blues') +
		scale_color_manual(values = c('#c6dbef', '#6baed6', '#08519c')) +
		scale_size_manual(values = c(1, 0.5)) +
		guides(color = guide_legend(override.aes = list(size = 0.75)))
} else {
	p = ggplot(aDt) +
		geom_step(aes(x = time, y = surv, color = factor(2 - snp)), size = 0.75) +
		labs(x = 'Age (y)', y = 'Survival', color = 'Minor allele\ncount') +
		scale_color_manual(values = c('#c6dbef', '#6baed6', '#08519c'))
}
print(p)

############################################################

# exome_compare_pval
# exome_examples_cumhaz
#
# resultSpread = result %>%
#   mutate(negLog10Pval = -log10(pval)) %>%
#   select(phecode, snp, method, negLog10Pval) %>%
#   spread(key = method, value = negLog10Pval) %>%
#   inner_join(select(phecodeData, phecode, phenotype), by = 'phecode') %>%
#   mutate(facetTitle = sprintf('%s (%s)', phenotype, phecode)) %>%
#   filter((cox >= 3) | (logistic >= 3))
#
# lineSz = 0.5
# ptShape = 16
# ptSz = 0.75
# ptAlpha = 0.5
# stripSz = 8
#
# p1 = ggplot(resultSpread) +
#   facet_wrap(~ facetTitle, nrow = 2, scales = 'free') +
#   geom_abline(slope = 1, intercept = 0, color = 'darkgray', size = lineSz) +
#   geom_point(aes(x = logistic, y = cox), shape = ptShape, size = ptSz, alpha = ptAlpha) +
#   theme(strip.text = element_text(size = stripSz))
#
# p2 = ggplot(resultSpread) +
#   facet_wrap(~ facetTitle, nrow = 2, scales = 'free') +
#   geom_hline(yintercept = 0, color = 'darkgray', size = lineSz) +
#   geom_point(aes(x = (cox + logistic)/2, y = cox - logistic), shape = ptShape, size = ptSz, alpha = ptAlpha) +
#   theme(strip.text = element_text(size = stripSz))
#
# p = plot_grid(ggdraw() + draw_label(analysisDir),
#               plot_grid(p1, p2, ncol = 1, align = 'v'),
#               ncol = 1, rel_heights = c(0.05, 1))
#
# ggsave(file.path(resultDir, analysisDir, 'compare_pval.pdf'), plot = p, width = 6, height = 8)
#
#
# survminer::ggcoxzph(cox.zph(coxFit))
#
# r = result %>%
#   filter(method == 'cox') %>%
#   group_by(phecode) %>%
#   arrange(pval) %>%
#   filter(row_number() == 1) %>%
#   ungroup() %>%
#   inner_join(phecodeData, by = 'phecode')
#
# registerDoParallel(cores = 1)
# denovoCalc = FALSE
#
# if (denovoCalc) {
#   d = foreach(ii = 1:nrow(r), .combine = rbind) %dopar% {
#     phenoDataNow = semi_join(phenoData, r[ii,], by = 'phecode')
#
#     if (r$whichSex[ii] == 'both') {
#       runCoxphNow = runCoxphSexBoth
#     } else {
#       runCoxphNow = runCoxphSexOne}
#
#     inputCox = makeCoxInput(phenoDataNow, gridData, genoData$genoFull, minEvents, buffer, r$whichSex[ii])
#     inputCox = addSnpToInput(inputCox, genoData$genoFull, r$snp[ii])
#     coxFit = runCoxphNow(inputCox)
#
#     dNow = makeSurvPlotDataframe(coxFit, inputCox, 'snp')
#     dNow = cbind(dNow, r[ii,])}
#
#   d = d %>%
#     mutate(value = 2 - value, # express as minor allele count
#            pheTitle = sprintf('%s (%s)', phenotype, phecode))
#   saveRDS(d, file.path(resultDir, analysisDir, 'examples_surv.rds'))
#
# } else {
#   d = readRDS(file.path(resultDir, analysisDir, 'examples_surv.rds'))}
#
#
# p = ggplot(d) +
#   facet_wrap(~ pheTitle + snp, nrow = 2, scales = 'free_y') +
#   geom_step(aes(x = time, y = 1 - surv, color = factor(value, levels = 2:0)), size = 1) +
#   scale_color_manual(values = c('#225ea8', '#41b6c4', '#a1dab4')) +
#   labs(title = analysisDir, x = 'Age (y)', y = 'Adjusted cumulative hazard', color = 'Minor allele\ncount') +
#   theme(strip.text = element_text(size = 8))
#
# ggsave(file.path(resultDir, analysisDir, 'examples_cumhaz.pdf'), plot = p, width = 7, height = 4.5)
#
#
# dRatio0 = d %>%
#   filter(value == 0) %>%
#   rename(surv0 = surv) %>%
#   select(phecode, snp, time, surv0)
#
# dRatio = d %>%
#   filter(value != 0) %>%
#   inner_join(dRatio0, by = c('phecode', 'snp', 'time')) %>%
#   mutate(ratioCumhaz = (1 - surv) / (1 - surv0))
#
# p = ggplot(dRatio) +
#   facet_wrap(~ pheTitle + snp, nrow = 2, scales = 'free_y') +
#   geom_step(aes(x = time, y = ratioCumhaz, color = factor(value, levels = 2:1)), size = 1) +
#   scale_color_manual(values = c('#225ea8', '#41b6c4', '#a1dab4')) +
#   labs(x = 'Age (y)', y = 'Ratio of adjusted cumulative hazards',
#        color = 'Minor allele\ncount', title = analysisDir) +
#   theme(strip.text = element_text(size = 8))
#
# ggsave(file.path(resultDir, analysisDir, 'examples_ratio_cumhaz.pdf'), plot = p, width = 7, height = 4.5)
