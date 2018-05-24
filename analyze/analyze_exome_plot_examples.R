source(file.path('analyze', 'analyze_setup.R'))

analysisDir = 'exome_full'
# analysisDir = 'exome_test_nPC2'
# analysisDir = 'exome_pilot'

filePrefix = 'exome'

############################################################
# ensures that plotting uses exactly the same data used in the analysis,
# including gridData, phenoData, snpData, and all parameters and functions

load(file.path(resultDir, analysisDir, sprintf('%s_workspace.Rdata', filePrefix)))

gwasFilenames = list.files(file.path(resultDir, analysisDir), pattern = '\\.csv\\.gz$')
gwasFilenames = c(gwasFilenames[grepl('phe250p[1-2]\\.csv\\.gz', gwasFilenames)],
						gwasFilenames[grepl('phe290p11\\.csv\\.gz', gwasFilenames)],
						gwasFilenames[grepl('phe335\\.csv\\.gz', gwasFilenames)])

genoData = readRDS(file.path(procDir, 'exome_genotype_data.rds'))

############################################################

result = foreach(filename = gwasFilenames, .combine = rbind) %do% {
	setDT(read_csv(file.path(resultDir, analysisDir, filename), col_types = '??????c'))}

result = merge(result, snpData, by = 'snp', sort = FALSE)
result = merge(result, phecodeData[, .(phecode, phenotype)], by = 'phecode')

############################################################

minGrids = 50
phenoData = setDT(read_csv(file.path(procDir, sprintf('%s_phenotype_data.csv.gz', filePrefix)),
									col_types = 'ccD'))
phenoData = phenoData[, if (length(unique(grid)) >= minGrids) .SD, by = phecode]
phenoData = phenoData[order(grid, phecode, entry_date), .(grid, phecode, entry_date)]

phenoData = merge(phenoData, gridData, by = 'grid')
phenoData[, age := time_length(entry_date - dob, 'years')]
phenoData = phenoData[age <= maxAgeAtEvent,]




phecodeNow = '250.1'
snpNow = 'rs9275495'

# phecodeNow = '250.2'
# snpNow = 'rs7903146'

# phecodeNow = '290.11'
# snpNow = 'rs769449'

whichSexNow = phecodeData[phecode == phecodeNow, whichSex]
phenoDataNow = phenoData[phecode == phecodeNow, .(grid, age)]

inputBase = makeInput(phenoDataNow, gridData, genoData$genoFull, whichSexNow, minEvents, buffer)
inputNow = addSnpToInput(inputBase, genoData$genoFull, snpNow)

coxStr = getCoxStr(whichSexNow, 0)
a = survfit(as.formula(coxStr), data = inputNow)

p = survminer::ggsurvplot(a)
pdf('example_phe250p1_rs9275495.pdf', width = 5, height = 4)
p
dev.off()


############################################################

# exome_compare_pval
# exome_examples_cumhaz

# resultSpread = result %>%
# 	mutate(negLog10Pval = -log10(pval)) %>%
# 	select(phecode, snp, method, negLog10Pval) %>%
# 	spread(key = method, value = negLog10Pval) %>%
# 	inner_join(select(phecodeData, phecode, phenotype), by = 'phecode') %>%
# 	mutate(facetTitle = sprintf('%s (%s)', phenotype, phecode)) %>%
# 	filter((cox >= 3) | (logistic >= 3))
#
# lineSz = 0.5
# ptShape = 16
# ptSz = 0.75
# ptAlpha = 0.5
# stripSz = 8
#
# p1 = ggplot(resultSpread) +
# 	facet_wrap(~ facetTitle, nrow = 2, scales = 'free') +
# 	geom_abline(slope = 1, intercept = 0, color = 'darkgray', size = lineSz) +
# 	geom_point(aes(x = logistic, y = cox), shape = ptShape, size = ptSz, alpha = ptAlpha) +
# 	theme(strip.text = element_text(size = stripSz))
#
# p2 = ggplot(resultSpread) +
# 	facet_wrap(~ facetTitle, nrow = 2, scales = 'free') +
# 	geom_hline(yintercept = 0, color = 'darkgray', size = lineSz) +
# 	geom_point(aes(x = (cox + logistic)/2, y = cox - logistic), shape = ptShape, size = ptSz, alpha = ptAlpha) +
# 	theme(strip.text = element_text(size = stripSz))
#
# p = plot_grid(ggdraw() + draw_label(analysisDir),
# 				  plot_grid(p1, p2, ncol = 1, align = 'v'),
# 				  ncol = 1, rel_heights = c(0.05, 1))
#
# ggsave(file.path(resultDir, analysisDir, 'compare_pval.pdf'), plot = p, width = 6, height = 8)

############################################################

# survminer::ggcoxzph(cox.zph(coxFit))

# r = result %>%
# 	filter(method == 'cox') %>%
# 	group_by(phecode) %>%
# 	arrange(pval) %>%
# 	filter(row_number() == 1) %>%
# 	ungroup() %>%
# 	inner_join(phecodeData, by = 'phecode')
#
# registerDoParallel(cores = 1)
# denovoCalc = FALSE
#
# if (denovoCalc) {
# 	d = foreach(ii = 1:nrow(r), .combine = rbind) %dopar% {
# 		phenoDataNow = semi_join(phenoData, r[ii,], by = 'phecode')
#
# 		if (r$whichSex[ii] == 'both') {
# 			runCoxphNow = runCoxphSexBoth
# 		} else {
# 			runCoxphNow = runCoxphSexOne}
#
# 		inputCox = makeCoxInput(phenoDataNow, gridData, genoData$genoFull, minEvents, buffer, r$whichSex[ii])
# 		inputCox = addSnpToInput(inputCox, genoData$genoFull, r$snp[ii])
# 		coxFit = runCoxphNow(inputCox)
#
# 		dNow = makeSurvPlotDataframe(coxFit, inputCox, 'snp')
# 		dNow = cbind(dNow, r[ii,])}
#
# 	d = d %>%
# 		mutate(value = 2 - value, # express as minor allele count
# 				 pheTitle = sprintf('%s (%s)', phenotype, phecode))
# 	saveRDS(d, file.path(resultDir, analysisDir, 'examples_surv.rds'))
#
# } else {
# 	d = readRDS(file.path(resultDir, analysisDir, 'examples_surv.rds'))}

############################################################

# p = ggplot(d) +
# 	facet_wrap(~ pheTitle + snp, nrow = 2, scales = 'free_y') +
# 	geom_step(aes(x = time, y = 1 - surv, color = factor(value, levels = 2:0)), size = 1) +
# 	scale_color_manual(values = c('#225ea8', '#41b6c4', '#a1dab4')) +
# 	labs(title = analysisDir, x = 'Age (y)', y = 'Adjusted cumulative hazard', color = 'Minor allele\ncount') +
# 	theme(strip.text = element_text(size = 8))
#
# ggsave(file.path(resultDir, analysisDir, 'examples_cumhaz.pdf'), plot = p, width = 7, height = 4.5)

############################################################

# dRatio0 = d %>%
# 	filter(value == 0) %>%
# 	rename(surv0 = surv) %>%
# 	select(phecode, snp, time, surv0)
#
# dRatio = d %>%
# 	filter(value != 0) %>%
# 	inner_join(dRatio0, by = c('phecode', 'snp', 'time')) %>%
# 	mutate(ratioCumhaz = (1 - surv) / (1 - surv0))
#
# p = ggplot(dRatio) +
# 	facet_wrap(~ pheTitle + snp, nrow = 2, scales = 'free_y') +
# 	geom_step(aes(x = time, y = ratioCumhaz, color = factor(value, levels = 2:1)), size = 1) +
# 	scale_color_manual(values = c('#225ea8', '#41b6c4', '#a1dab4')) +
# 	labs(x = 'Age (y)', y = 'Ratio of adjusted cumulative hazards',
# 		  color = 'Minor allele\ncount', title = analysisDir) +
# 	theme(strip.text = element_text(size = 8))
#
# ggsave(file.path(resultDir, analysisDir, 'examples_ratio_cumhaz.pdf'), plot = p, width = 7, height = 4.5)
