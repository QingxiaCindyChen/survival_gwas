# exome_compare_pval
# exome_examples_cumhaz

############################################################

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
