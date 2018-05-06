source(file.path('analyze', 'analyze_setup.R'))

procDir = 'processed'
resultDir = 'results'
analysisDir = 'exome_minEvents2_pilot'

minEvents = 2
buffer = 1 # years, used for cox regression

phecodeData = read_csv(file.path(procDir, 'phecode_definitions1.2.csv'), col_types = 'ccccii') %>%
	rename(phecode = jd_code,
			 phenotype = jd_string,
			 controlExcludeRange = jd_control_exclude_range,
			 whichSex = sex) %>%
	mutate(whichSex = tolower(ifelse(is.na(whichSex), 'both', whichSex)))

minRecLen = 0 # years
gridData = read_csv(file.path(procDir, 'grid_data.csv.gz'), col_types = 'ccDDD')
gridData = gridData %>%
	mutate(first_age = time_length(first_entry_date - dob, 'years'),
			 last_age = time_length(last_entry_date - dob, 'years'),
			 rec_len = last_age - first_age) %>%
	filter(first_age >= 0,
			 rec_len >= minRecLen)

minGrids = 50
phenoData = read_csv(file.path(procDir, 'phenotype_data.csv.gz'), col_types = 'cDc') %>%
	group_by(phecode) %>%
	filter(n_distinct(grid) >= minGrids) %>%
	ungroup()

genoData = readRDS(file.path(procDir, 'genotype_data_exome.rds'))
snpData = read_csv(file.path(procDir, 'exome_map.csv.gz'), col_types = cols()) %>%
	transmute(snp = snp.name, chr = chromosome, pos = position)

#######################################################

files = list.files(file.path(resultDir, analysisDir), pattern = '\\.csv\\.gz$')
resultTmp = foreach(fileNow = files, .combine = rbind) %do% {
	read_csv(file.path(resultDir, analysisDir, fileNow), col_types = 'cc????c')}

result = resultTmp %>%
	mutate(pval = ifelse(pval==0, 2 * pnorm(-abs(z)), pval)) %>%
	inner_join(snpData, by = 'snp')

#######################################################

r = result %>%
	distinct(method, phecode) %>%
	inner_join(select(phecodeData, phecode, phenotype), by = 'phecode')

done = foreach(ii = 1:nrow(r)) %do% {
	resultNow = semi_join(result, r[ii,])
	plotTitle = sprintf('%s (%s), %s regression', r$phenotype[ii], r$phecode[ii], r$method[ii])
	pheFileText = gsub('.', 'p', r$phecode[ii], fixed = TRUE)

	filename = sprintf('phe%s_%s_man.pdf', pheFileText, r$method[ii])
	pdf(file.path(resultDir, analysisDir, filename), width = 6, height = 4)
	manhattan(resultNow, p = 'pval', snp = 'snp', chr = 'chr', bp = 'pos', main = plotTitle)
	dev.off()

	filename = sprintf('phe%s_%s_qq.pdf', pheFileText, r$method[ii])
	pdf(file.path(resultDir, analysisDir, filename), width = 6, height = 4)
	qq(resultNow$pval, main = plotTitle)
	dev.off()
}

#######################################################

resultSpread = result %>%
	mutate(negLog10Pval = -log10(pval)) %>%
	select(phecode, snp, method, negLog10Pval) %>%
	spread(key = method, value = negLog10Pval) %>%
	inner_join(select(phecodeData, phecode, phenotype), by = 'phecode') %>%
	mutate(facetTitle = sprintf('%s (%s)', phenotype, phecode)) %>%
	filter((cox >= 3) | (logistic >= 3))

lineSz = 0.5
ptShape = 16
ptSz = 0.75
ptAlpha = 0.5
stripSz = 8

p1 = ggplot(resultSpread) +
	facet_wrap(~ facetTitle, nrow = 2, scales = 'free') +
	geom_abline(slope = 1, intercept = 0, color = 'darkgray', size = lineSz) +
	geom_point(aes(x = logistic, y = cox), shape = ptShape, size = ptSz, alpha = ptAlpha) +
	theme(strip.text = element_text(size = stripSz))

p2 = ggplot(resultSpread) +
	facet_wrap(~ facetTitle, nrow = 2, scales = 'free') +
	geom_hline(yintercept = 0, color = 'darkgray', size = lineSz) +
	geom_point(aes(x = (cox + logistic)/2, y = cox - logistic), shape = ptShape, size = ptSz, alpha = ptAlpha) +
	theme(strip.text = element_text(size = stripSz))

p = plot_grid(ggdraw() + draw_label(analysisDir),
				  plot_grid(p1, p2, ncol = 1, align = 'v'),
				  ncol = 1, rel_heights = c(0.05, 1))

ggsave(file.path(resultDir, analysisDir, 'compare_pval.pdf'), plot = p, width = 6, height = 8)

########################################

r = result %>%
	filter(method == 'cox') %>%
	group_by(phecode) %>%
	arrange(pval) %>%
	filter(row_number() == 1) %>%
	ungroup() %>%
	inner_join(phecodeData, by = 'phecode')

registerDoParallel(cores = 1)

d = foreach(ii = 1:nrow(r), .combine = rbind) %dopar% {
	phenoDataNow = semi_join(phenoData, r[ii,], by = 'phecode')

	coxInput = makeCoxInput(phenoDataNow, gridData, genoData$genoFull, minEvents, buffer, r$whichSex[ii])
	coxInput$snp = as(genoData$genoFull$genotypes[coxInput$grid, r$snp[ii]], 'numeric')[,1]
	coxInput = filter(coxInput, !is.na(snp))

	coxFit = runCoxph(coxInput, genoData$genoFull, r$snp[ii], r$whichSex[ii])
	dNow = makeSurvPlotDataframe(coxFit, coxInput, 'snp')
	dNow = cbind(dNow, r[ii,])
}

d = mutate(d, pheTitle = sprintf('%s (%s)', phenotype, phecode))

p = ggplot(d) +
	facet_wrap(~ pheTitle + snp, nrow = 2, scales = 'free_y') +
	geom_step(aes(x = time, y = 1 - surv, color = factor(2 - value, levels = 2:0)), size = 1) +
	scale_color_manual(values = c('#225ea8', '#41b6c4', '#a1dab4')) +
	labs(title = analysisDir, x = 'Age (y)', y = 'Adjusted cumulative hazard', color = 'Minor allele\ncount') +
	theme(strip.text = element_text(size = 8))

ggsave(file.path(resultDir, analysisDir, 'cumhaz_examples.pdf'), plot = p, width = 7, height = 4.5)

# survminer::ggcoxzph(cox.zph(coxFit))
