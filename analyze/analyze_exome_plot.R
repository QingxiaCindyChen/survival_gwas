source(file.path('analyze', 'analyze_setup.R'))

# analysisDir = 'exome_full'
analysisDir = 'exome_test'
# analysisDir = 'exome_pilot'

filePrefix = 'exome'

############################################################
# ensures that plotting uses exactly the same data used in the analysis,
# including gridData, phenoData, snpData, and all parameters and functions

load(file.path(resultDir, analysisDir, sprintf('%s_workspace.Rdata', filePrefix)))

# genoData = readRDS(file.path(procDir, 'exome_genotype_data.rds'))

############################################################

result = foreach(filename = phecodeDataKeep$filename, .combine = rbind) %do% {
	setDT(read_csv(file.path(resultDir, analysisDir, filename), col_types = '??????c'))
}

result = merge(result, snpData, by = 'snp', sort = FALSE)
result = merge(result, phecodeData[, .(phecode, phenotype)], by = 'phecode')

############################################################

r = unique(result[, .(phecode, phenotype, method)])

done = foreach(ii = 1:nrow(r)) %do% {
	resultNow = merge(result, r[ii,], by = colnames(r))

	plotTitle = sprintf('%s (%s), %s regression', r$phenotype[ii], r$phecode[ii], r$method[ii])
	pheFileText = gsub('.', 'p', r$phecode[ii], fixed = TRUE)

	filename = sprintf('%s_phe%s_%s_man.pdf', filePrefix, pheFileText, r$method[ii])
	pdf(file.path(resultDir, analysisDir, filename), width = 6, height = 4)
	manhattan(resultNow, p = 'pval', snp = 'snp', chr = 'chr', bp = 'pos', main = plotTitle)
	dev.off()

	filename = sprintf('%s_phe%s_%s_qq.pdf', filePrefix, pheFileText, r$method[ii])
	pdf(file.path(resultDir, analysisDir, filename), width = 6, height = 4)
	qq(resultNow$pval, main = plotTitle)
	dev.off()
}

############################################################

ptShp = 16
ptSz = 0.5
# ptShp = 21
# ptSz = 1
ptAlph = 0.2

lnSz = 0.5
lnCol = 'gray'

# hazard ratios are extremely similar to odds ratios
maxPval = 1e-5
resultSig = result[, if (mean(-log(pval)) >= -log(maxPval)) .SD, by = .(phecode, snp)]

# coefficients are based on major allele counts
resultSig[, c('logRatio', 'negLogPval') := .(log2(exp(-coef)), -log10(pval))]

resultSigEffect = dcast(resultSig, phecode + snp ~ method, value.var = 'logRatio')

p = ggplot(resultSigEffect) +
	geom_abline(slope = 1, intercept = 0, color = lnCol, size = lnSz) +
	geom_point(aes(x = logistic, y = cox), shape = 16, size = 0.5, alpha = ptAlph) +
	labs(title = 'Effect size', x = 'log2(hazard ratio)', y = 'log2(odds ratio)')

paramList = list(col = 'white', fill = 'darkgray', size = 0.25)
p1 = ggExtra::ggMarginal(p, type = 'histogram', binwidth = 0.15, boundary = 0,
								 xparams = paramList, yparams = paramList)

# p1 = ggplot(resultSigEffect) +
# 	geom_hline(yintercept = 0, color = lnCol, size = lnSz) +
# 	geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
# 	geom_smooth(aes(x = (logistic + cox) / 2, y = cox - logistic), size = 0.5, method = 'loess') +
# 	labs(title = 'Effect size', x = 'log2(hazard ratio * odds ratio) / 2',
# 		  y = 'log2(hazard ratio / odds ratio)') +
# 	scale_x_continuous(limits = c(-4, 4)) + scale_y_continuous(limits = c(-0.1, 0.1))

cor(resultSigEffect$logistic, resultSigEffect$cox, use = 'na.or.complete')

resultSigPval = dcast(resultSig, phecode + snp ~ method, value.var = 'negLogPval')
p2 = ggplot(resultSigPval) +
	geom_hline(yintercept = 0, color = lnCol, size = lnSz) +
	geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
	geom_smooth(aes(x = (logistic + cox) / 2, y = cox - logistic), size = 0.5, method = 'loess') +
	labs(title = '-log10(p)')

# standard errors are slightly lower for cox
resultSigSe = dcast(resultSig, phecode + snp ~ method, value.var = 'se')
p3 = ggplot(resultSigSe) +
	geom_histogram(aes(x = cox - logistic), binwidth = 0.001, boundary = 0,
						fill = 'darkgray', color = 'white') +
	labs(title = 'Standard error') +
	scale_x_continuous(limits = c(-0.015, 0.005))

# genomic inflation factors are similar or slightly higher
resultTmp = result[, .(lambdaMed = median(z^2) / qchisq(0.5, 1)), by = .(phecode, method)]
resultLambda = dcast(resultTmp, phecode ~ method, value.var = 'lambdaMed')

p4 = ggplot(resultLambda) +
	geom_hline(yintercept = 0, color = lnCol, size = lnSz) +
	geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic), shape = ptShp, size = ptSz, alpha = ptAlph) +
	labs(title = 'Lambda median') +
	scale_x_continuous(limits = c(0.6, NA))

p = plot_grid(p1, p2, p3, p4, align = 'hv', axis = 'tb', nrow = 2)
# p = ggdraw() +
# 	draw_plot(p1, x = 0, y = 0.5, width = 0.5, height = 0.5) +
# 	draw_plot(p2, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
# 	draw_plot(p3, x = 0, y = 0, width = 0.5, height = 0.5) +
# 	draw_plot(p4, x = 0.5, y = 0, width = 0.5, height = 0.5)
ggsave(file.path(resultDir, 'exome_test.pdf'), plot = p, width = 6, height = 5.75)



result1 = dcast(result, phecode + snp ~ method, value.var = 'pval')
result2 = result1[, .(r = cor(-log10(logistic), -log10(cox))), by = phecode]
p = ggplot(result2) +
	# geom_step(aes(x = r, y = 1 - ..y..), stat = 'ecdf')
	stat_ecdf(aes(x = r), pad = FALSE)
p



# a = merge(resultSigPval[((cox + logistic)/2 < 10) & (cox - logistic > 1),], phecodeData, by = 'phecode')
# a = merge(resultSigPval[logistic - cox > 1,], phecodeData, by = 'phecode')
a = merge(resultSigPval[cox - logistic > 1,], phecodeData, by = 'phecode')
a = merge(a, snpData, by = 'snp')
a = a[order(cox - logistic, decreasing = TRUE),]


