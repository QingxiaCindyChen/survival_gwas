source(file.path('analyze', 'analyze_setup.R'))

analysisDir = 'exome_full'
# analysisDir = 'exome_test_nPC2'
filePrefix = 'exome'

############################################################
# ensures that plotting uses exactly the same data used in the analysis,
# including gridData, phenoData, snpData, and all parameters and functions

load(file.path(resultDir, analysisDir, sprintf('%s_workspace.Rdata', filePrefix)))

gwasFilenames = list.files(file.path(resultDir, analysisDir), pattern = '\\.csv\\.gz$')

# genoData = readRDS(file.path(procDir, 'exome_genotype_data.rds'))

############################################################

result = foreach(filename = gwasFilenames, .combine = rbind) %do% {
	read_csv(file.path(resultDir, analysisDir, filename), col_types = '??????c')}

result = merge(setDT(result), snpData, by = 'snp', sort = FALSE)
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

# hazard ratios are extremely similar to odds ratios
maxPval = 1e-5
resultSig = result[, if (mean(-log(pval)) >= -log(maxPval)) .SD, by = .(phecode, snp)]

# coefficients are based on major allele counts
resultSig[, c('logRatio', 'negLogPval') := .(log2(exp(-coef)), -log10(pval))]

resultSigEffect = dcast(resultSig, phecode + snp ~ method, value.var = 'logRatio')
p1 = ggplot(resultSigEffect) +
	geom_abline(slope = 1, intercept = 0, color = 'lightgray', size = 0.5) +
	geom_point(aes(x = logistic, y = cox), shape = 16, size = 0.5, alpha = 0.5) +
	labs(x = 'Log2 odds ratio (logistic)', y = 'Log2 hazard ratio (cox)', title = 'Effect size')

cor(resultSigEffect$logistic, resultSigEffect$cox, use = 'na.or.complete')

resultSigPval = dcast(resultSig, phecode + snp ~ method, value.var = 'negLogPval')
p2 = ggplot(resultSigPval) +
	geom_hline(yintercept = 0, color = 'darkgray', size = 0.5) +
	geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic), shape = 16, size = 0.75, alpha = 0.5) +
	geom_smooth(aes(x = (logistic + cox) / 2, y = cox - logistic), size = 0.5, method = 'loess') +
	labs(title = '-log10(p)')

# standard errors are slightly lower for cox
resultSigSe = dcast(resultSig, phecode + snp ~ method, value.var = 'se')
p3 = ggplot(resultSigSe) +
	geom_histogram(aes(x = cox - logistic), binwidth = 0.001, boundary = 0,
						fill = 'darkgray', color = 'white') +
	scale_x_continuous(limits = c(-0.02, 0.01)) +
	labs(title = 'Standard error')

# genomic inflation factors are similar or slightly higher
resultTmp = result[, .(lambdaMed = median(z^2) / qchisq(0.5, 1)), by = .(phecode, method)]
resultLambda = dcast(resultTmp, phecode ~ method, value.var = 'lambdaMed')

p4 = ggplot(resultLambda) +
	geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
	geom_point(aes(x = (logistic + cox) / 2, y = cox - logistic), shape = 16, size = 0.75, alpha = 0.5) +
	labs(title = 'Lambda median')

p = plot_grid(p1, p2, p3, p4, align = 'hv', nrow = 2)
ggsave(file.path(resultDir, 'exome_full_prelim.pdf'), plot = p, width = 6, height = 5.75)
