library('readr')
library('dplyr')
library('cowplot')
library('doParallel')
library('qqman')

eb = element_blank()
theme_set(theme_light() +
			 	theme(axis.text=element_text(color='black'), strip.text=element_text(color='black'),
			 			panel.grid.minor=eb, legend.margin=margin(t=0, r=0, b=0, l=0, unit='cm')))

procDir = 'processed'
resultDir = 'results'

snpInfo = read_csv(file.path(procDir, 'exome_map.csv.gz'), col_types=cols()) %>%
	transmute(snpName = snp.name, chr = chromosome, pos = position)

phenotypes = c('alzheimers', 'atrial_fibrillation', 'gout',
					'multiple_sclerosis', 'prostate_cancer', 'rheumatoid_arthritis')

#######################################################

idx = 1:6

coxDf = foreach(phenotypeNow = phenotypes[idx], .combine = rbind) %do% {
	read_csv(file.path(resultDir, sprintf('surv2/surv2_%s.csv.gz', phenotypeNow)), col_types=cols()) %>%
		mutate(phenotype = phenotypeNow)}

a = inner_join(coxDf, snpInfo, by='snpName') %>%
	mutate(pval = ifelse(pval==0, 1e-20, pval))

for (phenotypeNow in phenotypes[idx]) {
	pdf(file.path(resultDir, sprintf('surv2_man_%s.pdf', phenotypeNow)), width=6, height=4)
	manhattan(a %>% filter(phenotype==phenotypeNow), p='pval', snp='snpName', chr='chr', bp='pos',
				 main = phenotypeNow, ylim=c(0, 21))
	dev.off()}

for (phenotypeNow in phenotypes[idx]) {
	pdf(file.path(resultDir, sprintf('surv2_qq_%s.pdf', phenotypeNow)), width=4, height=4)
	qq(a %>% filter(phenotype==phenotypeNow) %>% .$pval, main = phenotypeNow)
	dev.off()}

#######################################################

idx = 1:6

glmDf = foreach(phenotypeNow = phenotypes[idx], .combine = rbind) %do% {
	read_csv(file.path(resultDir, sprintf('logistic2/logistic2_%s.csv.gz', phenotypeNow)), col_types=cols()) %>%
		mutate(phenotype = phenotypeNow)}

a = inner_join(glmDf, snpInfo, by='snpName')

for (phenotypeNow in phenotypes[idx]) {
	pdf(file.path(resultDir, sprintf('logistic2_man_%s.pdf', phenotypeNow)), width=6, height=4)
	manhattan(a %>% filter(phenotype==phenotypeNow), p='pval', snp='snpName', chr='chr', bp='pos',
				 main = phenotypeNow)
	dev.off()}

for (phenotypeNow in phenotypes[idx]) {
	pdf(file.path(resultDir, sprintf('logistic2_qq_%s.pdf', phenotypeNow)), width=4, height=4)
	qq(a %>% filter(phenotype==phenotypeNow) %>% .$pval, main = phenotypeNow)
	dev.off()}

#######################################################

# clmDf = foreach(phenotypeNow = phenotypes, .combine = rbind) %do% {
# 	read_csv(file.path(resultDir, sprintf('cumord3_%s.csv.gz', phenotypeNow)), col_types=cols()) %>%
# 		mutate(phenotype = phenotypeNow)}
#
# a = inner_join(clmDf, snpInfo, by='snpName') %>%
# 	filter(!is.na(pval))
#
# for (phenotypeNow in phenotypes) {
# 	pdf(file.path(resultDir, sprintf('cumord3_man_%s.pdf', phenotypeNow)), width=6, height=4)
# 	manhattan(a %>% filter(phenotype==phenotypeNow), p='pval', snp='snpName', chr='chr', bp='pos',
# 				 main = phenotypeNow)
# 	dev.off()}
#
# for (phenotypeNow in phenotypes) {
# 	pdf(file.path(resultDir, sprintf('cumord3_qq_%s.pdf', phenotypeNow)), width=4, height=4)
# 	qq(a %>% filter(phenotype==phenotypeNow) %>% .$pval, main = phenotypeNow)
# 	dev.off()}

#######################################################

df = inner_join(coxDf, glmDf, by=c('phenotype', 'snpName')) %>%
	transmute(phenotype = phenotype, snp = snpName, pvalCox = pval.x, pvalGlm = pval.y) %>%
	mutate(pvalCox = ifelse(pvalCox==0, 1e-20, pvalCox)) %>%
	filter(pvalCox <= 1e-3 | pvalGlm <= 1e-3)

p = ggplot(df) +
	facet_wrap(~ phenotype, nrow=2, scales='free') +
	geom_abline(slope=1, intercept=0, color='darkgray', size=0.5) +
	geom_point(aes(x=-log10(pvalGlm), y=-log10(pvalCox)), shape=16, size=0.75, alpha=0.5)

ggsave(file.path(resultDir, 'compare_pval_surv2_logistic2.pdf'), plot=p, width=8, height=5.75)

#######################################################

# df = inner_join(coxDf, clmDf, by=c('phenotype', 'snpName')) %>%
# 	transmute(phenotype = phenotype, snp = snpName, pvalCox = pval.x, pvalClm = pval.y) %>%
# 	mutate(pvalCox = ifelse(pvalCox==0, 1e-20, pvalCox)) %>%
# 	filter(pvalCox <= 1e-3 | pvalClm <= 1e-3)
#
# p = ggplot(df) +
# 	facet_wrap(~ phenotype, nrow=2, scales='free') +
# 	geom_abline(slope=1, intercept=0, color='darkgray', size=0.5) +
# 	geom_point(aes(x=-log10(pvalClm), y=-log10(pvalCox)), shape=16, size=0.75, alpha=0.5) +
# 	scale_x_continuous(breaks=seq(0, 30, 3)) + #limits=c(2, NA),
# 	scale_y_continuous(breaks=seq(0, 30, 3))
#
# ggsave(file.path(resultDir, 'compare_pval_surv2_cumord3.pdf'), plot=p, width=8, height=5.75)
