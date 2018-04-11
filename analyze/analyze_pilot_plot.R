library('readr')
library('dplyr')
library('cowplot')
library('doParallel')
library('qqman')

eb = element_blank()
theme_set(theme_light() +
			 	theme(axis.text=element_text(color='black'), strip.text=element_text(color='black'),
			 			panel.grid.minor=eb, legend.margin=margin(t=0, r=0, b=0, l=0, unit='cm')))

snpInfo = read_csv(file.path('processed', 'exome_map.csv'), col_types=cols()) %>%
	transmute(snpName = snp.name, chr = chromosome, pos = position)

phenotypeNames = c('alzheimers', 'atrial_fibrillation', 'gout',
						 'multiple_sclerosis', 'prostate_cancer', 'rheumatoid_arthritis')

resultFolder = 'results'

#######################################################

coxDf = foreach(phenotypeName=phenotypeNames[1], .combine=rbind) %do% {
	read_csv(file.path(resultFolder, sprintf('surv2/surv2_%s.csv', phenotypeName)), col_types=cols()) %>%
		mutate(phenotypeName = phenotypeName)}

a = inner_join(coxDf, snpInfo, by='snpName') %>%
	mutate(pval = ifelse(pval==0, 1e-20, pval))

for (phenotypeNow in phenotypeNames[1]) {
	pdf(file.path(resultFolder, sprintf('surv1_man_%s.pdf', phenotypeNow)), width=6, height=4)
	manhattan(a %>% filter(phenotypeName==phenotypeNow), p='pval', snp='snpName', chr='chr', bp='pos',
				 main = phenotypeNow, ylim=c(0, 21))
	dev.off()
}

for (phenotypeNow in phenotypeNames[1]) {
	pdf(file.path(resultFolder, sprintf('surv1_qq_%s.pdf', phenotypeNow)), width=4, height=4)
	qq(a %>% filter(phenotypeName==phenotypeNow) %>% .$pval, main = phenotypeNow)
	dev.off()
}

#######################################################

glmDf = foreach(phenotypeName=phenotypeNames, .combine=rbind) %do% {
	read_csv(file.path(resultFolder, sprintf('logistic2_reclen/logistic2_reclen_%s.csv', phenotypeName)),
				col_types=cols()) %>%
		mutate(phenotypeName = phenotypeName)}

a = inner_join(glmDf, snpInfo, by='snpName')

for (phenotypeNow in phenotypeNames) {
	pdf(file.path(resultFolder, sprintf('logistic2_reclen_man_%s.pdf', phenotypeNow)), width=6, height=4)
	manhattan(a %>% filter(phenotypeName==phenotypeNow), p='pval', snp='snpName', chr='chr', bp='pos',
				 main = phenotypeNow)
	dev.off()
}

for (phenotypeNow in phenotypeNames) {
	pdf(file.path(resultFolder, sprintf('logistic2_reclen_qq_%s.pdf', phenotypeNow)), width=4, height=4)
	qq(a %>% filter(phenotypeName==phenotypeNow) %>% .$pval, main = phenotypeNow)
	dev.off()
}

#######################################################

# clmDf = foreach(phenotypeName=phenotypeNames, .combine=rbind) %do% {
# 	read_csv(file.path(resultFolder, sprintf('cumord3_%s.csv', phenotypeName)), col_types=cols()) %>%
# 		mutate(phenotypeName = phenotypeName)}
#
# a = inner_join(clmDf, snpInfo, by='snpName') %>%
# 	filter(!is.na(pval))
#
# for (phenotypeNow in phenotypeNames) {
# 	pdf(file.path(resultFolder, sprintf('cumord3_man_%s.pdf', phenotypeNow)), width=6, height=4)
# 	manhattan(a %>% filter(phenotypeName==phenotypeNow), p='pval', snp='snpName', chr='chr', bp='pos',
# 				 main = phenotypeNow)
# 	dev.off()
# }
#
# for (phenotypeNow in phenotypeNames) {
# 	pdf(file.path(resultFolder, sprintf('cumord3_qq_%s.pdf', phenotypeNow)), width=4, height=4)
# 	qq(a %>% filter(phenotypeName==phenotypeNow) %>% .$pval, main = phenotypeNow)
# 	dev.off()
# }

#######################################################

df = inner_join(coxDf, glmDf, by=c('phenotypeName', 'snpName')) %>%
	transmute(phenotype = phenotypeName, snp = snpName, pvalCox = pval.x, pvalGlm = pval.y) %>%
	mutate(pvalCox = ifelse(pvalCox==0, 1e-20, pvalCox)) %>%
	filter(pvalCox <= 1e-3 | pvalGlm <= 1e-3)

p = ggplot(df) +
	facet_wrap(~ phenotype, nrow=2, scales='free') +
	geom_abline(slope=1, intercept=0, color='darkgray', size=0.5) +
	geom_point(aes(x=-log10(pvalGlm), y=-log10(pvalCox)), shape=16, size=0.75, alpha=0.5) +
	scale_x_continuous(breaks=seq(0, 30, 3)) + #limits=c(2, NA),
	scale_y_continuous(breaks=seq(0, 30, 3))

ggsave(file.path(resultFolder, 'compare_pval_surv1_logistic2_reclen.pdf'), plot=p, width=8, height=5.75)

#######################################################

# df = inner_join(coxDf, clmDf, by=c('phenotypeName', 'snpName')) %>%
# 	transmute(phenotype = phenotypeName, snp = snpName, pvalCox = pval.x, pvalClm = pval.y) %>%
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
# ggsave(file.path(resultFolder, 'compare_pval_surv1_cumord3.pdf'), plot=p, width=8, height=5.75)
