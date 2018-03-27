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

#######################################################

coxDf = foreach(phenotypeName=phenotypeNames, .combine=rbind) %do% {
	read_csv(file.path('results', sprintf('count3_%s.csv', phenotypeName)), col_types=cols()) %>%
		mutate(phenotypeName = phenotypeName)}

a = inner_join(coxDf, snpInfo, by='snpName') %>%
	mutate(pval = ifelse(pval==0, 1e-20, pval))

for (phenotypeNow in phenotypeNames) {
	pdf(file.path('results', sprintf('count3_man_%s.pdf', phenotypeNow)), width=6, height=4)
	manhattan(a %>% filter(phenotypeName==phenotypeNow), p='pval', snp='snpName', chr='chr', bp='pos',
				 main = phenotypeNow, ylim=c(0, 21))
	dev.off()
}

for (phenotypeNow in phenotypeNames) {
	pdf(file.path('results', sprintf('count3_qq_%s.pdf', phenotypeNow)), width=4, height=4)
	qq(a %>% filter(phenotypeName==phenotypeNow) %>% .$pval, main = phenotypeNow)
	dev.off()
}

#######################################################

glmDf = foreach(phenotypeName=phenotypeNames, .combine=rbind) %do% {
	read_csv(file.path('results', sprintf('naive2_%s.csv', phenotypeName)), col_types=cols()) %>%
		mutate(phenotypeName = phenotypeName)}

a = inner_join(glmDf, snpInfo, by='snpName')

# for (phenotypeNow in phenotypeNames) {
# 	pdf(file.path('results', sprintf('naive2_man_%s.pdf', phenotypeNow)), width=6, height=4)
# 	manhattan(a %>% filter(phenotypeName==phenotypeNow), p='pval', snp='snpName', chr='chr', bp='pos',
# 				 main = phenotypeNow)
# 	dev.off()
# }
#
# for (phenotypeNow in phenotypeNames) {
# 	pdf(file.path('results', sprintf('naive2_qq_%s.pdf', phenotypeNow)), width=4, height=4)
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

ggsave(file.path('results', 'compare_pval_count3_naive2.pdf'), plot=p, width=8, height=5.75)
