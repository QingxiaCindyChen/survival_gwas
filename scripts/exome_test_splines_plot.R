source(file.path('analyze', 'analyze_setup.R'))

load(file.path('results/exome_test_splines.Rdata'))

############################################################

ptShp = 16
ptSz = 0.5
ptAlph = 0.2
lnSz = 0.5
lnCol = 'gray'

dPval = dcast(d1[pval <= 1e-1], phecode + snp ~ ageDf, value.var = 'pval')
colnames(dPval)[3:ncol(dPval)] = paste0('ageDf', colnames(dPval)[3:ncol(dPval)])

pPval = ggplot(dPval) +
  geom_point(aes(x = -log10(ageDf3), y = -log10(ageDf4)), shape = 16, size = 0.5, alpha = ptAlph) +
  labs(title = 'P-value for genotype', x = '3 d.f. for last_age', y = '4 d.f. for last_age')

############################################################

dBeta = dcast(d1[pval <= 1e-4], phecode + snp ~ ageDf, value.var = 'beta')
colnames(dBeta)[3:ncol(dBeta)] = paste0('ageDf', colnames(dBeta)[3:ncol(dBeta)])

pBetaTmp = ggplot(dBeta) +
  geom_abline(slope = 1, intercept = 0, color = lnCol, size = lnSz) +
  geom_point(aes(x = -ageDf3, y = -ageDf4), shape = 16, size = 0.5, alpha = ptAlph) +
  labs(title = 'Beta for genotype', x = '3 d.f. for last_age', y = '4 d.f. for last_age')

paramList = list(fill = 'gray', col = 'black', size = 0.25)
pBeta = ggExtra::ggMarginal(pBetaTmp, type = 'histogram', binwidth = 0.1, boundary = 0,
                            xparams = paramList, yparams = paramList)

############################################################

dVifGeno = d2[variable == 'snp']
dVifGeno[, ageDfText := sprintf('%d d.f. for last_age', ageDf)]
dVifGeno[, ageDfSz := factor(ifelse(ageDf == 3, 0.6, 0.4))]

# pGeno = ggplot(dVifGeno) +
#   facet_wrap(~ ageDfText, ncol = 1) +
#   geom_histogram(aes(x = vif), boundary = 1, binwidth = 0.0005,
#                  fill = 'gray', color = 'black', size = 0.25) +
#   labs(x = 'VIF for genotype', y = 'Number of tested associations') +
#   scale_x_continuous(limits = c(1, 1.01))

pGeno = ggplot(dVifGeno) +
  stat_ecdf(aes(x = vif, color = factor(ageDf), size = ageDfSz), pad = FALSE) +
  scale_color_brewer(type = 'qual', palette = 'Paired') +
  scale_size_discrete(guide = FALSE) +
  labs(x = 'VIF for genotype', y = 'Empirical CDF', color = 'd.f. for\nlast_age')
# pGeno

############################################################

dVifAge = d2[startsWith(variable, 'last_age')]
dVifAge[, ageDfText := sprintf('%d d.f. for last_age', ageDf)]

# pAge = ggplot(dVifAge) +
#   facet_wrap(~ ageDfText, ncol = 1) +
#   geom_histogram(aes(x = log10(vif)), boundary = 1, binwidth = 0.25,
#                  fill = 'gray', color = 'black', size = 0.25) +
#   labs(x = 'log10 VIF for last_age*', y = 'Number of tested associations') +
#   scale_x_continuous(limits = c(0, 5))

pAge = ggplot(dVifAge) +
  stat_ecdf(aes(x = log10(vif), color = factor(ageDf)), pad = FALSE) +
  scale_color_brewer(type = 'qual', palette = 'Paired') +
  labs(x = 'log10 VIF for last_age*', y = 'Empirical CDF', color = 'd.f. for\nlast_age')


maxVif = 2e6
dTmp1 = dVifAge[, .(vif = max(vif)), by = .(phecode, snpName, ageDf)]
dTmp2 = dTmp1[(vif >= maxVif), .N, by = .(ageDf, phecode)]

############################################################

p = plot_grid(pPval, pBeta, pGeno, pAge, align = 'hv', axis = 'tb', nrow = 2)
ggsave(file.path(resultDir, 'exome_test_splines.pdf'), plot = p, width = 7, height = 7)
