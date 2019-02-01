source(file.path('scripts', 'setup_regression.R'))

plotDir = 'results/mega/20181030_104520/plots'
pcRaw = readRDS('processed/mega/pc_data_raw.rds')

nPc = 10
pcSummary = data.table(pc = 1:nPc, eigenval = pcRaw$eigenval[1:nPc])
pcData = data.table(pcRaw$eigenvect)
colnames(pcData) = paste0('PC', 1:ncol(pcData))

pScree = ggplot(pcSummary) +
  geom_point(aes(x = pc, y = eigenval), size = 2) +
  labs(x = 'Principal component', y = 'Eigenvalue') +
  scale_x_continuous(breaks = seq(1, 10, 2))

sz = 0.4
sp = 16
alph = 0.3

p12 = ggplot(pcData) +
  geom_point(aes(x = PC1, y = PC2), size = sz, shape = sp, alpha = alph)

p = plot_grid(pScree, p12, nrow = 1, labels = 'AUTO', rel_widths = )
ggsave(file.path(plotDir, 'pc_data.pdf'), plot = p, width = 6, height = 2.8)
