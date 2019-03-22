library('cowplot')
library('data.table')

theme_set(theme_classic() +
            theme(axis.text = element_text(color = 'black'),
                  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))

resultDir = 'results/sim'
simData = fread(file.path(resultDir, 'TruncOutput.txt'))

mean(simData$TCoxTPR)
mean(simData$LogitTPR)

ttestResult = t.test(simData$TCoxTPR, simData$LogitTPR, paired = TRUE)

1 - mean(simData$TCoxTNR)
1 - mean(simData$LogitTNR)

p = ggplot(simData) +
  geom_abline(slope = 1, intercept = 0, size = 0.5, color = 'gray') +
  geom_count(aes(x = LogitTPR, y = TCoxTPR), shape = 16, alpha = 0.8) +
  scale_size_area(max_size = 2) +
  labs(x = 'TPR of logistic regression', y = 'TPR of Cox regression',
       size = 'Number of\nsimulations')

ggsave(file.path(resultDir, 'sim_tpr.pdf'),
       plot = p, width = 3.5, height = 2.5)
