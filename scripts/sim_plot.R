library('cowplot')
library('data.table')

theme_set(theme_light() +
            theme(axis.text = element_text(color = 'black'),
                  strip.text = element_text(color = 'black'),
                  panel.grid.minor = element_blank(),
                  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))

resultDir = 'results/sim'
simData = fread(file.path(resultDir, 'TruncOutput.txt'))

mean(simData$TCoxTPR)
mean(simData$LogitTPR)

mean(simData$TCoxTNR)
mean(simData$LogitTNR)

p = ggplot(simData) +
  geom_abline(slope = 1, intercept = 0, size = 0.5, color = 'gray') +
  geom_count(aes(x = LogitTPR, y = TCoxTPR), shape = 16, alpha = 0.8) +
  scale_size_area(max_size = 3) +
  labs(x = 'TPR of logistic regression', y = 'TPR of Cox regression',
       size = 'Number of\nsimulations')

ggsave(file.path(resultDir, 'sim_tpr.pdf'),
       plot = p, width = 4, height = 3)
