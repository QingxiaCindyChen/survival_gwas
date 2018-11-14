source(file.path('scripts', 'setup_regression.R'))
library('viridis')

d = readRDS('results/exome/pilot/examples_surv.rds')

p = ggplot(d) +
  facet_wrap(~ pheTitle + snp, nrow = 2, scales = 'free_y') +
  geom_step(aes(x = time, y = 1 - surv, color = factor(value, levels = 2:0)), size = 1) +
  labs(x = 'Age (y)', y = 'Adjusted cumulative hazard', color = 'Minor allele\ncount') +
  scale_color_viridis(direction = 1, discrete = TRUE)

ggsave(file.path('results/exome/pilot/examples_survival_ashg_poster.pdf'),
       plot = p, width = 7, height = 4.25)

