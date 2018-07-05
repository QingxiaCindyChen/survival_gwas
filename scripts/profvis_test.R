library('profvis')

pf = profvis({
  source('scripts/run_regression.R')
})

saveRDS(pf, 'results/exome_profvis.rds')

# pf = readRDS('results/exome_profvis.rds')
# print(pf)
