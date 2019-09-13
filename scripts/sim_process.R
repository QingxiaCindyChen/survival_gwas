library('data.table')

load('cindy_1/PPRCoxT0_hybrid.Rdata')
fwrite(ResultAll, 'results/sim/PPRCoxT0.csv')

load('cindy_1/PPRLogitT0_hybrid.Rdata')
fwrite(ResultAll, 'results/sim/PPRLogitT0.csv')

load('cindy_1/NPRCoxT.Rdata')
fwrite(ResultAll, 'results/sim/NPRCoxT.csv')

load('cindy_1/NPRLogitT.Rdata')
ResultAll = as.data.table(ResultAll)
setnames(ResultAll, c('EventRate', 'NPRLogit0.01', 'NPRLogit0.05', 'NPRLogit0.1',
                      'NPRCox0.01', 'NPRCox0.05', 'NPRCox0.1'))
fwrite(ResultAll, 'results/sim/NPRLogitT.csv')
