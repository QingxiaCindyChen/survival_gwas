library('cowplot')
library('data.table')
library('ggplot2')

theme_set(theme_classic() +
            theme(axis.text = element_text(color = 'black'),
                  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))

resultDir = file.path('results', 'sim')

############################################################

simDataLogit = fread(file.path(resultDir, 'NPRLogitT.csv'))
simDataLogit[, trueModel := 'logistic']
simDataCox = fread(file.path(resultDir, 'NPRCoxT.csv'))
simDataCox[, trueModel := 'Cox']

simDataNeg = rbind(simDataLogit, simDataCox, fill = TRUE)
simDataNeg[, simIdx := 1:.N]

dNeg = melt(simDataNeg, id.vars = c('trueModel', 'simIdx'),
         measure.vars = colnames(simDataNeg)[2:7],
         variable.name = 'raw',
         value.name = 'tnr', variable.factor = FALSE)

dNeg[, fpr := 1 - tnr]
dNeg[, regressionModel := 'logistic']
dNeg[startsWith(raw, 'NPRCox'), regressionModel := 'Cox']

dNeg[, pvalCutoff := 0.01]
dNeg[endsWith(raw, '0.05'), pvalCutoff := 0.05]
dNeg[endsWith(raw, '0.1'), pvalCutoff := 0.1]

dNegSummary = dNeg[, .(meanFpr = mean(fpr)),
                   by = .(trueModel, regressionModel, pvalCutoff)]
setorderv(dNegSummary, c('trueModel', 'regressionModel', 'pvalCutoff'))

dNegSummary[, .(diffFpr = diff(rev(meanFpr))),
            by = .(trueModel, pvalCutoff)]

############################################################

simDataLogit = fread(file.path(resultDir, 'PPRLogitT0.csv'))
simDataLogit[, trueModel := 'logistic']
simDataCox = fread(file.path(resultDir, 'PPRCoxT0.csv'))
simDataCox[, trueModel := 'Cox']

simData = rbind(simDataLogit, simDataCox, fill = TRUE)
simData[, simIdx := 1:.N]

simData[, .(meanTprCox = mean(TPRCox0.05),
            meanTprLogit = mean(TPRLogit0.05),
            meanTprSeq = mean(TPRHybrid0.05)),
        by = trueModel]

d = melt(simData, id.vars = c('trueModel', 'simIdx'),
         measure.vars = colnames(simData)[2:10],
         variable.name = 'raw',
         value.name = 'tpr', variable.factor = FALSE)

d[, regressionModel := 'logistic']
d[startsWith(raw, 'TPRCox'), regressionModel := 'Cox']
d[startsWith(raw, 'TPRHybrid'), regressionModel := 'sequential']

d[, trueModelLabel := factor(trueModel, c('logistic', 'Cox'),
                             paste('Simulation model:', c('logistic', 'Cox')))]
d[, regressionModel := factor(regressionModel, c('logistic', 'Cox', 'sequential'))]

d[, pvalCutoff := 0.01]
d[endsWith(raw, '0.05'), pvalCutoff := 0.05]
d[endsWith(raw, '0.1'), pvalCutoff := 0.1]

setorderv(d, c('trueModelLabel', 'simIdx', 'regressionModel'))
d1 = d[, .(res = list(t.test(tpr[regressionModel == 'Cox'] - tpr[regressionModel == 'logistic']))),
       by = .(trueModelLabel, pvalCutoff)]

d1[, c('tpr_diff', 'ci_lo', 'ci_up') := list(res[[1]]$estimate, res[[1]]$conf.int[1],
                                             res[[1]]$conf.int[2]),
   by = 1:nrow(d1)]

############################################################

p1 = ggplot(d) +
  facet_wrap(~ trueModelLabel, nrow = 1) +
  geom_boxplot(aes(x = factor(pvalCutoff), y = tpr, fill = regressionModel),
               outlier.shape = NA, outlier.size = 0.7, alpha = 0.8, width = 0.6,
               position = position_dodge(width = 0.8)) +
  scale_y_continuous(limits = c(0.6, 0.93)) +
  scale_fill_manual(values = c('#33a02c', '#1f78b4', '#a6cee3')) +
  labs(x = 'Bonferroni-adjusted p-value cutoff', y = 'True positive rate',
       fill = 'Regression\nmodel')

p2 = ggplot(d1) +
  facet_wrap(~ trueModelLabel, nrow = 1) +
  geom_errorbar(aes(x = factor(pvalCutoff), ymin = ci_lo, ymax = ci_up), width = 0.2) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = 'Bonferroni-adjusted p-value cutoff', y = expression(TPR[Cox]~-~TPR[logistic]))

p = plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 0.85),
              align = 'v', axis = 'tblr', labels = 'AUTO')
ggsave(file.path(resultDir, 'sim_tpr.pdf'), plot = p, width = 5.25, height = 4.5)
