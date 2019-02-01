library('BEDMatrix')
library('cowplot')
library('data.table')
library('doParallel')
library('readr')
library('lubridate')
library('qqman')
library('survival')
library('yaml')

theme_set(theme_light() +
            theme(axis.text = element_text(color = 'black'),
                  strip.text = element_text(color = 'black'),
                  panel.grid.minor = element_blank(),
                  legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')))

############################################################

pcData = read_csv('processed/mega/pc_data.csv.gz')
setDT(pcData)

covarData = read_delim('params/mega/MEGA3_forPheWAS_01_Ws.cov', delim = ' ')
setDT(covarData)
covarData = covarData[, .(grid = FID, b1, b2, b3, b4, b5, b6,
                          b7, b8, b9, b10, b11, b12, b13, b14)]
covarData1 = melt(covarData, id.vars = 'grid', variable.factor = FALSE)
covarData2 = covarData1[, .(batch = variable[value == 1]), by = grid]

############################################################

plot(pcRaw$eigenval[1:10]^2 / sum(pcRaw$eigenval^2))
plot(pcRaw$eigenval[1:20])

d = merge(pcData, covarData2, by = 'grid')
d[, batchAgg := ifelse(batch %in% c('b13', 'b14'), 'b12', batch)]

p12 = ggplot(d) +
  geom_point(aes(x = PC1, y = PC2, color = batchAgg),
             shape = 16, size = 0.5, alpha = 0.5) +
  scale_color_brewer(type = 'qual', palette = 'Paired')

p13 = ggplot(d) +
  geom_point(aes(x = PC1, y = PC3, color = batchAgg),
             shape = 16, size = 0.5, alpha = 0.5) +
  scale_color_brewer(type = 'qual', palette = 'Paired')

p24 = ggplot(d) +
  geom_point(aes(x = PC2, y = PC4, color = batchAgg),
             shape = 16, size = 0.5, alpha = 0.5) +
  scale_color_brewer(type = 'qual', palette = 'Paired')

p25 = ggplot(d) +
  geom_point(aes(x = PC2, y = PC5, color = batchAgg),
             shape = 16, size = 0.5, alpha = 0.5) +
  scale_color_brewer(type = 'qual', palette = 'Paired')

p = plot_grid(p12, p13, p24, p25, nrow = 2)
p
