library('tidyverse')
library('cowplot')
library('lubridate')
library('survival')
library('survminer')

#filepath = file.path('Dropbox (VUMC)/', 'JakeColab_CountingGWAS', 'counting_gwas', 'analyze')
filepath = file.path('..')

gridInfo = read_csv(file.path(filepath, 'processed', 'grid_info.csv'), col_types='ccDTT')
colnames(gridInfo) = tolower(colnames(gridInfo))
gridInfo = gridInfo %>%
	mutate(gender = gender_epic,
			 first_entry_date = as.Date(first_entry_date),
			 last_entry_date = as.Date(last_entry_date),
			 first_age = time_length(first_entry_date - dob, 'years'),
			 last_age = time_length(last_entry_date - dob, 'years')) %>%
	filter(last_age - first_age >= 1)

# filter(age_at_first_event >= 0, age_at_last_event <= 90)

phenoRaw = read_csv(file.path(filepath, 'processed', 'atrial_fibrillation.csv'), col_types='ccTi')
colnames(phenoRaw) = tolower(colnames(phenoRaw))
phenoRaw$entry_date = as.Date(phenoRaw$entry_date)

makeDf = function(d) {
	tt = c(d$first_age[1], d$age[!is.na(d$age)], d$last_age[1])
	df = tibble(time1 = tt[1:(length(tt)-1)],
					time2 = tt[2:length(tt)],
					status = c(rep_len(1, length(tt)-2), 0))}

pheno = phenoRaw %>%
	distinct(grid, entry_date) %>%
	right_join(gridInfo, by='grid') %>%
	mutate(age = time_length(entry_date - dob, 'years')) %>%
	select(grid, gender, age, first_age, last_age) %>%
	group_by(grid) %>%
	do(makeDf(.)) %>%
	inner_join(select(gridInfo, grid, gender, dob), by='grid')

phenoTmp = pheno %>%
	filter(time1 < time2)

#Cox proportional hazard model
#Hazard function - h(t|theta,xi) = h0(t)exp(t(xi)*theta) - baseline hazard (h0(t)) is assumed constant, and in a proportional hazard model, is left unspecified in this model
#We assume nothing about the baseline hazard, and that the predictor (sex) produces a proportional change in hazard regardless of time
m1 = coxph(Surv(time1, time2, status) ~ gender + cluster(grid), data=phenoTmp)
summary(m1) #exp(coef) = 1.45 - hazard multiplied by 1.45 for males - note that if cluster(grid is removed, result is the same)

# Create the new data and plot
sex_df = with(phenoTmp, data.frame(gender = c('M', 'F'), time2 = rep(mean(time2, na.rm = TRUE), 2)))
preds = survfit(m1, newdata = sex_df)
#Very slow plotting function
survivalplot = ggsurvplot(preds, data = sex_df, conf.int = TRUE, legend.labs=c("Sex=M", "Sex=F"),
           ggtheme = theme_minimal())
survivalplot$plot + xlab('Age') + ylab('Onset Probability') + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.title=element_text(size=14),legend.text=element_text(size=12),title=element_text(size=14,face="bold"),plot.title = element_text(hjust=0.5))



#How does the coefficient for this model compare to other types - Poisson process
library(mgcv)

#############################
#Some of this is taken from:
#http://data.princeton.edu/wws509/notes/c7.pdf - Survival Models
#https://github.com/paul-buerkner/brms/issues/230 - Bayesian approach to survival analysis
#http://dustintran.com/blog/survival-analysis-counting-processes-and-cox-models - Equivalency of coxph and poisson
#############################

#Modeling a poisson process - if baseline hazard is assumed constant, this takes the form p(censoredi|ti, xi) = ti*h0*exp(xi*theta)
#Essentially we are modeling the status of the censor at time ti given the time and predictor variable
#Model below: probability of censor status given a time offset (fixes coefficient to 1) plus time and covariate
phenoTmp$interval = phenoTmp$time2 - phenoTmp$time1

#It can be shown the poisson takes the form log(h) = log(h0) + log(t) + x
#log of time interval is a normalized fixed offset, and time2 is a linear function of time (this could be made into e.g. a spline)
poissonmodel <- gam(status ~ offset(log(interval)) + time2 + gender, data = phenoTmp, family = 'poisson')
#Hazard ratio for gender in the poisson model - 1.49 - pretty close to the coxph

#Model a spline to the time - essentially adjusting the intercept of the baseline hazard
poissonmodel2 <- gam(status ~ offset(log(interval)) + s(time2) + gender, data = phenoTmp, family = 'poisson')
#Hazard ratio for gender in the poisson model - 1.44 - almost identical to the coxph
