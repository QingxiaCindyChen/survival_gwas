library('tidyverse')
library('cowplot')
library('lubridate')
library('survival')

gridInfo = read_csv(file.path('processed', 'grid_info.csv'), col_types='ccDTT')
colnames(gridInfo) = tolower(colnames(gridInfo))
gridInfo = gridInfo %>%
	mutate(gender = gender_epic,
			 first_entry_date = as.Date(first_entry_date),
			 last_entry_date = as.Date(last_entry_date),
			 first_age = time_length(first_entry_date - dob, 'years'),
			 last_age = time_length(last_entry_date - dob, 'years')) %>%
	filter(last_age - first_age >= 1)

# filter(age_at_first_event >= 0, age_at_last_event <= 90)

phenoRaw = read_csv(file.path('processed', 'atrial_fibrillation.csv'), col_types='ccTi')
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

m1 = coxph(Surv(time1, time2, status) ~ gender + cluster(grid), data=phenoTmp)
summary(m1)

