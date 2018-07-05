library('tidyverse')
library('cowplot')
library('lubridate')
library('survival')
library('snpStats')
library('survminer')

#filepath = file.path('Dropbox (VUMC)/', 'JakeColab_CountingGWAS', 'counting_gwas', 'analyze')
filepath = file.path('..')

buffer = 1 # years
minMaf = 0.01
nEvents = 3

# load snp data
exome = read.plink(file.path(filepath, 'genotype_data', 'exome', 'Exome_GRID_Euro'))
exomeSummary = col.summary(exome$genotypes)
idx = exomeSummary$MAF >= minMaf
snpInfo = read_csv(file.path(filepath, 'processed', 'exome_map.csv'), col_types=cols())
snpNames = intersect(colnames(exome$genotypes)[idx], snpInfo$snp.name[snpInfo$chromosome <= 22])

# load grid data
gridInfo = read_csv(file.path(filepath, 'processed', 'grid_info.csv'), col_types='ccDTT')
colnames(gridInfo) = tolower(colnames(gridInfo))
gridInfo = gridInfo %>%
  rename(gender = gender_epic) %>%
  mutate(first_entry_date = as.Date(first_entry_date),
         last_entry_date = as.Date(last_entry_date),
         first_age = time_length(first_entry_date - dob, 'years'),
         last_age = time_length(last_entry_date - dob, 'years'))

makeDf = function(d) {
  tt = c(d$first_age[1], d$age[!is.na(d$age)], d$last_age[1])
  df = tibble(age1 = tt[1:(length(tt)-1)],
              age2 = tt[2:length(tt)],
              status = c(rep_len(1, length(tt)-2), 0))}

# load phenotype data
phenoRaw = read_csv(file.path(filepath, 'processed', 'pheno_alzheimers.csv'), col_types='ccT')
colnames(phenoRaw) = tolower(colnames(phenoRaw))
phenoRaw$entry_date = as.Date(phenoRaw$entry_date)
  
# prepare phenotype data for analysis
pheno = phenoRaw %>%
  right_join(gridInfo, by='grid') %>%
  mutate(age = time_length(entry_date - dob, 'years')) %>%
  distinct(grid, first_age, last_age, age) %>%
  group_by(grid, first_age, last_age) %>%
  arrange(age) %>%
  filter(row_number() <= nEvents) %>%
  ungroup()

# run survival analysis
coxBase = tibble(grid = exome$fam$pedigree, sex = exome$fam$sex) %>%
  inner_join(pheno, by='grid') %>%
  group_by(grid, sex) %>%
  do(makeDf(.)) %>%
  ungroup() %>%
  select(grid, sex, age1, age2, status)
  
coxBase = coxBase %>%
  filter(age1 < age2)

#read in results already processed by Jake
results = read_csv(file.path(filepath, 'results', 'surv2_alzheimers.csv'))
snpName = results$snpName[1]

#Take the top hit for plotting purposes
coxInput = coxBase

coxInput$snp = as(exome$genotypes[coxBase$grid, snpName], 'numeric')[,1]
coxFit = coxph(Surv(age1, age2, status) ~ snp + sex + cluster(grid), data=coxInput)

# Create the new data for SNP status plotting
snp_df = with(coxInput,
              data.frame(
                         sex = rep(c('1', '2'), 3),
                         snp = rep(c(0, 1, 2), 2),
                         age2 = rep(mean(age2, na.rm = TRUE), 6))) #Using reference average age

#Predictions on dummy df
preds = survfit(coxFit, newdata = snp_df)

#Create survival plot
survivalplot = ggsurvplot(preds, data = coxInput, conf.int = TRUE, 
                          legend.labs = c('Male-0', 'Female-1', 'Male-2', 'Female-0', 'Male-1', 'Female-2'),
                                          ggtheme = theme_minimal())

survivalplotout = survivalplot$plot + xlab('Age') + ylab('Onset Probability') + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.title=element_text(size=14),legend.text=element_text(size=12),title=element_text(size=14,face="bold"),plot.title = element_text(hjust=0.5))

ggsave(filename = "Example_SurvivalPlot_Alz.pdf", plot = survivalplotout, path = file.path(filepath, 'results'), width = 6, height = 5, units = c('in'))
