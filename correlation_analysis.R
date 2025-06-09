setwd("~/Desktop/PhD_Project_related/TCGA EXPLORATORY ANALYSIS")


library(tidyverse)
library(skimr)
library(finalfit)
library(rstatix)
library(ggpubr)
library(GGally)
library(plotly)


## Importing the data
read_clin = read.csv(file="clinical_data_cleaned.csv",header = T,stringsAsFactors = F)


## Taming the data
read_clin = read_clin %>%
  mutate_if(is.character,as.factor) %>%
  mutate(patient_id = as.character(patient_id),
         age = as.integer(age),
         year_diagnose = as.integer(year_diagnose))

## check
glimpse(read_clin)


## The dependent variable
table(read_clin$vital_status)


## Numeric variables vs vital status
cols_numeric = read_clin %>%
  select_if(is.numeric) %>%
  names

read_clin_numeric = read_clin %>%
  select(one_of(cols_numeric,"vital_status"))

read_clin_numeric=read_clin_numeric[complete.cases(read_clin_numeric),]

## This plot a correlation matrix of all the continuous variables for alive and dead subjects
ggpairs(read_clin_numeric,columns = cols_numeric,
        title = "Correlation matrix",
        mapping = aes(colour = vital_status),
        upper = list(combo = wrap("box_no_facet", alpha=0.1),
                     continuous = wrap("cor",size=2,alignPercent = 0.8)),
        lower = list(continuous = wrap("smooth",alpha=0.2,smooth = 0.2))) +
  theme(panel.background = element_rect(colour = "black",size=0.5,fill="white"),
        panel.grid.major = element_blank())


## Run multiple t tests on vital_status

## Transform the data into long format

## Put all variables in the same column except "vital_status", the grouping variable

levels(read_clin_numeric$vital_status) = c("Dead","Alive","Not_reported")


## Convert to long tidyverse format

read_clin_numeric.long = read_clin_numeric %>%
  pivot_longer(-vital_status, names_to = "variables", values_to="value")

read_clin_numeric.long = read_clin_numeric %>%
  pivot_longer(-vital_status, names_to = "variables", values_to = "value")

read_clin_numeric.long = read_clin_numeric.long[!is.na(read_clin_numeric.long$value),]
read_clin_numeric.long$value.log = log2(read_clin_numeric.long$value+1)

read_clin_numeric.long %>% slice_sample(n=6) %>% knitr::kable(.)


## Now group the data by variables and compare over survival status
## Adjust the p-values and add significance levels

read_clin_numeric.long.no.na = read_clin_numeric.long %>%
  filter(value.log != "NA")

stat.test = read_clin_numeric.long.no.na %>%
  group_by(variables) %>%
  t_test(value ~ vital_status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test %>% knitr::kable(.)



## Create the plot on logscale
myplot = ggboxplot(
  read_clin_numeric.long.no.na,x="vital_status", y= "value.log",
  fill="vital_status", palette = "npg", legend="none",
  ggtheme = theme_pubr(border = TRUE)
) + facet_wrap(~variables)

## Add statistical test p-values
## OBS: different p-values over vaule vs log value
stat.test = stat.test %>% add_xy_position(x="vital_status")
myplot + stat_pvalue_manual(stat.test, label="p.adj.signif")



## Group the data by variables and do a graph for each variable
graphs = read_clin_numeric.long.no.na %>%
  group_by(variables) %>%
  doo(
    ~ggboxplot(
      data=., x="vital_status", y="value",
      fill="vital_status",palette = "npg", legend="none",
      ggtheme = theme_pubr()
    ) +
      geom_jitter(width = 0.05, alpha=0.2, color="orange"),
    result="plots"
  )

graphs %>% knitr::kable(.)



## Summary for continuous variables
explanatory_num = read_clin %>%
  select(-c(vital_status,days_to_diagnosis)) %>%
  select_if(is.numeric) %>%
  names
## Based on plot above, it seems that there are several correlations

dependent = "vital_status"

## Threw error due to not enough labels.
## now it does after removing days_to_diagnosis variable
table_num = read_clin %>%
  summary_factorlist(dependent, explanatory_num,p=TRUE,
                     add_dependent_label = TRUE,na_include = TRUE)
## This table just nicely summarises various metrics from table_num for alive, dead and not reported samples
## Also it returns p values for if there there is a significant difference between these variables across groups
knitr::kable(table_num, row.names = FALSE, align = c("l","l","r","r","r"))


## Correlation Matrix
## Pearson's or Spearman
corr_num = read_clin %>%
  select_if(is.numeric) %>%
  drop_na()

## Check correlatiuon between variables
cor_matrix = cor(corr_num, method = "spearman")
cor_matrix = round(cor_matrix,2)

## There aren't any large correlations with highest being -0.29


## Catergoricla varibles vs vital_status
explanatory_char = read_clin %>%
  select(-vital_status) %>%
  select_if(is.factor) %>%
  names

dependent = 'vital_status'

vars_use = c("tissue_or_organ_of_origin",
             "primary_diagnosis",
             "prior_malignancy",
             "classification_of_tumor",
             "race","gender","ethnicity")

## summary_factorlist
'''
function that takes a single dependent variable with a vector of explanatory 
variable names (continuous or categorical variables) to produce a summary table.
'''
table_char = read_clin %>%
  summary_factorlist(dependent,vars_use,p=TRUE,
                     add_dependent_label = TRUE, na_include = TRUE)


## Dropping levels with narrow distribution
## Group some levels or drop one when grouping is not possible
## These variables were missing from the dataset
read_clin2 = read_clin %>%
  mutate(neoplasm_stg = fct_collapse(tumor_stg, 'T1-T2' = c('T1','T2'), 'T3-T4'=c('T3','T4')),
         tumor_stg = fct_collapse(histology_grd, 'G1-G2' = c('G1','G2')))





  
  
