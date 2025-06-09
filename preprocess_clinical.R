setwd("~/Desktop/PhD_Project_related/TCGA EXPLORATORY ANALYSIS")

library(dplyr)
library(finalfit)

clinical = read.csv(file="TCGA-READ_clinical.csv",header = T,stringsAsFactors = F)

## Cleaning data
NA_fifty = dim(clinical)[1] / 2 ## determining number of samples

NA_sum = colSums(is.na(clinical)) ## number of NA values per variable of DF
NA_sum = as.data.frame(NA_sum) ## converting to data frame
NA_sum = tibble::rownames_to_column(NA_sum, "variables") ## just setting rowname as the column named variables

NA_sum = NA_sum %>% filter(NA_sum < NA_fifty) ## keeping only those variables with less than 50% missing

clin_clean = clinical %>% select(one_of(NA_sum$variables)) ## filtered to keep those variables with less than 50% missing


## Removing duplicate observations
clin_clean0 = clin_clean %>% distinct_at('submitter_id', .keep_all = TRUE) ## df of variables with less than 50% missing

library(skimr)
## Remove numeric variables with unique observations
## Below command just summarises number of missing samples,mean, sd and percentile values
clin_clean_numeric_summary=clin_clean0 %>%
  select_if(is.numeric) %>%
  skim()


## remove character variables with unique observations
## summary of missing and complete rate for categorical variables
clin_clean_summary=clin_clean0 %>%
  select_if(is.character) %>%
  skim()


## Remove variables with similar information
not_missing_character = subset(clin_clean_summary,clin_clean_summary$n_missing==0)

intersect(colnames(clin_clean0),not_missing_character$skim_variable)
## Now make data frame with low missing or not missing variables 
## 25 columns from clin_clean0 are still in clin_clean0_1
clin_clean0_1 = clin_clean0 %>% 
  select(c(not_missing_character$skim_variable,clin_clean_numeric_summary$skim_variable))

intersect(colnames(clin_clean0_1),colnames(clin_clean0))

## Skipping the changing variable name part cause variables names
## are in R format with underscores


## Taming the data
## Use lubridate for dates
## basically just renaming the columns andconverting to the type of object
## Example renaming year_of_diagnosis as year_diagnose and converting to integer
clin_clean0_2 = clin_clean0_1 %>%
  mutate_if(is.character,as.factor) %>%
  mutate(patient_id = as.character(submitter_id),
         age = as.integer(age_at_diagnosis),
         year_diagnose = as.integer(year_of_diagnosis))


class(clin_clean0_1$submitter_id) ## character
class(clin_clean0_2$submitter_id) ## factor

class(clin_clean0_1$age_at_diagnosis) ## integer
class(clin_clean0_2$age_at_diagnosis) ## integer

class(clin_clean0_1$year_of_diagnosis) ## integer
class(clin_clean0_2$year_of_diagnosis) ## integer


## Checking NA pattern, such as MCAR, MAR MNAR

clin_clean0_2 %>%
  missing_plot()

clin_clean0_2_miss_glimpse=missing_glimpse(clin_clean0_2) %>% knitr::kable(.)
## clin_clean0_2_miss_glimpse is of class knitr_kable giving number and percentage of missing values of columns in table

## Checking numeric variables
clin_clean0_2 %>%
  select_if(is.numeric) %>%
  summary()


ggplot(clin_clean0_2,aes(age)) + geom_histogram(bins=20, alpha=0.8, color="red")

ggplot(clin_clean0_2,aes(year_diagnose)) + geom_histogram(bins=20, alpha=0.8, color="green")


ggplot(clin_clean0_2,aes(age_at_diagnosis)) + geom_histogram(bins=20, alpha=0.8, color="blue")


## Checking categorical variables
clin_clean0_2 %>%
  select_if(is.factor) %>%
  summary()


## aggregating levels
## They seem distinct enough to not aggregate
## This commmand gives a summary of number of missing values per categorical and factor variable
## Additionally it also gives summary of mean, sd and percentile for numeric variables
skim(clin_clean0_2)


write.csv(clin_clean0_2,file="clinical_data_cleaned.csv",col.names = T)
