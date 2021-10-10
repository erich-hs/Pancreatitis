pacman::p_load(tidyverse, olsrr, forecast, corrr, caret, GGally,
               lmtest, car, rsample, class, lime, reshape2, ggpubr, usethis)

setwd('~/R/DANA-4830/Assignment/Pancreatitis')

lm_dfs <- list()
?read.csv
lm_dfs[['Original_df']] <- read.csv('data/original_df.csv', stringsAsFactors = TRUE)
lm_dfs[['MICE_df']] <- read.csv('data/MICE_df.csv', stringsAsFactors = TRUE)
lm_dfs[['Amelia_df']] <- read.csv('data/Amelia_df.csv', stringsAsFactors = TRUE)
lm_dfs[['missForest_df']] <- read.csv('data/missForest_df.csv', stringsAsFactors = TRUE)

##### Linear Regression #####
# Scenario 2: Do the hereditary information have any impact in the Severity of the disease and in the subclinical examination?
#  •	Dependent Variable
# o	ls_diem_apache_t0 – APACHE 2 score, disease severity
#
#  •	Independent Variable
# o	Age
# o	Gender
# o	ts_giadinh – Historical problems
# o	ts_benhmat – Gallblader Problems
# o	ts_ruou – Drinking Problems
# o	ts_dtd – Diabetes Problem

resp_var = 'ls_diem_apache_t0'
indep_var = c('ts_giadinh', 'ts_benhmat', 'ts_ruou', 'ts_dtd')

attach(lm_dfs)

# MICE
sc2_mice_model <- lm(ls_diem_apache_t0 ~ ts_giadinh + ts_benhmat + ts_ruou + ts_dtd, data = MICE_df)
summary(sc2_mice_model)

# Amelia
sc2_amelia_model <- lm(ls_diem_apache_t0 ~ ts_giadinh + ts_benhmat + ts_ruou + ts_dtd, data = Amelia_df)
summary(sc2_amelia_model)

# missForest
sc2_missForest_model <- lm(ls_diem_apache_t0 ~ ts_giadinh + ts_benhmat + ts_ruou + ts_dtd, data = missForest_df)
summary(sc2_missForest_model)

# Scenario 2 doesn't make sense

# Scenario 3: What factors most influence the duration of the patient in the hospital?
#  •	Dependent Variable
# o	rv_ngaydt – Duration in hospital in Days
#
#  •	Independent Variable
# o	Age
# o	Gender
# o	Cls – group
# o	Dt – group
# o	Ts – group


