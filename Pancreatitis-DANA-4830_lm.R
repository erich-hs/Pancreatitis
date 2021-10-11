pacman::p_load(tidyverse, olsrr, forecast, corrr, caret, GGally,
               lmtest, car, rsample, class, lime, reshape2, ggpubr, usethis)

setwd('~/R/DANA-4830/Assignment/Pancreatitis')

lm_dfs <- list()

lm_dfs[['Original_df']] <- read.csv('data/original_df.csv', stringsAsFactors = TRUE)
lm_dfs[['MICE_df']] <- read.csv('data/MICE_df.csv', stringsAsFactors = TRUE)
lm_dfs[['Amelia_df']] <- read.csv('data/Amelia_df.csv', stringsAsFactors = TRUE)
lm_dfs[['missForest_df']] <- read.csv('data/missForest_df.csv', stringsAsFactors = TRUE)

attach(lm_dfs)

##### Linear Regression #####
# Scenario 1: Number of pex treatment varies according to the exam results, age, gender, and disease severity?
#  •	Dependent Variable
# o	dt_pex_lan - number of PEX treatment
#  •	Independent Variable
# o	Age
# o	Gender
# o	ls_diem_apache_t0 – APACHE 2 score, disease severity
# o	Cls – group - the group of variables with the subclinical examination results.

# MICE
resp_var = 'dt_pex_lan'
indep_var = c('Age', 'Gender', 'ls_diem_apache_t0', names(Original_df[ , grepl("cls_", names(Original_df))]))
sc1_mice_model <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = MICE_df)
summary(sc1_mice_model)

indep_var = c('Gender', 'cls_hh_hct_t72', 'cls_hh_aptt_t0', 'cls_hh_aptt_t6', 'cls_hh_fib_t30',
              'cls_sh_chol_t0', 'cls_sh_tri_t6', 'cls_sh_tri_t72', 'cls_sh_na_t30', 'cls_km_paco2_t72',
              'cls_km_pao2_t54', 'cls_km_hco3_t6', 'cls_km_hco3_t30', 'cls_km_hco3_t54', 'cls_km_pf_t0')
sc1_mice_model2 <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = MICE_df)
summary(sc1_mice_model2)

indep_var = c('Gender', 'cls_hh_hct_t72', 'cls_hh_aptt_t6', 'cls_hh_fib_t30',
              'cls_sh_chol_t0', 'cls_sh_tri_t6', 'cls_sh_tri_t72', 'cls_sh_na_t30', 'cls_km_paco2_t72',
              'cls_km_pao2_t54', 'cls_km_hco3_t6', 'cls_km_hco3_t30')
sc1_mice_model3 <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = MICE_df)
mice_sc1 <- summary(sc1_mice_model3)

# Amelia
resp_var = 'dt_pex_lan'
indep_var = c('Age', 'Gender', 'ls_diem_apache_t0', names(Original_df[ , grepl("cls_", names(Original_df))]))
sc1_amelia_model <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = Amelia_df)
summary(sc1_amelia_model)

indep_var = c('Gender', 'cls_sh_ure_t30', 'cls_sh_cre_t6', 'cls_sh_ck_t0', 'cls_sh_na_t30',
              'cls_sh_ka_t0', 'cls_km_paco2_t72')
sc1_amelia_model2 <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = Amelia_df)
summary(sc1_amelia_model2)

# Adding best performed variables from MICE model 3
indep_var = c('Gender', 'cls_sh_ure_t30', 'cls_sh_na_t30', 'cls_sh_ka_t0', 'cls_hh_aptt_t6',
              'cls_hh_fib_t30', 'cls_hh_fib_t72', 'cls_sh_chol_t0', 'cls_sh_tri_t6', 'cls_km_pao2_t54')
sc1_amelia_model3 <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = Amelia_df)
summary(sc1_amelia_model3)

indep_var = c('Gender', 'cls_sh_na_t30', 'cls_sh_ka_t0',
              'cls_hh_fib_t30', 'cls_hh_fib_t72', 'cls_sh_chol_t0', 'cls_km_pao2_t54')
sc1_amelia_model4 <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = Amelia_df)
summary(sc1_amelia_model4)

indep_var = c('Gender', 'cls_sh_na_t30', 'cls_sh_ka_t0', 'cls_hh_fib_t72', 'cls_sh_chol_t0', 'cls_km_pao2_t54')
sc1_amelia_model5 <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = Amelia_df)
amelia_sc1 <- summary(sc1_amelia_model5)

# missForest
resp_var = 'dt_pex_lan'
indep_var = c('Age', 'Gender', 'ls_diem_apache_t0', names(Original_df[ , grepl("cls_", names(Original_df))]))
sc1_forest_model <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = missForest_df)
summary(sc1_forest_model)

indep_var = c('cls_hh_hct_t0', 'cls_hh_hc_t6', 'cls_hh_pt_t30', 'cls_hh_aptt_t72', 'cls_hh_fib_t72',
              'cls_sh_ure_t6', 'cls_sh_na_t30', 'cls_sh_ka_t0', 'cls_km_paco2_t30', 'cls_km_lac_t6')
sc1_forest_model2 <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = missForest_df)
summary(sc1_forest_model2)

# Adding best performed variables from MICE model 3
indep_var = c('cls_hh_fib_t30', 'cls_sh_chol_t0', 'cls_sh_na_t30', 'cls_km_paco2_t72',
              'cls_hh_aptt_t72', 'cls_hh_fib_t72')
sc1_forest_model3 <- lm(paste('dt_pex_lan ~ ', paste(indep_var, collapse = "+"), sep = ""), data = missForest_df)
missForest_sc1 <- summary(sc1_forest_model3)

# Removing underperforming models
rm(sc1_amelia_model, sc1_amelia_model2, sc1_amelia_model3, sc1_amelia_model4,
   sc1_mice_model, sc1_mice_model2, sc1_forest_model, sc1_forest_model2)

Adj.R.Squared <- c(mice_sc1$adj.r.squared, amelia_sc1$adj.r.squared, missForest_sc1$adj.r.squared)
Indep.Variables <- c(mice_sc1$fstatistic['numdf'], amelia_sc1$fstatistic['numdf'], missForest_sc1$fstatistic['numdf'])
P.Values <- c(pf(mice_sc1$fstatistic['value'], mice_sc1$fstatistic['numdf'], mice_sc1$fstatistic['dendf'], lower.tail = FALSE),
              pf(amelia_sc1$fstatistic['value'], amelia_sc1$fstatistic['numdf'], amelia_sc1$fstatistic['dendf'], lower.tail = FALSE),
              pf(missForest_sc1$fstatistic['value'], missForest_sc1$fstatistic['numdf'], missForest_sc1$fstatistic['dendf'], lower.tail = FALSE))

scenario1 <- data.frame(Model = c('MICE', 'Amelia', 'missForest'), Adj.R.Squared, Indep.Variables, P.Values)
scenario1


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

# MICE
sc3_mice_model <- lm(rv_ngaydt ~ ., data = MICE_df[, -124])
summary(sc3_mice_model)

indep_var = c('ls_diem_ranson_t0', 'cls_hh_tc_t30', 'cls_hh_pt_t0', 'cls_hh_aptt_t30', 'cls_hh_fib_t6', 
              'cls_sh_cre_t72', 'cls_sh_ck_t0', 'cls_sh_chol_t0', 'cls_sh_amy_t0', 'cls_km_pao2_t72',
              'cls_km_lac_t30', 'dt_dich_vao_t48', 'dt_dich_ra_t48', 'dt_dich_bilan_t24', 'dt_nhin_ngay')
sc3_mice_model2 <- lm(paste('rv_ngaydt ~ ', paste(indep_var, collapse = "+"), sep = ""), data = MICE_df)
summary(sc3_mice_model2)

indep_var = c('cls_hh_tc_t30', 'cls_hh_aptt_t30', 'cls_sh_ck_t0', 'dt_dich_vao_t48', 'dt_nhin_ngay')
sc3_mice_model3 <- lm(paste('rv_ngaydt ~ ', paste(indep_var, collapse = "+"), sep = ""), data = MICE_df)
mice_sc3 <- summary(sc3_mice_model3)

# Amelia
sc3_amelia_model <- lm(rv_ngaydt ~ ., data = Amelia_df[, -117])
summary(sc3_amelia_model)

indep_var = c('cls_sh_tri_t6', 'cls_sh_na_t30', 'cls_km_paco2_t54', 'cls_km_pao2_t6', 'dt_nhin_ngay')
sc3_amelia_model2 <- lm(paste('rv_ngaydt ~ ', paste(indep_var, collapse = "+"), sep = ""), data = Amelia_df)
summary(sc3_amelia_model2)

sc3_amelia_model3 <- lm(rv_ngaydt ~ cls_sh_tri_t6 + dt_nhin_ngay, data = Amelia_df)
amelia_sc3 <- summary(sc3_amelia_model3)

# missForest
sc3_forest_model <- lm(rv_ngaydt ~ ., data = missForest_df[, -124])
summary(sc3_forest_model)

indep_var = c('stomachache', 'ls_diem_ranson_t0', 'cls_sh_tri_t6', 'cls_km_hco3_t54', 'dt_dich_ra_t24')
sc3_forest_model2 <- lm(paste('rv_ngaydt ~ ', paste(indep_var, collapse = "+"), sep = ""), data = missForest_df)
summary(sc3_forest_model2)

sc3_forest_model3 <- lm(rv_ngaydt ~ ls_diem_ranson_t0 + cls_sh_tri_t6 + cls_km_hco3_t54, data = missForest_df)
missForest_sc3 <- summary(sc3_forest_model3)

# Removing underperforming models
rm(sc3_amelia_model, sc3_amelia_model2, sc3_forest_model, sc3_forest_model2, sc3_mice_model, sc3_mice_model2)

Adj.R.Squared <- c(mice_sc3$adj.r.squared, amelia_sc3$adj.r.squared, missForest_sc3$adj.r.squared)
Indep.Variables <- c(mice_sc3$fstatistic['numdf'], amelia_sc3$fstatistic['numdf'], missForest_sc3$fstatistic['numdf'])
P.Values <- c(pf(mice_sc3$fstatistic['value'], mice_sc3$fstatistic['numdf'], mice_sc3$fstatistic['dendf'], lower.tail = FALSE),
              pf(amelia_sc3$fstatistic['value'], amelia_sc3$fstatistic['numdf'], amelia_sc3$fstatistic['dendf'], lower.tail = FALSE),
              pf(missForest_sc3$fstatistic['value'], missForest_sc3$fstatistic['numdf'], missForest_sc3$fstatistic['dendf'], lower.tail = FALSE))

scenario3 <- data.frame(Model = c('MICE', 'Amelia', 'missForest'), Adj.R.Squared, Indep.Variables, P.Values)
scenario3
