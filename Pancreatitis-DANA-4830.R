pacman::p_load(tidyverse, olsrr, forecast, corrr, caret, GGally,
               lmtest, car, rsample, class, lime, reshape2, ggpubr, usethis)

setwd('~/R/DANA-4830/Assignment/Pancreatitis')

df <- read.csv('data/APNotCleaned.csv')
wdf <- df ## Setting working dataframe for data screening @ reset point @

##### Data Screening #####
wdf <- rename(wdf, ID = names(df)[1],
              vv_reason1 = vv_reason_1,
              dt_pex_ranson_s_lan1 = dt_pex_ranson_s_lan3,
              stomachache = daubung,
              vomiting = non,
              dead = kq,
              complication = bcxa)
wdf <- wdf %>% mutate(across(where(is.character), str_trim))

# Removing empty variables
wdf <- select(wdf, -c(CLS_S2, CLS_S1, CLS_S0, DT_PE0))

options(max.print=200000)
summary(wdf)

##### Categorical variables #####
table(wdf$Age)
table(wdf$Gender)
wdf$Gender <- factor(wdf$Gender,
                     levels = c("Nam", "Nu"),
                     labels = c("Male", "Female"))

table(wdf$vv_reason1) # Hospital admission reason 1
wdf$vv_reason1 <- factor(wdf$vv_reason1,
                         levels = c("dau bung"),
                         labels = c("Stomachache"))

table(wdf$vv_reason2) # Hospital admission reason 2
wdf$vv_reason2 <- factor(wdf$vv_reason2,
                         levels = c("dau bung thuong vi"),
                         labels = c("Pain Relief"))

table(wdf$vv_reason3) # Hospital admission reason 3
wdf$vv_reason3 <- factor(wdf$vv_reason3,
                         levels = c("buon non", "non"),
                         labels = c("Nausea", "Nothing?"))

table(wdf$vv_Others) # Hospital admission reason 4
# Grouped similar categories
wdf$vv_Others <- factor(wdf$vv_Others,
                        levels = c("Dau bung man suon (P)", "dau bung quanh ron", "ha suon=man suon", 
                                   "kho tho", "vtc tang triglycerid", "VTC tang triglycerid , gian dai be than do soi NQ"),
                        labels = c("Abdominal Pain", "Abdominal Pain", "Lower Ribs", 
                                   "Shortness of Breath", "VTC Increase Triglycerides", "VTC Increase Triglycerides"))

table(wdf$ts_giadinh) # Hereditary information
wdf$ts_giadinh <- factor(wdf$ts_giadinh,
                         levels = c("co", "khong"),
                         labels = c("Yes", "No"))

table(wdf$details_ts_giadinh) # Hereditary information breakdown
wdf$details_ts_giadinh <- factor(wdf$details_ts_giadinh,
                                 levels = c("rl lipid", "rl lipid", "RLCH lipid", "RLCH lipid mau", "rlmm cach 2 nam"),
                                 labels = c("RL Lipid", "RL Lipid", "RLCH Lipid", "RLCH Rapid Lipid Methabolism",
                                            "RLMM 2 Years Apart"))

table(wdf$ts_benhmat) # Gallbladder problem
wdf$ts_benhmat <- factor(wdf$ts_benhmat,
                         levels = c("Co", "khong"),
                         labels = c("Yes", "No"))

table(wdf$ts_ruou) # Drinking problem
wdf$ts_ruou <- factor(wdf$ts_ruou,
                         levels = c("Co", "Khong"),
                         labels = c("Yes", "No"))

table(wdf$ts_dtd) # Diabetes problem
wdf$ts_dtd <- factor(wdf$ts_dtd,
                         levels = c("Co", "Khong"),
                         labels = c("Yes", "No"))

table(wdf$ts_vtc) # Cholecystitis problem
wdf$ts_vtc <- factor(wdf$ts_vtc,
                         levels = c("Co", "Khong"),
                         labels = c("Yes", "No"))

table(wdf$ts_vtc_lancuoi) # Last detection of cholecystitis problem
# Variable that needs major adjustment
# For "Years away" responses, it will be assumed from 2021
# Changed to numeric variable - considered corresponding year value of each date
wdf$ts_vtc_lancuoi <- dplyr::recode(wdf$ts_vtc_lancuoi, !!!
                                    setNames(c(2009, 2012, 2013, 2014, 2015, 2016, 2016, 2016, 2017, 2018, 2019, 2018,
                                               2016, 2020, 2021, 2011, 2019, 2021, 2018, 2021, 2021, 2016,
                                               2021, 2020, 2016, 2019, 2016, 2017,
                                               2018, 2017, 2019, 2017, 2018, 2017,
                                               2017, 2015, 2017, 2017, 2015, 2021),
                                             c("2009", "2012", "2013", "2014", "2015", "2016",
                                               "2016-01-14", "2016-02-07", "2017", "2018", "Apr-19", "Aug-18",
                                               "cach  5 nam", "cach 1 nam", "cach 1 thang", "cach 10 nam",
                                               "cach 2 nam", "cach 2 thang", "cach 3 nam", "cach 3 thang", "cach 4 thang",
                                               "cach 5 nam", "cach 8 tuan", "cach1 nam", "Feb-16", "Feb-19", "Jan-16", "Jan-17",
                                               "Jan-18", "Jul-17", "Jun-19", "Mar-17", "Mar-18", "May-17",
                                               "Nov-17", "Sep-15", "Sep-17", "T3/2017",
                                               "T8/2015", "truoc 1 thang")
                                             ))

table(wdf$stomachache) # Stomachache
wdf$stomachache <- factor(wdf$stomachache,
                      levels = c(0, 1),
                      labels = c("No", "Yes"))
as.numeric(wdf$stomachache)

table(wdf$vomiting) # Vomiting
wdf$vomiting <- factor(wdf$vomiting,
                  levels = c(0, 1),
                  labels = c("No", "Yes"))
as.numeric(wdf$vomiting)
wdf$vomiting[is.na(wdf$vomiting)] <- 'No' # Assigning 'No' to missing values

table(wdf$ls_cn_bidaitien) # Symptoms of defecation
wdf$ls_cn_bidaitien <- factor(wdf$ls_cn_bidaitien,
                              levels = c("khong", "t0", "T0", "t30", "T30"),
                              labels = c("No", "T0", "T0", "T30", "T30"))

table(wdf$ls_cn_ialong) # Symptoms of diarrhea
wdf$ls_cn_ialong <- factor(wdf$ls_cn_ialong,
                           levels = c("khong", "t0", "T0", "t6", "t96"),
                           labels = c("No", "T0", "T0", "T6", "T96"))

table(wdf$ls_tht_bungchuong) # Symptoms of abdominal distension
wdf$ls_tht_bungchuong <- factor(wdf$ls_tht_bungchuong,
                                levels = c("khong", "t0", "T0", "t30", "T30", "t6", "T6", "t96"),
                                labels = c("No", "T0", "T0", "T30", "T30", "T6", "T6", "T96"))

table(wdf$ls_tt_lungsuon) # Symptoms of painful pressure throughout the abdomen
# "Yes" responses considered = "T0"
wdf$ls_tt_lungsuon <- factor(wdf$ls_tt_lungsuon,
                             levels = c("co", "khong", "t0", "T0"),
                             labels = c("T0", "No", "T0", "T0"))

table(wdf$ls_tn_ha_t6) # Blood Pressure
wdf$ls_tn_ha_t6 <- factor(trimws(wdf$ls_tn_ha_t6, which = c("right")), exclude = "")

table(wdf$cls_sa_tuy_t0) # Pancreas Ultrasound
# Removed responses such as "No", "Hard to see", "unobservable", etc.
wdf$cls_sa_tuy_t0 <- factor(wdf$cls_sa_tuy_t0,
                            levels = c("dich quanh tuy", "tang kich thuoc tham nhieu", "tham nhiem dau tuy",
                                       "tham nhieu phu", "vtc", "VTC", "vtc hoai tu", "vtc phu", "VTC phu"),
                            labels = c("Ditch around the marrow", "Increase the size of the bruise a lot", "Optional oil infiltration",
                                       "Take a lot of sub", "VTC", "VTC", "VTC", "VTC", "VTC"))

table(wdf$cls_sa_dichob_t0) # Abdominal fluid ultrasound at T0
wdf$cls_sa_dichob_t0 <- factor(wdf$cls_sa_dichob_t0,
                               levels = c("Co", "Khong"),
                               labels = c("Yes", "No"))

table(wdf$cls_sa_mat_t0) # Bladder ultrasound at T0
wdf$cls_sa_mat_t0 <- factor(wdf$cls_sa_mat_t0,
                            levels = c("bt", "polyp tu"),
                            labels = c("BT", "Polyps"))

table(wdf$cls_sa_ketluan_t0) # Ultrasound Conclusion
wdf$cls_sa_ketluan_t0[8] # Problematic spelling with unwanted characters
wdf$cls_sa_ketluan_t0 <- factor(wdf$cls_sa_ketluan_t0,
                                levels = c("TD VTC", "TDVTC", "vtc", "VTC", "vtc hoai tu", "VTC hoai tu", wdf$cls_sa_ketluan_t0[8],
                                           "vtc phu", "VTC phu", "vtc the phu", "VTC the phu"),
                                labels = c("TD VTC", "TD VTC", "VTC", "VTC", "Necrotic VTC", "Necrotic VTC", "VTC",
                                           "Sub VTC", "Sub VTC", "VTC Cover", "VTC Cover"))

table(wdf$cls_ct_tuy_lan1) # Pancreas Computer Tomography
wdf$cls_ct_tuy_lan1 <- factor(wdf$cls_ct_tuy_lan1,
                              levels = c("vtc", "VTC", "vtc phu", "VTC phu", "VTC Phu", "VTC phu ne", "vtc the phu", "VTC the phu"),
                              labels = c("VTC", "VTC", "Sub VTC", "Sub VTC", "Sub VTC", "Sub VTC", "Sub VTC", "Sub VTC"))

table(wdf$cls_ct_dichob_lan1) # Abdominal fluid Computer Tomography
wdf$cls_ct_dichob_lan1 <- factor(wdf$cls_ct_dichob_lan1,
                                 levels = c("ci", "co", "Co", "khong", "Khong"),
                                 labels = c("Yes", "Yes", "Yes", "No", "No"))

table(wdf$cls_ct_balthazar_lan1) # Balthazar score with Computer Tomography
# A: normal pancreas: 0
# B: enlargement of pancreas: 1
# C: inflammatory changes in pancreas and peripancreatic fat: 2
# D: ill-defined single peripancreatic fluid collection: 3
# E: two or more poorly defined peripancreatic fluid collections: 4
# https://radiopaedia.org/articles/balthazar-score
wdf$cls_ct_balthazar_lan1 <- factor(wdf$cls_ct_balthazar_lan1,
                                    levels = c("A", "b", "c", "C", "d", "D", "e", "E"),
                                    labels = c("A", "B", "C", "C", "D", "D", "E", "E"))

## cls_sh_bil: The level of bilirubin in the blood may increase if the common bile duct is blocked due to inflammation in the pancreas.
table(wdf$cls_sh_bil_t0) # Variable needs clarification before treating - Direct and Indirect bilirubin
table(wdf$cls_sh_bil_t6) # Variable needs clarification before treating
table(wdf$cls_sh_bil_t30) # Variable needs clarification before treating
table(wdf$cls_sh_bil_t72) # Variable needs clarification before treating

## cls_sh_gan: ALT or AST levels more than three times the upper limit of normal indicates gallstones as the cause of acute pancreatitis.
table(wdf$cls_sh_gan_t0) # Variable needs clarification before treating
table(wdf$cls_sh_gan_t6) # Variable needs clarification before treating
table(wdf$cls_sh_gan_t30) # Variable needs clarification before treating

## cls_sh_ca: Calci total
table(wdf$cls_sh_ca_t0) # Variable needs clarification before treating

table(wdf$pex)
wdf$pex <- factor(wdf$pex,
                  levels = c(0, 1),
                  labels = c('Standard Treatment', 'PEX Treatment'))

table(wdf$dead) ## Song = ???

table(wdf$complication)
wdf$complication[is.na(wdf$complication)] <- 0 # Assigning 0 to NAs to convert to a binary variable


##### Numeric Variables #####
#Changing ID 122 dt_pex_hdl_t_lan1 to NA and column from categorical to numeric
wdf$dt_pex_hdl_t_lan1[122] <- NA
wdf$dt_pex_hdl_t_lan1 <- as.numeric(wdf$dt_pex_hdl_t_lan1)

#Changing ID's dt_pex_ldl_t_lan1 to NA and column from categorical to numeric
table(wdf$dt_pex_ldl_t_lan1)
wdf$dt_pex_ldl_t_lan1[wdf$dt_pex_ldl_t_lan1 == 'duc'] <- NA
wdf$dt_pex_ldl_t_lan1[wdf$dt_pex_ldl_t_lan1 == 'ht duc'] <- NA
wdf$dt_pex_ldl_t_lan1 <- as.numeric(wdf$dt_pex_ldl_t_lan1)

#Changing ID's dt_pex_sauvv to its respective and column from categorical to numeric
wdf$dt_pex_sauvv[24] <- 9
wdf$dt_pex_sauvv[86] <- 5
wdf$dt_pex_sauvv[138] <- 23
wdf$dt_pex_sauvv[39] <- 22
wdf$dt_pex_sauvv[6] <- 15
wdf$dt_pex_sauvv[66] <- 12
wdf$dt_pex_sauvv[139] <- 10
wdf$dt_pex_sauvv[67] <- 39
wdf$dt_pex_sauvv[20] <- 19
wdf$dt_pex_sauvv[122] <- 17
wdf$dt_pex_sauvv[58] <- 16
wdf$dt_pex_sauvv[150] <- 16
wdf$dt_pex_sauvv[75] <- 14
wdf$dt_pex_sauvv[52] <- 11
wdf$dt_pex_sauvv[125] <- 11
wdf$dt_pex_sauvv[25] <- 10
wdf$dt_pex_sauvv[56] <- 0.25
wdf$dt_pex_sauvv[35] <- 9
wdf$dt_pex_sauvv[121] <- 8
wdf$dt_pex_sauvv <- as.numeric(wdf$dt_pex_sauvv)

table(wdf$rv_ngaydt) # Duration of hospital stay in days

#Box Plot
ggplot(wdf, aes(Gender, rv_ngaydt)) + geom_boxplot()  + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun.y = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Duration in Hospital (days)", y = "Duration (days)")

#replacing the duration of the hospital in days for this observation from 2 to 3. Since this patient has all the records for the exams of 72 hours. 
wdf$rv_ngaydt[112] <- 3

table(wdf$ts_ruou_nam) # Drinking problem - Unit unknown
#Box Plot
ggplot(wdf, aes(Gender, ts_ruou_nam)) + geom_boxplot()+ stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Drinking Problem")


table(wdf$ts_ruou_nam_ml) # Drinking problem - mililiters?
#Box Plot
ggplot(wdf, aes(Gender, ts_ruou_nam_ml)) + geom_boxplot()+ stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Drinking Problem (ml)", y = "Drinking Problem (ml)")


table(wdf$ls_tt_alob_t0) # Symptom Abdominal pressure at t0
ggplot(wdf, aes(Gender, y=ls_tt_alob_t0)) +  geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Symptom Abdominal pressure at t0", y = "Abdominal pressure at t0")


table(wdf$ls_tt_bmi_t0) # Symptom BMI measure t0
ggplot(wdf, aes(Gender, ls_tt_bmi_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Symptom BMI measure t0", y = "BMI measure t0")
# Clearly an outlier at y = 258, the correct input should be 25.8 since its normal range is 18.5 - 24.9
wdf$ls_tt_bmi_t0[148] <- 25.8

#Cleaned Boxplot
ggplot(wdf, aes(Gender, ls_tt_bmi_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Symptom BMI measure t0 (Clean)", y = "BMI measure t0")
#The outliers in this case mean that the patient has 30.0 or higher BMI, it falls within the obese range.


table(wdf$ls_tn_mach_t0) # Symptom Heartbeat at t0
ggplot(wdf, aes(Gender, ls_tn_mach_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Symptom Heartbeat at t0", y = "Heartbeat at t0")
#Outlier 158, is from a 50 year old man, high heartbeat rate for his age, although possible.


table(wdf$ls_tn_nhiet_t0) # Body Temperature
ggplot(wdf, aes(Gender, ls_tn_nhiet_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Symptom Body Temperature at t0", y = "Body Temperaturet at t0")
# Clearly an outlier at y = 366 and 3.7

wdf$ls_tn_nhiet_t0[164] <- 36.6
wdf$ls_tn_nhiet_t0[61] <- 37

ggplot(wdf, aes(Gender, ls_tn_nhiet_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Symptom Body Temperature at t0 (Clean)", y = "Body Temperaturet at t0")


table(wdf$ls_tn_spo2_t0) # Saturation of peripheral oxygen
ggplot(wdf, aes(Gender, ls_tn_spo2_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Saturation of peripheral oxygen at t0", y = "Saturation of peripheral oxygen")
#Values under 90 are considered low, but possible, outliers kept

table(wdf$ls_tn_cvp_t0) # Center venus pressure
ggplot(wdf, aes(Gender, ls_tn_cvp_t0)) + geom_boxplot()+ stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="Center venus pressure at t0", y = "Center venus pressure")
# Possible outliers at y = 99. Expected range 3 - 10 mmHg. 
#Drop this variable?
#wdf = subset(wdf,select = -c(ls_tn_cvp_t0))

table(wdf$ls_diem_apache_t0) # APACHE 2 score at t0
# Range 0 to 71: higher number = more severe condition
ggplot(wdf, aes(Gender, ls_diem_apache_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="APACHE 2 score at t0", y = "APACHE 2 score")
#16 is a possible value. Maintained


table(wdf$ls_diem_ranson_t0) # RANSON score at T0
# Range 0 to 8:
# A score of 3 or greater predicts severe acute pancreatitis and possible mortality
# If the score < 3, severe pancreatitis is unlikely
ggplot(wdf, aes(Gender, ls_diem_ranson_t0)) + geom_boxplot()+ stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title=" RANSON score at t0", y = " RANSON score")
#Possible values, outliers maintained.


table(wdf$ls_diem_ct_t0) # CTSI score at T0
# Range 0 to 10: higher number = higher pancreatic necrosis
ggplot(wdf, aes(Gender, ls_diem_ct_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="CTSI score at t0", y = "CTSI score")
#Possible values, outliers maintained.


table(wdf$ls_diem_imrie_t0) # IMRIE score at T0
# Range 0 to 8: higher number = more severe pancreatitis
ggplot(wdf, aes(Gender, ls_diem_imrie_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="IMRIE score at t0", y = "IMRIE score")
#Possible values, outliers maintained.


table(wdf$ls_diem_sofa_t0) # SOFA score at T0
# Range 0 to 24: higher number = higher organ failure assessment - higher chance of mortality
ggplot(wdf, aes(Gender, ls_diem_sofa_t0)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="SOFA score at t0", y = "SOFA score")
#Possible values, outliers maintained.


table(wdf$cls_ct_ctscore_lan1) # CTSI score with Computer Tomography
# Range 0 to 10: higher number = higher pancreatic necrosis
wdf$cls_ct_ctscore_lan1 <- as.numeric(wdf$cls_ct_ctscore_lan1)
ggplot(wdf, aes(Gender, cls_ct_ctscore_lan1)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="CTSI score with Computer Tomography", y = "CTSI score")
# Possible outlier at y = 23 and ID59 = e = 4
wdf$cls_ct_ctscore_lan1[59] <- 4
wdf$cls_ct_ctscore_lan1[132] <- NA

#Clean Boxplot
ggplot(wdf, aes(Gender, cls_ct_ctscore_lan1)) + geom_boxplot() + stat_summary(
  aes(label = round(stat(y), 1)), geom = "text", 
  fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
  hjust = -1
) + labs(title="CTSI score with Computer Tomography (Clean)", y = "CTSI score")
#Possible values, outliers maintained.

# cls_hh_aptt_t6 outlier at ID = 151 removed
wdf$cls_hh_aptt_t6[151] <- NA

# cls_sh_ka_tn6: tn6 meaning - Drop?
table(wdf$cls_sh_ka_tn6)
ggplot(wdf, aes(Gender, cls_sh_ka_tn6)) + geom_boxplot()


##### Naniar NA plots #####
### Subsetting ###
dfs <- list() # Defining master list of dataframes

## Demographics
dfs[['demo']] <- wdf[ , c(1:8, 18, 19, 188:190)]

## ts, ls, cls, dt
dfs[['ts']] <- mutate(wdf[ , grepl("ts_", names(wdf))], dfs$demo[c(1:3, 11, 13)]) # Patient Medical History

dfs[['ls']] <- mutate(wdf[ , grepl("\\<ls_", names(wdf))], dfs$demo[c(1:3, 11, 13)]) # Patient Clinical Symptoms

dfs[['cls']] <- mutate(wdf[ , grepl("cls_", names(wdf))], dfs$demo[c(1:3, 11, 13)]) # Patient Sub Clinical Examination

dfs[['dt']] <- mutate(wdf[ , grepl("dt_", names(wdf))], dfs$demo[c(1:3, 11, 13)]) # Patient Treatment

library(naniar)

gg_miss_var(dfs$demo, facet = pex, show_pct = TRUE) +
  labs(y = 'Missing [%]')

gg_miss_var(dfs$ts[ , -c((length(dfs$ts) - 4):(length(dfs$ts) - 1))], facet = pex, show_pct = TRUE) +
  labs(y = 'Missing [%]')

gg_miss_var(dfs$ls[ , -c((length(dfs$ls) - 4):(length(dfs$ls) - 1))], facet = pex, show_pct = TRUE) +
  labs(y = 'Missing [%]')

gg_miss_var(dfs$cls[ , -c((length(dfs$cls) - 4):(length(dfs$cls) - 1))], facet = pex, show_pct = TRUE) +
  labs(y = 'Missing [%]')

gg_miss_var(dfs$dt[ , -c((length(dfs$dt) - 4):(length(dfs$dt) - 1))], facet = pex, show_pct = TRUE) +
  labs(y = 'Missing [%]')

vis_miss(dfs$cls[ , grepl("_lan1", names(dfs$cls))], show_perc_col = FALSE) +
  coord_flip() +
  ggtitle('Missing values - Subclinical Examination | Computer Tomography') +
  facet_grid(dfs$cls$pex) +
  theme_bw()

vis_miss(dfs$cls[ , grepl("_t0", names(dfs$cls))], show_perc_col = FALSE) +
  coord_flip() +
  ggtitle('Missing values - Subclinical Examination t0') +
  facet_grid(dfs$cls$pex) +
  theme_bw()

vis_miss(dfs$cls[ , grepl("_t6", names(dfs$cls))], show_perc_col = FALSE) +
  coord_flip() +
  ggtitle('Missing values - Subclinical Examination t6') +
  facet_grid(dfs$cls$pex) +
  theme_bw()

vis_miss(dfs$cls[ , grepl("_t30", names(dfs$cls))], show_perc_col = FALSE) +
  coord_flip() +
  ggtitle('Missing values - Subclinical Examination t30') +
  facet_grid(dfs$cls$pex) +
  theme_bw()

vis_miss(dfs$cls[ , grepl("_t54", names(dfs$cls))], show_perc_col = FALSE) +
  coord_flip() +
  ggtitle('Missing values - Subclinical Examination t54') +
  facet_grid(dfs$cls$pex) +
  theme_bw()

vis_miss(dfs$cls[ , grepl("_t72", names(dfs$cls))], show_perc_col = FALSE) +
  coord_flip() +
  ggtitle('Missing values - Subclinical Examination t72') +
  facet_grid(dfs$cls$pex) +
  theme_bw()

## Row-wise NA visualization
gg_miss_fct(dfs$demo, fct = ID)
gg_miss_fct(dfs$ts, fct = ID)
gg_miss_fct(dfs$ls, fct = ID)
gg_miss_fct(dfs$cls, fct = ID)
gg_miss_fct(dfs$dt, fct = ID)

## Random Missing and Non at Random Missing
scatter_fun = function(x, y, x_name, y_name) {
  ggplot(dfs$dt,
         aes(x, y)) +
    geom_miss_point() +
    theme_bw() +
    xlab(x_name) +
    ylab(y_name)
}

scatter_fun(dfs$dt$dt_dich_vao_t24, dfs$dt$dt_dich_bilan_t24,
            'Fluid Intake t24', 'Fluid Balance t24')
scatter_fun(dfs$dt$dt_dich_ra_t24, dfs$dt$dt_dich_bilan_t24,
            'Fluid Output t24', 'Fluid Balance t24')

scatter_fun(dfs$dt$dt_dich_vao_t48, dfs$dt$dt_dich_bilan_t48,
            'Fluid Intake t48', 'Fluid Balance t48')
scatter_fun(dfs$dt$dt_dich_ra_t48, dfs$dt$dt_dich_bilan_t48,
            'Fluid Output t48', 'Fluid Balance t48')

scatter_fun(dfs$dt$dt_dich_vao_t72, dfs$dt$dt_dich_bilan_t72,
            'Fluid Intake t72', 'Fluid Balance t72')
scatter_fun(dfs$dt$dt_dich_ra_t72, dfs$dt$dt_dich_bilan_t72,
            'Fluid Output t72', 'Fluid Balance t72')

gg_miss_upset(dfs$dt[1:9], nsets = 9)


#### Box plots for Dt
coldt <- colnames(dfs$dt)
for(i in coldt){
  print(ggplot(wdf,aes(Gender, wdf[ ,i])) + geom_boxplot()+ stat_summary(
    aes(label = round(stat(y), 1)), geom = "text", 
    fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
    hjust = -1
  ) + labs(y =i))
}


#### Box plots for cls
colcls <- colnames(dfs$cls)
for(i in colcls){
  print(ggplot(wdf,aes(Gender, wdf[ ,i])) + geom_boxplot()+ stat_summary(
    aes(label = round(stat(y), 1)), geom = "text", 
    fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
    hjust = -1
  ) + labs(y =i))
}


#### Box plots for ls
colls <- colnames(dfs$ls)
for(i in colls){
  print(ggplot(wdf,aes(Gender, wdf[ ,i])) + geom_boxplot()+ stat_summary(
    aes(label = round(stat(y), 1)), geom = "text", 
    fun = function(y) { o <- boxplot.stats(y)$out; if(length(o) == 0) NA else o },
    hjust = -1
  ) + labs(y =i))
}

boxplot(select(dfs$dt, c(1:9)))
boxplot(select(dfs$dt, c(10:12)))
boxplot(select(dfs$dt, c(14:17)))
boxplot(select(dfs$dt, c(19:30)))


##### Dataset for Regression Model and MAR/MCAR imputation #####
## Variables to remove
rm_var1 <- c('dt_pex_ranson_s_lan1', # Only 2 observation
            'vv_Others', 'vv_reason1', 'vv_reason2', 'vv_reason3', # Categorical variables without meaningful information
            'details_ts_giadinh',
            'ts_vtc', 'ts_vtc_lancuoi',
            'ls_cn_bidaitien', 'ls_cn_ialong', 'ls_tht_bungchuong', 'ls_tt_lungsuon',
            'cls_sa_tuy_t0', 'cls_sa_mat_t0', 'cls_sa_ketluan_t0',
            'cls_ct_tuy_lan1', 'cls_ct_balthazar_lan1',
            'cls_sh_bil_t0', 'cls_sh_bil_t6', 'cls_sh_bil_t30', 'cls_sh_bil_t72',
            'cls_sh_gan_t0', 'cls_sh_gan_t6', 'cls_sh_gan_t30',
            'cls_sh_gan_t0', 'cls_sh_gan_t6', 'cls_sh_gan_t30',
            'cls_sh_ca_t0',
            'dead', 'ls_tn_cvp_t0', 'ls_tn_ha_t6') # Variables with inconsistent range of values

rm_var2 <- c('dt_pex_ngaybenh', "dt_pex_tri_t_lan1", "dt_pex_tri_s_lan1", "dt_pex_chol_t_lan1", "dt_pex_chol_s_lan1",
             "dt_pex_ldl_t_lan1", "dt_pex_ldl_s_lan1", "dt_pex_hdl_t_lan1", "dt_pex_apache_t_lan1", "dt_pex_apache_s_lan1",
             "dt_pex_ranson_t_lan1", "dt_pex_ranson_s_lan1", "dt_pex_imrie_t_lan1", "dt_pex_imrie_s_lan1",
             "dt_pex_balthazar_t_lan1", "dt_pex_balthazar_s_lan1", "dt_pex_sofa_t_lan1", "dt_pex_sofa_s_lan1",
             "dt_pex_alob_t_lan1", "dt_pex_alob_s_lan1") # Removing dt_pex_ variables with MNAR observations for standard treatment patients

### Dealing with MNAR (Missing Not At Random) observations ###
# Assigning 0 to t54 and t72 variables, for individuals that stayed less than 3 days in the hospital
x <- wdf$ID[wdf$rv_ngaydt < 3]
for (i in x) {
  wdf[wdf$ID == i, grepl("_t54", names(wdf))] <- 0
}

for (i in x) {
  wdf[wdf$ID == i, grepl("_t72", names(wdf))] <- 0
}

# Assigning 0 to patients that informed 'No' for drinking problem variable 'ts_ruou'
wdf$ts_ruou[is.na(wdf$ts_ruou)] <- 'No' # Assigning 'No' to missing values on 'ts_ruou'
wdf$ts_ruou_nam[wdf$ts_ruou == 'No'] <- 0 # Assigning 0 to every patient without declared drinking problem
wdf$ts_ruou_nam_ml[wdf$ts_ruou == 'No'] <- 0 # Assigning 0 to every patient without declared drinking problem

# Verifying variables triglycerids for PEX and cls clinical examinations to fill some of the MNAR
wdf$cls_sh_tri_t6
wdf$dt_pex_tri_s_lan1

# Assigning 'PEX Treatment' to patients that has values on variable dt_pex_lan
wdf$pex[!is.na(wdf$dt_pex_lan)] <- 'PEX Treatment'

# Fixing incorrect treatments for patients ID 23, 92, 31, 115
wdf$pex[wdf$ID == 23] <- 'Standard Treatment'
wdf$pex[wdf$ID == 92] <- 'Standard Treatment'
wdf$pex[wdf$ID == 31] <- 'PEX Treatment'
wdf$pex[wdf$ID == 115] <- 'PEX Treatment'

# Assigning 0 to variable dt_pex_lan, for patients that has 'Standard Treatment' on pex variable
wdf$dt_pex_lan[wdf$pex == 'Standard Treatment'] <- 0

x <- wdf$ID[is.na(wdf$cls_sh_tri_t6)]
for (i in x) {
  wdf$cls_sh_tri_t6[wdf$ID == i] <- wdf[wdf$ID == i, 'dt_pex_tri_s_lan1']
}

x <- wdf$ID[is.na(wdf$cls_sh_tri_t0)]
for (i in x) {
  wdf$cls_sh_tri_t0[wdf$ID == i] <- wdf[wdf$ID == i, 'dt_pex_tri_t_lan1']
}

fdf <- select(wdf, -c(rm_var1), -c(rm_var2))
names(fdf)
str(fdf)
table(fdf$cls_sh_ca_t0)
vis_miss(fdf)

pMiss <- function(x){round(sum(is.na(x))/length(x)*100, 2)}
perc_miss <- sapply(fdf, pMiss) # Identifying missing values percentage
perc_miss_75 <- names(perc_miss[perc_miss > 75]) # Vector of variables with missing values > 75%

# Final dataset to assign missing values
fdf <- select(fdf, -c(perc_miss_75))

##### Packages for Missing Values imputation #####
pacman::p_load(mice, Amelia, missForest)
vis_miss(fdf, show_perc_col = FALSE) +
  coord_flip() +
  facet_grid(fdf$pex) +
  theme_bw()
# 35.2% of missing values to enter mice, Amelia, and missForest packages

### MICE
cpus <- parallel::detectCores() # Verifying number of available CPUs
MiceRun <- parlmice(fdf[-1], m = 5, maxit = 5, cluster.seed = 123, core = cpus)
summary(MiceRun)
mice_df <- complete(MiceRun, 5)
mice_df$ID <- fdf$ID


### Amelia
# Subsets
ts_var <- grepl("ts_", names(fdf))
ls_var <- grepl("\\<ls_", names(fdf))
cls_var <- grepl("cls_", names(fdf))
dt_var <- grepl("dt_", names(fdf))

amelia_var <- function(df) {
  nominal_var <<- names(df[, sapply(df, is.factor)])
  ordinal_var <<- names(df[, sapply(df, is.integer)])
}

cpus <- parallel::detectCores() # Verifying number of available CPUs

# ts_ variables
amelia_var(fdf[ , ts_var])
str(fdf[ , ts_var])
AmeliaRun <- amelia(fdf[, ts_var], m = 5, p2s = 2, noms = nominal_var, #ords = ordinal_var,
                    seed = 123, parallel = 'snow', ncpus = cpus,
                    empri = 0.8*nrow(fdf[-1]))
amelia_ts <- AmeliaRun$imputations[[5]]

# ls_ variables
amelia_var(fdf[ , ls_var])
str(fdf[ , ls_var])
AmeliaRun <- amelia(fdf[, ls_var], m = 5, p2s = 2, ords = ordinal_var,
                    seed = 123, parallel = 'snow', ncpus = cpus,
                    empri = 0.8*nrow(fdf[-1]))
amelia_ls <- AmeliaRun$imputations[[5]] # Individuals 23 and 92 with no information

# cls_ variables
# cls_ subset 1
amelia_var(fdf[ , cls_var][1:44])
str(fdf[ , cls_var][1:44])
AmeliaRun <- amelia(fdf[, cls_var][1:44], m = 5, p2s = 2, noms = nominal_var, ords = ordinal_var,
                    seed = 123, parallel = 'snow', ncpus = cpus,
                    empri = 0.8*nrow(fdf[-1]))
amelia_cls1 <- AmeliaRun$imputations[[5]]

# cls_ subset 2
amelia_var(fdf[ , cls_var][45:88])
str(fdf[ , cls_var][45:88])
AmeliaRun <- amelia(fdf[, cls_var][45:88], m = 5, p2s = 2,
                    seed = 123, parallel = 'snow', ncpus = cpus,
                    empri = 0.8*nrow(fdf[-1]))
amelia_cls2 <- AmeliaRun$imputations[[5]]

# dt_ variables
amelia_var(fdf[ , dt_var])
str(fdf[ , dt_var])
AmeliaRun <- amelia(fdf[, dt_var], m = 5, p2s = 2, ords = c(ordinal_var, 'dt_pex_lan'),
                    seed = 123, parallel = 'snow', ncpus = cpus,
                    empri = 0.8*nrow(fdf[-1]))
amelia_dt <- AmeliaRun$imputations[[5]]

# Amelia final DF
amelia_df <- cbind(amelia_ts, amelia_ls, amelia_cls1, amelia_cls2, amelia_dt)
amelia_df[ , c('ID', 'Age', 'Gender', 'rv_ngaydt', 'stomachache', 'vomiting', 'complication', 'pex')] <- 
  fdf[ , c('ID', 'Age', 'Gender', 'rv_ngaydt', 'stomachache', 'vomiting', 'complication', 'pex')]
rm(amelia_cls1, amelia_cls2, amelia_dt, amelia_ls, amelia_ts)


### missForest
cpus <- parallel::detectCores() # Verifying number of available CPUs
doParallel::registerDoParallel(cores = cpus) # Set based on number of CPU cores
doRNG::registerDoRNG(seed = 123)
ForestRun <- missForest(fdf[-1], parallelize = 'variables', verbose = TRUE)
forest_df <- ForestRun$ximp
forest_df$ID <- fdf$ID


##### Imputations comparison #####
### cls_hh_bc variables
# Original DF
vis_miss(fdf[ , c('cls_hh_bc_t0', 'cls_hh_bc_t6', 'cls_hh_bc_t30', 'cls_hh_bc_t72')])

cls.melt <- melt(fdf, id.vars = 'ID',
                 measure.vars = c('cls_hh_bc_t0', 'cls_hh_bc_t6', 'cls_hh_bc_t30', 'cls_hh_bc_t54', 'cls_hh_bc_t72'))
ggplot(cls.melt) +
  geom_boxplot(aes(ID, value, color = variable)) +
  labs(x = '', y = 'Exam Results', title = 'Boxplot for cls_hh_bc variables') +
  theme_bw()

ggplot(cls.melt) +
  geom_histogram(aes(value)) +
  facet_grid(cls.melt$variable) +
  theme_bw()

# Amelia DF
amelia_df$ID <- fdf$ID
cls.melt2 <- melt(amelia_df, id.vars = 'ID',
                  measure.vars = c('cls_hh_bc_t0', 'cls_hh_bc_t6', 'cls_hh_bc_t30', 'cls_hh_bc_t54', 'cls_hh_bc_t72'))
ggplot(cls.melt2) +
  geom_boxplot(aes(ID, value, color = variable)) +
  labs(x = '', y = 'Exam Results', title = 'Boxplot for cls_hh_bc variables') +
  theme_bw()

ggplot(cls.melt2) +
  geom_histogram(aes(value)) +
  facet_grid(cls.melt2$variable) +
  theme_bw()

# MICE DF
mice_df$ID <- fdf$ID
cls.melt3 <- melt(mice_df, id.vars = 'ID',
                  measure.vars = c('cls_hh_bc_t0', 'cls_hh_bc_t6', 'cls_hh_bc_t30', 'cls_hh_bc_t54', 'cls_hh_bc_t72'))

ggplot(cls.melt3) +
  geom_boxplot(aes(ID, value, color = variable)) +
  labs(x = '', y = 'Exam Results', title = 'Boxplot for cls_hh_bc variables') +
  theme_bw()

ggplot(cls.melt3) +
  geom_histogram(aes(value)) +
  facet_grid(cls.melt3$variable) +
  theme_bw()

# missForest DF
forest_df$ID <- fdf$ID
cls.melt4 <- melt(forest_df, id.vars = 'ID',
                  measure.vars = c('cls_hh_bc_t0', 'cls_hh_bc_t6', 'cls_hh_bc_t30', 'cls_hh_bc_t54', 'cls_hh_bc_t72'))

ggplot(cls.melt4) +
  geom_boxplot(aes(ID, value, color = variable)) +
  labs(x = '', y = 'Exam Results', title = 'Boxplot for cls_hh_bc variables') +
  theme_bw()

ggplot(cls.melt4) +
  geom_histogram(aes(value)) +
  facet_grid(cls.melt4$variable) +
  theme_bw()

# Removing melt dfs
rm(cls.melt, cls.melt2, cls.melt3, cls.melt4)

### dt_ variables
# Original DF
boxplot(fdf[, c(111:119)])
vis_miss(fdf[, c(111:119)])
vis_miss(fdf[ , c("dt_dich_vao_t24", "dt_dich_vao_t48", "dt_dich_vao_t72", "dt_dich_ra_t24", "dt_dich_ra_t48",
                  "dt_dich_ra_t72", "dt_dich_bilan_t24", "dt_dich_bilan_t48", "dt_dich_bilan_t72")])

cls.melt <- melt(fdf, id.vars = 'ID',
                 measure.vars = c("dt_dich_vao_t24", "dt_dich_vao_t48", "dt_dich_vao_t72", "dt_dich_ra_t24", "dt_dich_ra_t48",
                                  "dt_dich_ra_t72", "dt_dich_bilan_t48", "dt_dich_bilan_t72"))

ggplot(cls.melt) +
  geom_boxplot(aes(ID, value, color = variable)) +
  labs(x = '', y = 'Exam Results', title = 'Boxplot for cls_km_ variables') +
  theme_bw()

ggplot(cls.melt) +
  geom_histogram(aes(value)) +
  facet_grid(cls.melt$variable, shrink = FALSE) +
  theme_bw()

# Amelia DF
cls.melt2 <- melt(amelia_df, id.vars = 'ID',
                 measure.vars = c("dt_dich_vao_t24", "dt_dich_vao_t48", "dt_dich_vao_t72", "dt_dich_ra_t24", "dt_dich_ra_t48",
                                  "dt_dich_ra_t72", "dt_dich_bilan_t48", "dt_dich_bilan_t72"))

ggplot(cls.melt2) +
  geom_boxplot(aes(ID, value, color = variable)) +
  labs(x = '', y = 'Exam Results', title = 'Boxplot for cls_km_ variables') +
  theme_bw()

ggplot(cls.melt2) +
  geom_histogram(aes(value)) +
  facet_grid(cls.melt$variable, shrink = FALSE) +
  theme_bw()

# MICE DF
cls.melt3 <- melt(mice_df, id.vars = 'ID',
                  measure.vars = c("dt_dich_vao_t24", "dt_dich_vao_t48", "dt_dich_vao_t72", "dt_dich_ra_t24", "dt_dich_ra_t48",
                                   "dt_dich_ra_t72", "dt_dich_bilan_t48", "dt_dich_bilan_t72"))

ggplot(cls.melt3) +
  geom_boxplot(aes(ID, value, color = variable)) +
  labs(x = '', y = 'Exam Results', title = 'Boxplot for cls_km_ variables') +
  theme_bw()

ggplot(cls.melt3) +
  geom_histogram(aes(value)) +
  facet_grid(cls.melt$variable, shrink = FALSE) +
  theme_bw()

# missForest DF
cls.melt4 <- melt(forest_df, id.vars = 'ID',
                  measure.vars = c("dt_dich_vao_t24", "dt_dich_vao_t48", "dt_dich_vao_t72", "dt_dich_ra_t24", "dt_dich_ra_t48",
                                   "dt_dich_ra_t72", "dt_dich_bilan_t48", "dt_dich_bilan_t72"))

ggplot(cls.melt4) +
  geom_boxplot(aes(ID, value, color = variable)) +
  labs(x = '', y = 'Exam Results', title = 'Boxplot for cls_km_ variables') +
  theme_bw()

ggplot(cls.melt4) +
  geom_histogram(aes(value)) +
  facet_grid(cls.melt$variable, shrink = FALSE) +
  theme_bw()

boxplot(fdf[, c(86:95)])
vis_miss(fdf[, c(86:95)])

##### Linear Regression Datasets #####
lm_dfs <- list()
lm_dfs[['Original']] <- fdf
lm_dfs[['MICE']] <- mice_df
lm_dfs[['Amelia']] <- amelia_df
lm_dfs[['missForest']] <- forest_df
rm(mice_df, amelia_df, forest_df)

# Checking for final datasets % pecent
vis_miss(lm_dfs$Original)
vis_miss(lm_dfs$MICE)
vis_miss(lm_dfs$Amelia) # Individuals with ID = 23, and 92
vis_miss(lm_dfs$missForest)

# Dropping individuals ID = 23, and 92
lm_dfs[['Original']] <- lm_dfs[['Original']][-c(23, 92), ]
lm_dfs[['MICE']] <- lm_dfs[['MICE']][-c(23, 92), ]
lm_dfs[['Amelia']] <- lm_dfs[['Amelia']][-c(23, 92), ]
lm_dfs[['missForest']] <- lm_dfs[['missForest']][-c(23, 92), ]

write.csv(lm_dfs[['Original']], file = 'data/original_df.csv', row.names = FALSE)
write.csv(lm_dfs[['MICE']], file = 'data/MICE_df.csv', row.names = FALSE)
write.csv(lm_dfs[['Amelia']], file = 'data/Amelia_df.csv', row.names = FALSE)
write.csv(lm_dfs[['missForest']], file = 'data/missForest_df.csv', row.names = FALSE)
