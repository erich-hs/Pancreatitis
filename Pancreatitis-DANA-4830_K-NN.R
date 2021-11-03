pacman::p_load(tidyverse, olsrr, forecast, corrr, caret, GGally,
               lmtest, car, rsample, class, lime, reshape2, ggpubr, usethis)

setwd("~/Langara/DANA 4830 - 001/Assignment 2/Pancreatitis")

missForest_df <- read.csv('~/Langara/DANA 4830 - 001/Assignment 2/Pancreatitis/data/final_models/missForest_df.csv', stringsAsFactors = TRUE)

summary(missForest_df)

#Sub-setting the df by treatment
pext <- subset(missForest_df,missForest_df$pex == 'PEX Treatment')
stdt <- subset(missForest_df,missForest_df$pex == 'Standard Treatment')

#Computing summary statistics for each treatment group
summary(pext)
summary(stdt)

#Another way to compute summary statistics
library(dplyr)
group_by(missForest_df, pex) %>%
  summarise(
    count = n(),
    mean = mean(Age, na.rm = TRUE),
    sd = sd(Age, na.rm = TRUE)
)

#For loop testing, to compute all columns
numeric <- names(missForest_df[ , sapply(missForest_df, is.numeric)])
for (i in numeric) {
  print(i)
  print(group_by(missForest_df, pex) %>%
          summarise(
            count = n(),
            disp_mean = mean(missForest_df[[i]], na.rm = TRUE),
            disp_sd = sd(missForest_df[[i]], na.rm = TRUE)
          ))
}

#################### TEST
numeric <- names(missForest_df[ , sapply(missForest_df, is.numeric)])
for (i in numeric) {
  print(i)
  print(missForest_df %>% group_by(pex) %>%
          summarise(
            count = n(),
            disp_mean = mean(missForest_df[[i]], na.rm = TRUE),
            disp_sd = sd(missForest_df[[i]], na.rm = TRUE)
          ))
}


#### t-test to evaluate whether the means are different on the SCORES variables
### APACHE Score ####
##Assumption 1: Are the Samples independent? Yes, the two treatments are independent of each other 
##Assumption 2: Normal Distributed? 
with(missForest_df, shapiro.test(ls_diem_apache_t0[pex == "PEX Treatment"]))
with(missForest_df,shapiro.test(ls_diem_apache_t0[pex == "Standard Treatment"]))
#Pex Treatment p-value < 0.05, reject H0, data NOT normally distributed
#Std Treatment p-value < 0.05, reject H0, data NOT normally distributed
hist(pext$ls_diem_apache_t0)
hist(stdt$ls_diem_apache_t0)

##Assumption 3: Same variance of the two population?
res.ftest <- var.test(ls_diem_apache_t0 ~ pex, data = missForest_df)
res.ftest
#The p-value 0.2432 is greater than the significance level alpha = 0.05.
#In conclusion, there is no significant difference between the variances of the two sets of data.
#Therefore, we can use the classic t-test, ASSUMING NORMALITY, witch assume equality of the two variances.

##T-test -  Is there any significant difference between the two groups?
APACHEtest <- t.test(ls_diem_apache_t0 ~ pex, data = missForest_df, var.equal = TRUE)
APACHEtest
#The p-value of the test is 0.005851, which is less than the significance level alpha = 0.05. Reject Null Hypothesis. Means are different.
#We can conclude that pex patient's average APACHE Score is significantly different from std treatment's average APACHE Score with a p-value = 0.005851.

### RANSON Score ####
##Assumption 1: Are the Samples independent? Yes, the two treatments are independent of each other
##Assumption 2: Normal Distributed?
with(missForest_df, shapiro.test(ls_diem_ranson_t0[pex == "PEX Treatment"]))
with(missForest_df,shapiro.test(ls_diem_ranson_t0[pex == "Standard Treatment"]))
#Pex Treatment p-value < 0.05, reject H0, data NOT normally distributed
#Std Treatment p-value < 0.05, reject H0, data NOT normally distributed
hist(pext$ls_diem_ranson_t0)
hist(stdt$ls_diem_ranson_t0)

##Assumption 3: Same variance of the two population?
res.ftest <- var.test(ls_diem_ranson_t0 ~ pex, data = missForest_df)
res.ftest
#The p-value 0.3125 is greater than the significance level alpha = 0.05.
#In conclusion, there is no significant difference between the variances of the two sets of data.
#Therefore, we can use the classic t-test, ASSUMING NORMALITY, witch assume equality of the two variances.

##T-test -  Is there any significant difference between the two groups?
RANSONtest <- t.test(ls_diem_ranson_t0 ~ pex, data = missForest_df, var.equal = TRUE)
RANSONtest

#The p-value of the test is 0.05065, which is slighlty greater than the significance level alpha = 0.05. However, Reject Null Hypothesis. Means are slightly different.
#We can conclude that pex patient's average Ranson Score is slightly statistically significant different from std treatment's average Ranson Score
#with a p-value = 0.05065.



### CTSI Score ####
##Assumption 1: Are the Samples independent? Yes, the two treatments are independent of each other
##Assumption 2: Normal Distributed?
with(missForest_df, shapiro.test(ls_diem_ct_t0[pex == "PEX Treatment"]))
with(missForest_df,shapiro.test(ls_diem_ct_t0[pex == "Standard Treatment"]))
#Pex Treatment p-value < 0.05, reject H0, data NOT normally distributed
#Std Treatment p-value < 0.05, reject H0, data NOT normally distributed
hist(pext$ls_diem_ct_t0)
hist(stdt$ls_diem_ct_t0)

##Assumption 3: Same variance of the two population?
res.ftest <- var.test(ls_diem_ct_t0 ~ pex, data = missForest_df)
res.ftest
#The p-value 0.8652 is greater than the significance level alpha = 0.05.
#In conclusion, there is no significant difference between the variances of the two sets of data.
#Therefore, we can use the classic t-test, ASSUMING NORMALITY, witch assume equality of the two variances.

##T-test -  Is there any significant difference between the two groups?
CTSItest <- t.test(ls_diem_ct_t0 ~ pex, data = missForest_df, var.equal = TRUE)
CTSItest

#The p-value of the test is 0.01543, which is less than the significance level alpha = 0.05.Reject Null Hypothesis. Means are different.
#We can conclude that pex patient's average CTSI Score is significantly different from std treatment's average CTSI Score
#with a p-value = 0.01543.


### Imrie Score ####
##Assumption 1: Are the Samples independent? Yes, the two treatments are independent of each other
##Assumption 2: Normal Distributed?
with(missForest_df, shapiro.test(ls_diem_imrie_t0[pex == "PEX Treatment"]))
with(missForest_df,shapiro.test(ls_diem_imrie_t0[pex == "Standard Treatment"]))
#Pex Treatment p-value < 0.05, reject H0, data NOT normally distributed
#Std Treatment p-value < 0.05, reject H0, data NOT normally distributed
hist(pext$ls_diem_imrie_t0)
hist(stdt$ls_diem_imrie_t0)

##Assumption 3: Same variance of the two population?
res.ftest <- var.test(ls_diem_imrie_t0 ~ pex, data = missForest_df)
res.ftest
#The p-value 0.1226 is greater than the significance level alpha = 0.05.
#In conclusion, there is no significant difference between the variances of the two sets of data.
#Therefore, we can use the classic t-test, ASSUMING NORMALITY, witch assume equality of the two variances.

##T-test -  Is there any significant difference between the two groups?
Imrietest <- t.test(ls_diem_imrie_t0 ~ pex, data = missForest_df, var.equal = TRUE)
Imrietest

#The p-value of the test is 0.09801, which is greater than the significance level alpha = 0.05.Accept Null Hypothesis. Means are not different.
#We can conclude that pex patient's average Imrie Score is not significantly different from std treatment's average Imrie Score
#with a p-value = 0.09801.


### SOFA Score ####
##Assumption 1: Are the Samples independent? Yes, the two treatments are independent of each other
##Assumption 2: Normal Distributed?
with(missForest_df, shapiro.test(ls_diem_sofa_t0[pex == "PEX Treatment"]))
with(missForest_df,shapiro.test(ls_diem_sofa_t0[pex == "Standard Treatment"]))
#Pex Treatment p-value < 0.05, reject H0, data NOT normally distributed
#Std Treatment p-value < 0.05, reject H0, data NOT normally distributed
hist(pext$ls_diem_sofa_t0)
hist(stdt$ls_diem_sofa_t0)

##Assumption 3: Same variance of the two population?
res.ftest <- var.test(ls_diem_sofa_t0 ~ pex, data = missForest_df)
res.ftest
#The p-value 0.9603 is greater than the significance level alpha = 0.05.
#In conclusion, there is no significant difference between the variances of the two sets of data.
#Therefore, we can use the classic t-test, ASSUMING NORMALITY, witch assume equality of the two variances.

##T-test -  Is there any significant difference between the two groups?
SOFAtest <- t.test(ls_diem_sofa_t0 ~ pex, data = missForest_df, var.equal = TRUE)
SOFAtest

#The p-value of the test is 0.6958, which is greater than the significance level alpha = 0.05.Accept Null Hypothesis. Means are not different.
#We can conclude that pex patient's average SOFA Score is not significantly different from std treatment's average SOFA Score
#with a p-value = 0.6958.