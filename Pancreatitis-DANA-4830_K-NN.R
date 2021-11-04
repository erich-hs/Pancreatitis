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
#For loop testing, to compute all columns
library(dplyr)
library(rlang)
numeric <- names(missForest_df[ , sapply(missForest_df, is.numeric)])
for (i in numeric) {
  print(i)
  print(group_by(missForest_df,pex) %>%
          summarise(
            count = n(),
            disp_mean = mean(!!sym(i), na.rm = TRUE),
            disp_sd = sd(!!sym(i), na.rm = TRUE)
          ))
}

#### Factor Analysis ####
library(psych)
library(GPArotation)

### Defining numeric variables only df ###
NumDF <- dplyr::select_if(missForest_df, is.numeric)

### Removing ID column ###
NumWDF <- NumDF[,c(2:114)]

### Defining number of Factors ###
nofactors = fa.parallel(NumWDF, fm="ml", fa="fa")
nofactors$fa.values #13 Factors recommended

sum(nofactors$fa.values > 1.0) # old kaiser criterion, 16 Factors
sum(nofactors$fa.values > .7) # new kaiser criterion, 18 Factors

### Models ###
## 13 Variables test ##
fit.one <- factanal(NumWDF,factors=13,rotation="oblimin")
print(fit.one)
print(fit.one, digits = 2, cutoff = .2, sort = TRUE)

model.one <- fa(NumWDF, 
                nfactors=13, 
                rotate = "oblimin", 
                fm = "ml")
fa.diagram(model.one)
model.one$Phi # Rotation
print(model.one, digits = 2, cutoff = .5, sort = TRUE) #13 Factors describes 47% of the variance

## 16 Variables test ##
fit.two <- factanal(NumWDF,factors=16,rotation="oblimin")
print(fit.two)
print(fit.two, digits = 2, cutoff = .2, sort = TRUE)

model.two <- fa(NumWDF, 
             nfactors=16, 
             rotate = "oblimin", 
             fm = "ml")
fa.diagram(model.two)
model.two$Phi # Rotation
print(model.two, digits = 2, cutoff = .5, sort = TRUE) #16 Factors describes 52% of the variance


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

#### Correlation tests
### Correlation of Blood Count Group ####
bloodcount <- missForest_df[,c(26:44)]
cor(bloodcount, method = "pearson")

### Correlation of Blood Gases Group ####
bloodgas <- missForest_df[,c(86:110)]
cor(bloodgas, method = "pearson")

### Correlation of Clynical Symptoms ####
clynsymp <- missForest_df[,c(13:17)]
cor(clynsymp, method = "pearson")

### Correlation of Coagulation Exams ####
coagexams <- missForest_df[,c(45:56)]
cor(coagexams, method = "pearson")

### Correlation of Disease Severity variables####
severity <- missForest_df[,c(18:22,25,66)]
cor(severity, method = "pearson")

### Correlation between CTSI Variables ###
CTSIcor <- cor(missForest_df$ls_diem_ct_t0,missForest_df$cls_ct_ctscore_lan1,method = "pearson")
CTSIcor

ggscatter(missForest_df, x = "ls_diem_ct_t0", y = "cls_ct_ctscore_lan1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "CTSI Score", ylab = "CTSI Tomography")

### Correlation of Inflamation Exams ####
inflaexams <- missForest_df[,c(57:65)]
cor(inflaexams, method = "pearson")

### Correlation of Pancreas Monitoring Exams ####
pancreasmonit <- missForest_df[,c(67:85)]
cor(pancreasmonit, method = "pearson")

### Correlation of Pex treatment exams ####
pexexams <- missForest_df[,c(121,122)]
cor(pexexams, method = "pearson")

### Correlation of STD treatment exams ####
stdexams <- missForest_df[,c(111:120)]
cor(stdexams, method = "pearson")
