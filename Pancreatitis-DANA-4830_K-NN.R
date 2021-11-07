pacman::p_load(tidyverse, olsrr, forecast, corrr, caret, GGally,
               lmtest, car, rsample, class, lime, reshape2, ggpubr, usethis)

setwd("~/Langara/DANA 4830 - 001/Assignment 2/Pancreatitis")

missForest_df <- read.csv('~/Langara/DANA 4830 - 001/Assignment 2/Pancreatitis/data/final_models/missForest_df.csv', stringsAsFactors = TRUE)

summary(missForest_df)

#Complication, was numeric, changed to Factor
table(missForest_df$complication)
missForest_df$complication <- as.factor(missForest_df$complication)

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

#### Reducing Number of Variables for T-test
### Defining number of Factors ###
nofactors = fa.parallel(NumWDF, fm="ml", fa="fa")
nofactors$fa.values #13 Factors recommended
sum(nofactors$fa.values > 1.0) # old kaiser criterion, 16 Factors
sum(nofactors$fa.values > .7) # new kaiser criterion, 18 Factors

### FA Models ####
## 13 Factors test ##
fit.one <- factanal(NumWDF,factors=13,rotation="oblimin")
print(fit.one, digits = 2, cutoff = .2, sort = TRUE) #reject hypothesis that 13 factors are perfect fit, p-value 3.47e^-124
#The null hypothesis, H0, is that the number of factors in the model, is sufficient to capture the full dimensionality
#of the data set. 
#We reject H0 if the p-value is less than 0.05.Such a result indicates that the number of factors is too small.
#If we do not reject H0,the result indicates that there are likely enough (or more than enough) factors capture the full dimensionality
#of the data set (Teetor 2011).
#When you accept the H0 you have an appropriate model.

model.one <- fa(NumWDF, 
                nfactors=13, 
                rotate = "oblimin", 
                fm = "ml")
fa.diagram(model.one)
model.one$Phi # Rotation
print(model.one, digits = 2, cutoff = .2, sort = TRUE) #13 Factors describes 47% of the variance


## 16 Factors test ##
fit.two <- factanal(NumWDF,factors=16,rotation="oblimin")
print(fit.two, digits = 2, cutoff = .2, sort = TRUE) #reject hypothesis that 16 factors are perfect fit, p-value 5.89e^-86

model.two <- fa(NumWDF, 
                nfactors=16, 
                rotate = "oblimin", 
                fm = "ml")
fa.diagram(model.two)
model.two$Phi # Rotation
print(model.two, digits = 2, cutoff = .2, sort = TRUE) #16 Factors describes 52% of the variance

## 18 Factors test ##
fit.three <- factanal(NumWDF,factors=18,rotation="oblimin")
print(fit.three, digits = 2, cutoff = .2, sort = TRUE) #reject hypothesis that 18 factors are perfect fit, p-value 9.08e^-68

model.three <- fa(NumWDF, 
                  nfactors=18, 
                  rotate = "oblimin", 
                  fm = "ml")
fa.diagram(model.three)
model.three$Phi # Rotation
print(model.three, digits = 2, cutoff = .2, sort = TRUE) #18 Factors describes 55% of the variance
#The 18 Factor model describes more variance of the Dataframe, but it also increase its complexity, thus
#the group decided to use the 13 Factor model.

### Refining the Factor Analysis model ####
uniq <- data.frame(model.one$uniquenesses)
uniq <- rownames_to_column(uniq, var = "Variable")
uniq[uniq[,2] >0.5, "Variable"]

### FACTOR ANALYSIS SCORE T-TEST ###
m1df <- data.frame(model.one$scores)
m1df$pex <- paste(missForest_df$pex)
t.test(ML1 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML2 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML3 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML4 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML5 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML6 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML7 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML8 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML9 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML10 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML11 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML12 ~ pex, data = m1df, var.equal = TRUE)
t.test(ML13 ~ pex, data = m1df, var.equal = TRUE)



####For loop t-test all numeric variables ####
##Assumption 1: Are the Samples independent? Yes, the two treatments are independent of each other 
##Assumption 2: Normal Distributed?
for (i in numeric) {
  print(i)
  print(with(missForest_df, shapiro.test(missForest_df[,i][pex == "PEX Treatment"])))
  hist(pext[,i],main=i)
}

for (i in numeric) {
  print(i)
  print(with(missForest_df,shapiro.test(missForest_df[,i][pex == "Standard Treatment"])))
  hist(stdt[,i],main=i)
}

##Assumption 3: Same variance of the two population?
for (i in numeric) {
  print(i)
  print(var.test(missForest_df[,i] ~ pex, data = missForest_df))
}


#### Master Table - Numeric ####
##Mean, Standard Deviation, Normality Check, Variance Check and T-test of all numeric variables.
tpvalues <- list()
normpvalues <- list()
varipvalues <- list()
pextmean <- list()
stdtmean <- list()
pextstd <- list()
stdtstd <- list()

for (i in numeric[2:113]) {#Numeric ID 2: 113 to remove the ID and X variables
  t <- t.test(missForest_df[,i] ~ pex, data = missForest_df, var.equal = TRUE)
  norm <- with(missForest_df, shapiro.test(missForest_df[,i]))
  varia <- var.test(missForest_df[,i] ~ pex, data = missForest_df)
  
  gmean <-group_by(missForest_df,pex) %>%
    summarise(
      disp_mean = mean(!!sym(i), na.rm = TRUE),
      disp_sd = sd(!!sym(i), na.rm = TRUE)
    )
  
  tpvalues[[i]] <- round(t$p.value,4)
  normpvalues[[i]] <- round(norm$p.value,4)
  varipvalues[[i]] <- round(varia$p.value,4)
  pextmean[[i]] <- round(gmean[1,2],2)
  stdtmean[[i]] <- round(gmean[2,2],2)
  pextstd[[i]] <- round(gmean[1,3],2)
  stdtstd[[i]] <- round(gmean[2,3],2)
}

tpvalues <- do.call("rbind", tpvalues)
normpvalues <- do.call("rbind", normpvalues)
varipvalues <- do.call("rbind", varipvalues)
pextmean <- do.call("rbind", pextmean)
stdtmean <- do.call("rbind", stdtmean)
pextstd <- do.call("rbind", pextstd)
stdtstd <- do.call("rbind", stdtstd)

Numericdf <- data.frame(pextmean,stdtmean,pextstd,stdtstd,normpvalues,varipvalues,tpvalues)
names(Numericdf) <- c("Pex T Mean","Standard T Mean","Pex T Stdev","Standard T Stdev","Shapiro-W p-values","Var-test p-values",
               "T-test p-values")
#Filtering all variables with t-test p-value <0.05
Numericdf <- Numericdf[Numericdf[,7] <0.05,]

#### Categorical ####
categorical <- names(missForest_df[ , sapply(missForest_df, is.factor)])
#Removing categorical variables ts_benhmat and stomachache, since they do not produce meaningful summary tables
chisqpvalues <- list()

for (i in categorical[-c(3,6)]) {
  print(i)
  chisqt <- chisq.test(missForest_df$pex, missForest_df[,i])
  print(chisqt$observed)
  print(chisqt$expected)
  print(chisqt$residuals)
  
  chisqpvalues[[i]] <-round(chisqt$p.value,4)
}
chisqpvalues <- do.call("rbind", chisqpvalues)

### LDA ####
##Assumption 1:	Data is normally distributed - We checked in the shapiro-wilk test that data is not normally distributed
##Histograms of the final dataframe
Rnames <- rownames(Numericdf)
for (i in Rnames) {
  hist(missForest_df[,i],main=i)
}

# Split the data into training (80%) and test set (20%)
LDAdf <- missForest_df[,c(Rnames,categorical)]

set.seed(123)
training.samples <- missForest_df$Species %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- iris[training.samples, ]
test.data <- iris[-training.samples, ]



#### t-test to evaluate whether the means are different on the SCORES variables####
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
