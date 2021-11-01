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
for (i in numeric) {
  print(group_by(missForest_df, pex) %>%
          summarise(
            count = n(),
            mean = mean(i, na.rm = TRUE),
            sd = sd(i, na.rm = TRUE)
          ))
}



###t-test to evaluate whether the means are different.
##Assumption 1: Are the Samples independent? Yes, the two treatments are independent of each other
##Assumption 2: Normal Distributed?
with(missForest_df, shapiro.test(Age[pex == "PEX Treatment"]))
with(missForest_df,shapiro.test(Age[pex == "Standard Treatment"]))
#p-value of 0.2843 and 0.01233, Fail to reject, data not normally distributed
hist(pext$Age)
hist(stdt$Age)

##Assumption 3: Same variance of the two population?
res.ftest <- var.test(Age ~ pex, data = missForest_df)
res.ftest
#The p-value 0.8367 is greater than the significance level alpha = 0.05.
#In conclusion, there is no significant difference between the variances of the two sets of data.
#Therefore, we can use the classic t-test witch assume equality of the two variances.

##T-test -  Is there any significant difference between the two groups?
Agettest <- t.test(Age ~ pex, data = missForest_df, var.equal = TRUE)
Agettest
#The p-value of the test is 0.02984, which is less than the significance level alpha = 0.05.
#We can conclude that pex patient's average age is significantly different from std treatment's average age with a p-value = 0.02984.