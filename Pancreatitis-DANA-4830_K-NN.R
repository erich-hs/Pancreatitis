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
