pacman::p_load(tidyverse, olsrr, forecast, corrr, caret, GGally,
               lmtest, car, rsample, class, lime, reshape2, ggpubr, usethis)

setwd("~/Langara/DANA 4830 - 001/Assignment 2/Pancreatitis")

missForest_df <- read.csv('data/final_models/missForest_df.csv', stringsAsFactors = TRUE)

summary(missForest_df)

#Complication was numeric, changed to Factor
table(missForest_df$complication)
missForest_df$complication <- as.factor(missForest_df$complication)

### Defining numeric variables only df ###
NumDF <- dplyr::select_if(missForest_df, is.numeric)
numeric <- names(missForest_df[ , sapply(missForest_df, is.numeric)])
### Removing ID column ###
NumWDF <- NumDF[,c(2:114)]

#### Master Table - Numeric ####
## Mean, Standard Deviation, Normality Check, Variance Check and T-test of all numeric variables.
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

### LDA  ####
library(MASS)
##Assumption 1:	Data is normally distributed - We checked in the shapiro-wilk test that data is not normally distributed
##Histograms of the final dataframe
Rnames <- rownames(Numericdf)

#Split the data into training (80%) and test set (20%)
LDAdf <- missForest_df[,c("ID",Rnames,categorical[-c(3,6)])]

set.seed(123)
training.samples <-LDAdf$pex %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- LDAdf[training.samples, ]
test.data <- LDAdf[-training.samples, ]

## Normalizing the data
# Estimate pre-processing parameters
preproc.param <- train.data[,-1] %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data[,-1])
test.transformed <- preproc.param %>% predict(test.data[,-1])

# Fit the model
ldamodel <- lda(pex~., data = train.transformed)
ldamodel
plot(ldamodel)

# Make predictions
predictions <- ldamodel %>% predict(test.transformed)
# Predicted probabilities of class membership.
predictions$posterior
# Linear discriminant - Shows the linear combination of predictor variables that are used to form the LDA decision rule
predictions$x

# Model accuracy - 96.87% chance of correctly classifying the patients in its groups.
mean(predictions$class==test.transformed$pex)

### Dataset with miss classification analysis ####
missForest_df2 <- missForest_df

# Reversing the qualitative manual adjust made in the Assignment2 for patients ID 23, 25, 31, 75, 92, 115, 121, and 122 
#ID 23 and 92 were completely blank, thus were removed
missForest_df2$pex[missForest_df2$ID == 25] <- 'Standard Treatment'
missForest_df2$pex[missForest_df2$ID == 31] <- 'Standard Treatment'
missForest_df2$pex[missForest_df2$ID == 75] <- 'Standard Treatment'
missForest_df2$pex[missForest_df2$ID == 115] <- 'Standard Treatment'
missForest_df2$pex[missForest_df2$ID == 121] <- 'Standard Treatment'
missForest_df2$pex[missForest_df2$ID == 122] <- 'Standard Treatment'

### LDA for the 2nd case with misclassified individuals ####
LDAdf2 <- missForest_df2[,c("ID",Rnames,categorical[-c(3,6)])]

missclass <- c(25, 31, 75, 115, 121, 122)
test.data2 <- filter(LDAdf2, ID %in% missclass)
`%notin%` <- Negate(`%in%`)
train.data2 <- filter(LDAdf2, ID %notin% missclass)

## Normalizing the data
# Estimate pre-processing parameters
preproc.param2 <- train.data2[,-1] %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed2 <- preproc.param %>% predict(train.data2[,-1])
test.transformed2 <- preproc.param %>% predict(test.data2[,-1])

# Fit the model
ldamodel2 <- lda(pex~., data = train.transformed2)
ldamodel2
plot(ldamodel2)

# Make predictions
predictions2 <- ldamodel2 %>% predict(test.transformed2)
# Predicted probabilities of class membership.
predictions2$posterior
# Linear discriminant - Shows the linear combination of predictor variables that are used to form the LDA decision rule
predictions2$x
# Model accuracy - 100% chance of correctly classifying the patients in its groups.
mean(predictions2$class==test.transformed2$pex)

### Removing complication variable
## Variable complication is heavily influencing the prediction of these 6 individuals. We are removing it on the following model
ldamodel3 <- lda(pex~., data = train.transformed2[, -40])
ldamodel3
plot(ldamodel3)
# Make predictions
predictions3 <- ldamodel3 %>% predict(test.transformed2[, -40])
# Predicted probabilities of class membership.
predictions3$posterior
# Linear discriminant - Shows the linear combination of predictor variables that are used to form the LDA decision rule
predictions3$x
# Model accuracy - 100% chance of correctly classifying the patients in its groups.
mean(predictions3$class==test.transformed2$pex)


#### LDA on the entire dataset ####
transformed4 <- preproc.param %>% predict(LDAdf2[,-1])
# Fit the model
ldamodel4 <- lda(pex~., data = transformed4[,-40])
ldamodel4
plot(ldamodel4)

# Make predictions
predictions4 <- ldamodel4 %>% predict(transformed4[,-40])
# Predicted probabilities of class membership.
round(predictions4$posterior, 3)
# Linear discriminant - Shows the linear combination of predictor variables that are used to form the LDA decision rule
predictions4$x
# Dataset accuracy according to model - 96.32%
mean(predictions4$class==transformed4$pex)

# Fitted dataset analysis according to model
fit.df <- data.frame(round(predictions4$posterior, 3))
fit.df <- cbind(LDAdf2$ID, LDAdf2$pex, fit.df)

# Identifying misclassified individuals
fit.df$misclassified <- ifelse(fit.df$PEX.Treatment > 0.5,
                               ifelse(fit.df$`LDAdf2$pex` == 'PEX Treatment', '', 'x'),
                               ifelse(fit.df$`LDAdf2$pex` == 'Standard Treatment', '', 'x'))
fit.df$ldapredicted <- ifelse(fit.df$PEX.Treatment > 0.5, 'PEX Treatment', 'Standard Treatment')
names(fit.df) <- c('ID', 'Original Treatment', 'LDA Pex %', 'LDA Standard %', 'Misclassified', 'LDA Predicted Treatment')
fit.df # 6 individuals misclassified according to model
lda.confusion.matrix <- table(fit.df[, c('Original Treatment', 'LDA Predicted Treatment')])
lda.confusion.matrix


#### KNN for Classification ####
set.seed(123)
knn.train <- LDAdf$pex %>%
  createDataPartition(p = 0.8, list = FALSE)
knn.train.data <- LDAdf[knn.train, -c(1, 41)]
knn.test.data <- LDAdf[-knn.train, -c(1, 41)]

knn.model1 <- train(pex~., data = knn.train.data, method = 'knn', # knn method for classification
                    trControl = trainControl('cv', number = 20), # 10-fold cross validation
                    preProcess = c('center', 'scale'), # Normalizing and centering the data
                    tuneLength = 50)
plot(knn.model1)
knn.model1$results
knn.model1$bestTune # Accuracy peak at k=31 neighbors
knn.predictions1 <- knn.model1 %>% predict(knn.test.data)
knn.predictions1
# Prediction accuracy - 90.63% for test dataset
mean(knn.predictions1==knn.test.data$pex)

# Classification for the entire dataset
knn.predictions2 <- knn.model1 %>% predict(LDAdf[, -c(1, 41)]) # Removing ID and complication variables
# Prediction accuracy - 91.41% for the entire dataset
mean(knn.predictions2==LDAdf$pex)

# Fitted dataset analysis according to model
fit.knn.df <- data.frame(LDAdf2$ID, LDAdf2$pex, knn.predictions2)
names(fit.knn.df) <- c('ID', 'Original Treatment', 'KNN Predicted Treatment')
fit.knn.df

# Identifying misclassified individuals
fit.knn.df$misclassified <- ifelse(fit.knn.df$`Original Treatment` == fit.knn.df$`KNN Predicted Treatment`, '', 'x')
fit.knn.df # 20 individuals misclassified according to KNN model
# Misclassified individuals according to KNN model:
# 25, 26, 31, 46, 49, 51, 61, 62, 65, 72, 75, 88, 96, 107, 112, 115, 118, 121, 122, and 145
knn.confusion.matrix <- table(fit.knn.df[, c('Original Treatment', 'KNN Predicted Treatment')])
knn.confusion.matrix

