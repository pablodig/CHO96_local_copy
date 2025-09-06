library(data.table)
library(tidyverse)
library(dplyr)
library(caret)
library(pls)
library(ggplot2)
library(glmnet)
library(ggpubr)
library(caretEnsemble)
library(RhpcBLASctl)
library(Matrix)
library(nnet)
library(kernlab)
library(gbm)
library(plyr)
library(mboost)
library(randomForest)
library(Cubist)


###### limit number of cores used in jabba
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(35)


### use cleaned up Tegel data from correlation analysis
allFeatures <- readRDS("Data/Tegel2020_allFeatures.RDS")
featureMatrix <- allFeatures$features_passSamples
proteins.dt_pass <- allFeatures$proteins.dt_pass
identical(featureMatrix$Cultivation, proteins.dt_pass$Cultivation)


# extract id columns
featureMatrix.id <- featureMatrix[,1:5] 
featureMatrix <- featureMatrix[, -c(1:5)]


### Some model methods can't handle 'complex' features names
# Renames columns of feature matrix to generic names (e.g. Feature_1, Feature_2, etc.)
featureName.map <- data.frame("feature.name"=c("amount.µg", colnames(featureMatrix)))
colnames(featureMatrix) <- paste("Feature_", seq(1:(ncol(featureMatrix))), sep="")
featureMatrix <- cbind("response.var"=log10(proteins.dt_pass$Protein.amount..µg.), featureMatrix)
featureName.map$feature.ID <- colnames(featureMatrix)

preProcess.MLdata <- function(ML.data) {
  ## Remove near-zero variance predictors
  nzv <- nearZeroVar(ML.data, saveMetrics = TRUE) #nzv metrics for all variables
  #View(nzv[nearZeroVar(ML.data),]) #predictors with near zero variance
  filteredVars <- as.data.frame(ML.data)[,-nearZeroVar(ML.data)]
  
  
  ## Identify incomplete features (features w/ NAs)
  NAcols_high <- colnames(filteredVars)[colSums(is.na(filteredVars)) > 0.1*nrow(filteredVars)] #(5-10% does not make a difference)
  NAcols_all <- colnames(filteredVars)[colSums(is.na(filteredVars)) > 0]
  rm.highNA <- filteredVars[, !(colnames(filteredVars) %in% NAcols_high)]
  rm.allNA <- filteredVars[, !(colnames(filteredVars) %in% NAcols_all)]
  
  ## Apply pre-processing to data
  ppData <- list()
  ppData[['knnImp']]<- preProcess(rm.highNA, method=c("knnImpute", "center", "scale")) %>% predict(rm.highNA) #, "YeoJohnson"
  ppData[['rm.allNA']] <- preProcess(rm.allNA, method=c("center", "scale")) %>% predict(rm.allNA) #, "YeoJohnson"
  
  return(ppData)
}

ppData <- preProcess.MLdata(featureMatrix)

## Remove highly correlated variables, make sure to keep MW
MW_ID <- featureName.map$feature.ID[featureName.map$feature.name=="sequence_Mol.Weight_MW..Da."]
rm.highCor <- function(data, cutoff) {
  features <- subset(data, select=-c(response.var))
  corMatrix <- cor(features)
  # want to keep MW
  corMF <- names(corMatrix[, MW_ID])[abs(corMatrix[, MW_ID]) > cutoff & corMatrix[, MW_ID] !=1]
  features <- features[, !(names(features) %in% corMF)] #drop columns that are highly correlated with MW
  # repeat
  corMatrix <- cor(features)
  highlyCor <- findCorrelation(corMatrix, cutoff = cutoff, exact=TRUE)
  data.rmCor <- cbind("response.var"=data$response.var, features[, -highlyCor])
  return(data.rmCor)
}

ppData <- lapply(ppData, rm.highCor, cutoff=.98)

generateModels <- function(chooseData) {
  ### NO data splitting.  Want to use ALL the data to develop a descriptive model (not predictive)
  train.data  <- chooseData
  
  ### Specify type of resampling: repeated 10-fold CV
  fitControl <- trainControl(
    method="repeatedcv",
    number=10, # 10-fold
    repeats=5, # repeated 5 times
    index=createMultiFolds(train.data$response.var, k = 10, times = 5), # 10-fold cv repeated 10 times
    savePredictions=T,
    classProbs = F
  )
  
  ### Generate series of models
  model_list = list()
  
  set.seed(123)
  model_list <- tryCatch(c(model_list,
                           list(elasticNet = caret::train(x=subset(train.data, select=-c(response.var)), y=train.data$response.var,
                                                          method = "glmnet",
                                                          trControl = fitControl,
                                                          tuneLength=30))),
                         error = function(cond){
                           message('elasticNet failed')
                           message(cond)
                           return(model_list)
                         })
  set.seed(123)
  model_list <- tryCatch(c(model_list,
                           list(PLS = caret::train(x=subset(train.data, select=-c(response.var)), y=train.data$response.var,
                                                   method = "pls",
                                                   trControl = fitControl,
                                                   tuneLength=30))),
                         error = function(cond){
                           message('PLS failed')
                           message(cond)
                           return(model_list)
                         })
  set.seed(123)
  model_list <- tryCatch(c(model_list,
                           list(avNNet = caret::train(x=subset(train.data, select=-c(response.var)), y=train.data$response.var,
                                                      method = "avNNet",
                                                      trControl = fitControl,
                                                      tuneLength=10))),
                         error = function(cond){
                           message('avNNet failed')
                           message(cond)
                           return(model_list)
                         })
  set.seed(123)
  model_list <- tryCatch(c(model_list,
                           list(svmRadial = caret::train(x=subset(train.data, select=-c(response.var)), y=train.data$response.var,
                                                         method = "svmRadial",
                                                         trControl = fitControl,
                                                         tuneLength=60))),
                         error = function(cond){
                           message('svmRadial failed')
                           message(cond)
                           return(model_list)
                         })
  set.seed(123)
  model_list <- tryCatch(c(model_list,
                           list(gbm = caret::train(x=subset(train.data, select=-c(response.var)), y=train.data$response.var,
                                                   method = "gbm",
                                                   trControl = fitControl,
                                                   tuneLength=6,
                                                   verbose = FALSE))),
                         error = function(cond){
                           message('gbm failed')
                           message(cond)
                           return(model_list)
                         })
  set.seed(123)
  model_list <- tryCatch(c(model_list,
                           list(glmboost = caret::train(x=subset(train.data, select=-c(response.var)), y=train.data$response.var,
                                                        method = "glmboost",
                                                        trControl = fitControl,
                                                        tuneLength=30))),
                         error = function(cond){
                           message('glmboost failed')
                           message(cond)
                           return(model_list)
                         })
  set.seed(123)
  model_list <- tryCatch(c(model_list,
                           list(randomForest = caret::train(x=subset(train.data, select=-c(response.var)), y=train.data$response.var,
                                                      method = "rf",
                                                      trControl = fitControl,
                                                      tuneLength=30))),
                         error = function(cond){
                           message('random forest failed')
                           message(cond)
                           return(model_list)
                         })
  set.seed(123)
  model_list <- tryCatch(c(model_list,
                           list(cubist = caret::train(x=subset(train.data, select=-c(response.var)), y=train.data$response.var,
                                                    method = "cubist",
                                                    trControl = fitControl,
                                                    tuneLength=15))),
                         error = function(cond){
                           message('cubist failed')
                           message(cond)
                           return(model_list)
                         })
  
  
  
  return(list("modelList"=model_list, "train.data"=train.data))
}

model_list <- lapply(ppData, generateModels)
model_list[["featureName.map"]] <- featureName.map #add feature name map to saved object
saveRDS(model_list, "Output/modelList_regression_descriptive_microgram.RDS")
