library(data.table)
library(tidyverse)
library(dplyr)
library(caret)
library(pls)
library(ggplot2)
library(glmnet)
library(Matrix)
library(ggpubr)
library(caretEnsemble)
library(klaR)
library(nnet)
library(kernlab)
library(gbm)
library(plyr)
library(mboost)
library(randomForest)


###### limit number of cores used in jabba
library(RhpcBLASctl)
RhpcBLASctl::blas_set_num_threads(15)


### use cleaned up Tegel data from correlation analysis
allFeatures <- readRDS("Data/Tegel2020_allFeatures.RDS")
featureMatrix <- allFeatures$features_fail0
proteins.dt_multiClass <- allFeatures$proteins.dt_fail0
identical(featureMatrix$Cultivation, proteins.dt_multiClass$Cultivation)

### simple 2 class classification: pass vs fail
proteins.dt_2Class <- proteins.dt_multiClass[grepl("Fail", proteins.dt_multiClass$Status), Status := "Fail"] #set every kind of fail to simple "Fail"
identical(featureMatrix$Cultivation, proteins.dt_2Class$Cultivation)

# extract id columns
featureMatrix.id <- featureMatrix[,1:5] 
featureMatrix <- featureMatrix[, -c(1:5)]

### Some model methods can't handle 'complex' features names
# Rename columns of feature matrix to generic names (e.g. Feature_1, Feature_2, etc.)
featureName.map <- data.frame("feature.name"=colnames(featureMatrix))
featureMatrix <- setNames(featureMatrix, paste("Feature_", seq(1:(ncol(featureMatrix))), sep=""))
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
  ppData[['knnImp']]<- preProcess(rm.highNA, method=c("knnImpute", "center", "scale")) %>% predict(rm.highNA)
  ppData[['rm.allNA']] <- preProcess(rm.allNA, method=c("center", "scale")) %>% predict(rm.allNA)
  
  return(ppData)
}

ppData <- preProcess.MLdata(featureMatrix)
ppData_2class <- lapply(ppData, function(x) {cbind(response.var=as.factor(proteins.dt_2Class$Status), x)})
ppData_multiClass <- lapply(ppData, function(x) {cbind(response.var=as.factor(str_replace_all(proteins.dt_multiClass$Status, " ", ".")), x)})
rm(ppData)

## Remove highly correlated variables
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

ppData_2class <- lapply(ppData_2class, rm.highCor, cutoff=0.98)
ppData_multiClass <- lapply(ppData_multiClass, rm.highCor, cutoff=0.98)


#####################
generateModels <- function(chooseData, xClass) {
  ### NO data splitting
  train.data  <- chooseData
  
  ### Specify type of resampling: repeated 10-fold CV
  fitControl <- trainControl(
    method="repeatedcv",
    number=10, # 10-fold
    repeats=5, # repeated 5 times
    index=createMultiFolds(train.data$response.var, k = 10, times = 5), # 10-fold cv repeated 10 times
    savePredictions=T,
    classProbs = T,
    if (xClass=="2Class") {
      summaryFunction=twoClassSummary
    } else if (xClass=="multiClass") {
      summaryFunction=multiClassSummary
    }
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
                           list(nb = caret::train(x=subset(train.data, select=-c(response.var)), y=train.data$response.var,
                                                      method = "nb",
                                                      trControl = fitControl,
                                                      tuneLength=15))),
                         error = function(cond){
                           message('naive bayes failed')
                           message(cond)
                           return(model_list)
                         })
  
  
  
  return(list("modelList"=model_list, "train.data"=train.data))
}
#####################

modelList.multiClass <- lapply(ppData_multiClass, generateModels, xClass="multiClass")
modelList.2Class <- lapply(ppData_2class, generateModels, xClass="2Class")
saveRDS(list("multiClass"=modelList.multiClass, "2Class"=modelList.2Class, "featureName.map"=featureName.map), "OutputmodelList_classification_descriptive.RDS")