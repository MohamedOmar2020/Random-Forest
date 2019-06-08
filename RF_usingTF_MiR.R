#################################################################################
### Mohamed Omar
### 18/4/2019
### GOAL: Creating a rondom forest classifier for bladder cancer progression (Non-muscle invasive << Muscle-invasive)
   ### Using: TF_MiR genes (Luigi)
#################################################################################

# Clean the work space
rm(list = ls())

## settng the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/RF_Classifier")

## Load necessary libraries
library(randomForest)
library(pROC)
library(ROCR)
library(caret)
library(limma)
library(genefilter)
## Load the data
load("./Objs/ProgressionDataGood.rda")

## Load TF_MiR genes
TF_MiR <- load("/Users/mohamedomar/Desktop/TF_MiR.rda")
keepGns <- intersect(as.vector(myTSPs), rownames(mixTrainMat))

### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat)[keepGns, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat)[keepGns, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))
################################################################################
################################################################################
################################################################################
### Creating the classifier using the training set



## Filter out any variables (genes) that are not expressed or do not have enough variance to be informative in classification. 
## We will first take the values and un-log2 them, then filter out any genes according to following criteria: (1) At least 20% of samples should have raw intensity greater than 100; (2) The coefficient of variation (sd/mean) is between 0.7 and 10.
#X <- usedTrainMat
#ffun <- filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))

#filt <- genefilter(2^X,ffun)
#filt_Data <- usedTrainMat[filt,]

## transpose the matrix
predictor_data <- t(usedTrainMat)
predictor_names <- c(as.vector(rownames(usedTrainMat))) #gene symbol
colnames(predictor_data) <- predictor_names

## Setting the variable we are trying to predict as our target variable. In this case, it is Progression status.
## train group here is just the column containing the phenptype of interest (Progression vs NoProgression) from the phenotype table

#usedTrainGroup <- ordered(usedTrainGroup, levels=c("NoProgression", "Progression"))
target <- usedTrainGroup

## Finally we run the RF algorithm. 
## NOTE: use an ODD number for ntree. This is because when the forest is used on test data, ties are broken randomly. Having an odd number of trees avoids this issue.
## Use down-sampling to attempt to compensate for unequal class-sizes (less progression than noProgression).
tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)

######
#Data <- cbind(predictor_data, usedTrainGroup)
#Data <- as.data.frame(Data)
#Data$usedTrainGroup <- as.factor(Data$usedTrainGroup)

## Using undersampling of the majority class to compensate for the unBalanced classes
#library(ROSE)
#set.seed(333)
#Over <- ovun.sample(usedTrainGroup ~., data = Data, method = "over", N=300)$data   ## AUC:68,65

#levels(Over$usedTrainGroup) <- c("NoProgression", "Progression")
#Over$usedTrainGroup <- ordered(Over$usedTrainGroup, levels=c("Progression", "NoProgression"))
#table(Over$usedTrainGroup)

# Tunning RF model (to find the best mtry)
set.seed(333)
tuneRF(x=predictor_data, y=target, plot = TRUE, improve = 0.1, ntreeTry = 1001)

set.seed(333)
rf_output <- randomForest(x =predictor_data, y=target,importance = TRUE, ntree = 1001, mtry =46 ,proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
print(rf_output)
plot(rf_output)
## RandomForest calculates an importance measures for each variable.
rf_importances <- importance(rf_output, scale=FALSE)

## Overview of the classifier's performance by generating a confusion table to allow calculation of sensitivity, specificity, accuracy, etc.
#confusion_train <- rf_output$confusion
#sensitivity_train <- (confusion_train[2,2]/(confusion_train[2,2]+confusion_train[1,2]))*100
#specificity_train <- (confusion_train[1,1]/(confusion_train[1,1]+confusion_train[2,1]))*100
#overall_error_train <- rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
#class1_error_train <- paste(rownames(confusion_train)[1]," error rate= ",confusion_train[1,3], sep="")
#class2_error_train <- paste(rownames(confusion_train)[2]," error rate= ",confusion_train[2,3], sep="")
#overall_accuracy_train <- 100-overall_error_train
################

## Predict in the training data
Prediction_Train <- predict(rf_output, predictor_data)
confusion_train <- confusionMatrix(Prediction_Train, usedTrainGroup, positive = "Progression")
confusion_train

## Create a representation of the top 50 variables categorized by importance.
pdf(file="varimp_TF_MiR.pdf")
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

## An MDS plot provides a sense of the separation of classes.
pdf(file="MDS_TF_MiR.pdf")
target_labels=as.vector(target)
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")
dev.off()



#train_pred <- predict(rf_output, newdata = Over, type = "vote")
#roc(Over$usedTrainGroup, train_pred[,2], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE)

## Create ROC curve and calculate AUC. We use the Progression/noProgression vote fractions as predictive variable. The ROC curve is generated by stepping through different thresholds for calling Progression vs noProgression.
predictions <- as.vector(rf_output$votes[,1])
pred <- prediction(predictions,target)
#First calculate the AUC value
perf_AUC_train <- performance(pred,"auc")
AUC_train <- perf_AUC_train@y.values[[1]]
#Then, plot the actual ROC curve
perf_ROC_train <- performance(pred,"tpr","fpr")
pdf(file="ROC_train_TF_MiR.pdf")
plot(perf_ROC_train, main="ROC plot train TF_MiR")
text(0.5,0.5,paste("AUC = ",format(AUC_train, digits=5, scientific=FALSE)))
dev.off()


################################################################################ 
################################################################################ 
## Testing the classifier using the testing set

## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(usedTestMat)
predictor_names2 <- c(as.vector(rownames(usedTestMat))) #gene symbol
colnames(predictor_data2) <- predictor_names2

#testGroup <- factor(do.call("c", groupProgression_test))
#levels(testGroup) <- c("Progression", "NoProgression")
#usedTestGroup <- ordered(usedTestGroup, levels=c("NoProgression", "Progression"))

## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(rf_output$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses <- predict(rf_output, predictor_data2, type="response")
RF_predictions_votes <- predict(rf_output, predictor_data2, type="vote")

## As before, the following lines will give an overview of the classifier's performance. This time instead of estimated performance from out-of-bag (OOB) testing we are measuring actual performance on an independent test set.
#confusion <- table(usedTestGroup, RF_predictions_responses)
#sensitivity_test <- (confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
#specificity_test <- (confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
#overall_error_test <- ((confusion[1,2]+confusion[2,1])/sum(confusion))*100
#overall_accuracy_test <- ((confusion[1,1]+confusion[2,2])/sum(confusion))*100
#class1_error_test <- confusion[1,2]/(confusion[1,1]+confusion[1,2])
#class2_error_test <- confusion[2,1]/(confusion[2,2]+confusion[2,1])
#####################

### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses, usedTestGroup, positive = "Progression")
confusion_test

## Create variables for the known target class and predicted class probabilities.
target <- usedTestGroup
progression_scores <- RF_predictions_votes[,"Progression"]

## Once again, we will create an ROC curve and calculate the area under it (AUC). Recall that we use the Progression/NoProgression vote fractions as predictive variable. The ROC curve is generated by stepping through different thresholds for calling Progression vs NoProgression.
pred <- prediction(progression_scores,target)
perf_AUC <- performance(pred,"auc")
AUC_test <- perf_AUC@y.values[[1]]
AUC_out <- paste("AUC=",AUC_test,sep="")

perf_ROC_test <- performance(pred,"tpr","fpr")
pdf(file="ROC_test_TF_MiR.pdf")
plot(perf_ROC_test, main="ROC plot TF_MiR")
text(0.5,0.5,paste("AUC = ",format(AUC_test, digits=5, scientific=FALSE)))
dev.off()

#########
## Save RF classifier
save(rf_output, file = "./Objs/RF_classifier_TF_MiR.rda")

## Save RF importances
save(rf_importances, file = "./Objs/RF_importances_TF_MiR.rda")

