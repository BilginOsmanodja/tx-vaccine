install.packages("aod")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("pastecs")
install.packages("psych")
install.packages("RColorBrewer")
install.packages("pROC")
install.packages("randomForest")
install.packages("glmnet")
install.packages("gbm")
install.packages("caret")
install.packages('lme4')
install.packages('class')
install.packages("mice")
install.packages("boot")
install.packages("tidyverse")
install.packages("janitor")
install.packages("plyr")
install.packages("mltools")
install.packages("data.table")

library(data.table)
library(mltools)
library(plyr)
library(janitor)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(aod)
library(ggplot2)
library(pastecs)
library(psych)
library(pROC)
library(randomForest)
library(glmnet)
library(gbm)
library(lme4)
library(caret)
library(class)
library(doParallel)
library(mice)
library(boot)

##### 1. Data Preparation Development Dataset
##### 2. Model Fitting on entire development dataset
##### 3. Internal Validation using Resampling
##### 4. Data Preparataion Validation Datasets
##### 5. External Validation point estimate and bootstrapping

##### 1. Data Preparation Development Dataset

# Read Development Dataset
set.seed(42)
data <- read.table(file = "E:\\Impfungen_modell\\Final\\Dataset_345_Vacc_Final.csv", sep=";", strip.white=TRUE, header=TRUE)

# Data preparation Development Dataset
data$PatientID <- as.numeric(data$PatientID)
data[data == "n.berechenbar"] <- 0
data$Albuminuria <- pmax(data$AlbuminKSU, data$AlbuminKUR)
data$CNI <- pmax(data$Tac, data$CyA)
data$DialysisYears[is.na(data$DialysisYears)] <- 0
data$DateOfVaccination <- as.Date(data$DateOfVaccination, format = "%d.%m.%Y")
data$BirthDate <- as.Date(data$BirthDate, format = "%d.%m.%Y")
data$DeathDate <- as.Date(data$DeathDate, format = "%d.%m.%Y")
data$LatestTxDate <- as.Date(data$LatestTxDate, format = "%d.%m.%Y")
data$TimeSinceLastVacc <- 0

# Compute TimeSinceLastVacc
akt_id = 65
for (i in 2:(nrow(data)-1)){
  if (data$PatientID[i+1] != akt_id){
    akt_id <- data$PatientID[i+1]
  }
  if ((data$PatientID[i]==data$PatientID[i-1]) & (data$VaccinationNo[i-1]==data$VaccinationNo[i]-1)) {
    data$TimeSinceLastVacc[i] <- difftime(data$DateOfVaccination[i], data$DateOfVaccination[i-1], units="days")
  }
}

# Select only 3rd or 4th Vaccination  with negative, low positive antibody level
data_model <- subset(data , (data$IgGpreDoV=="neg" | data$IgGpreDoV=="low pos") & (data$exclusion==0) & data$VaccinationNo>=3 & data$VaccinationNo<5)
data_model$VaccinationNo[data_model$VaccinationNo == 3] <- 1
data_model$VaccinationNo[data_model$VaccinationNo == 4] <- 0
names(data_model)[names(data_model) == "VaccinationNo"] <- "VaccinationNo_3"

# Define dataset without lymphocyte count (19 predictor variables) and with lymphocyte count (20 predictor variables)
data_model_2 <- data_model[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                             'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS', 
                             'GFR4vHPSCH', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
data_lympho <- data_model[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                             'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS', 
                             'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

# Data preparation
data_model_2$Sex[data_model_2$Sex == "m"] <- 0 # male = 0
data_model_2$Sex[data_model_2$Sex == "w"] <- 1 # female = 1
data_model_2$IgGpreDoV[data_model_2$IgGpreDoV == "neg"] <- 0
data_model_2$IgGpreDoV[data_model_2$IgGpreDoV == "low pos"] <- 1
data_model_2$IgGpostDoV[data_model_2$IgGpostDoV == "neg"] <- 0
data_model_2$IgGpostDoV[data_model_2$IgGpostDoV == "pos"] <- 1
data_model_2$IgGpostDoV[data_model_2$IgGpostDoV == "low pos"] <- 0
data_model_2$GFR4vHPSCH[data_model_2$GFR4vHPSCH == ">90"] <- 100 # eGFR > 90 -> 100
data_model_2$GFR4vHPSCH[data_model_2$GFR4vHPSCH == "<15"] <- 10 # eGFR < 15 -> 10
data_model_2[data_model_2 == "negativ"] <- 0
data_model_2$Sex <- as.numeric(data_model_2$Sex)
data_model_2$BMI <- as.numeric(data_model_2$BMI)
data_model_2$Albuminuria <- as.numeric(data_model_2$Albuminuria)
data_model_2$GFR4vHPSCH <- as.numeric(data_model_2$GFR4vHPSCH)
data_model_2$IgGpreDoV <- as.numeric(data_model_2$IgGpreDoV)
data_model_2$IgGpostDoV <- as.numeric(data_model_2$IgGpostDoV)
data_model_2$MPADose <- data_model_2$MPADose/1000 # MPA Dose in g
data_model_2$Albuminuria <- data_model_2$Albuminuria/1000 # Albumine-Creatinine Ratio in g/g

data_lympho$Sex[data_lympho$Sex == "m"] <- 0
data_lympho$Sex[data_lympho$Sex == "w"] <- 1
data_lympho$IgGpreDoV[data_lympho$IgGpreDoV == "neg"] <- 0
data_lympho$IgGpreDoV[data_lympho$IgGpreDoV == "low pos"] <- 1
data_lympho$IgGpostDoV[data_lympho$IgGpostDoV == "neg"] <- 0
data_lympho$IgGpostDoV[data_lympho$IgGpostDoV == "pos"] <- 1
data_lympho$IgGpostDoV[data_lympho$IgGpostDoV == "low pos"] <- 0
data_lympho$GFR4vHPSCH[data_lympho$GFR4vHPSCH == ">90"] <- 100
data_lympho$GFR4vHPSCH[data_lympho$GFR4vHPSCH == "<15"] <- 10
data_lympho[data_lympho == "negativ"] <- 0
data_lympho$Sex <- as.numeric(data_lympho$Sex)
data_lympho$BMI <- as.numeric(data_lympho$BMI)
data_lympho$LymphoEB <- as.numeric(data_lympho$LymphoEB)
data_lympho$Albuminuria <- as.numeric(data_lympho$Albuminuria)
data_lympho$GFR4vHPSCH <- as.numeric(data_lympho$GFR4vHPSCH)
data_lympho$IgGpreDoV <- as.numeric(data_lympho$IgGpreDoV)
data_lympho$MPADose <- data_lympho$MPADose/1000
data_lympho$Albuminuria <- data_lympho$Albuminuria/1000
data_lympho$PatientID <- as.character(data_lympho$PatientID) # to exclude from imputation
data_lympho$IgGpostDoV <- as.character(data_lympho$IgGpostDoV) # to exclude from imputation

# Single Imputation in Development Dataset for Testing Purposes -> Conclusion Complete Case Analysis
imp_dev <- mice(data_lympho, m=5, maxit = 50, method = 'pmm', seed = 42)
summary(imp_dev)
data_model_3_imp <- complete(imp_dev,1)
data_model_3_imp$IgGpostDoV <- as.integer(data_model_3_imp$IgGpostDoV)

data_model_3_imp <- data_model_3_imp[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                        'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS', 
                                        'GFR4vHPSCH', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

# Complete Case Analysis and Descriptive Stastistics
data_lympho$IgGpostDoV <- as.integer(data_lympho$IgGpostDoV) # 
data_model_3 <- na.omit(data_model_2) # Dataset without lymphocyte count as predictor variable
data_lympho_2 <- na.omit(data_lympho)
length(unique(data_lympho_2$PatientID)) # Show number of unique Patients
describe(data_lympho_2)
summary(data_lympho_2)

nrow(data_lympho_2[data_lympho_2$VaccinationNo_3==0 & data_lympho_2$IgGpostDoV==0, ])


# Remove PatientID as variable
data_model_3 <- data_model_3[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS', 
                                'GFR4vHPSCH', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
data_lympho_2 <- data_lympho_2[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                  'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS', 
                                  'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

# Check for normality -> Report Mean or Median respectively
#hist(data_lympho_2$Albuminuria) # Histogramm für jede Variable
#qqnorm(data_lympho_2$Albuminuria, pch = 1, frame = FALSE) # QQ Plot für jede Variable
#qqline(data_lympho_2$Albuminuria, col = "steelblue", lwd = 2) # QQ Plot - Linie

##### 2. Model Fitting on entire development dataset
#### 2.1 Logistic Regression with all predictor variables (only dataset with lymphocyte count)
lr_lympho <- glm(IgGpostDoV ~ ., data = data_lympho_2, family = "binomial")
summary(lr_lympho)
save(lr_lympho, file="E:\\Impfungen_modell\\Final\\LR_Revision.RData")
rm(lr_lympho)
load("E:\\Impfungen_modell\\Final\\LR_Revision.RData")
# Predict Scores to perform ROC Analysis and Calculate Thresholds

lr_score_lympho <- predict(lr_lympho, newdata = data_lympho_2)
lr_roc_lympho <- roc(data_lympho_2$IgGpostDoV ~ as.numeric(lr_score_lympho), plot = TRUE, print.auc = TRUE)
coords_lr <- coords(lr_roc_lympho, "best", best.method="closest.topleft", transpose = T)
cutoff_lr <- as.numeric(coords_lr[1])

#### 2.2 LASSO Regularized Logistic Regression (only dataset with lymphocyte count)
x_vars_lympho <- model.matrix(IgGpostDoV~. , data_lympho_2)[,-1]
y_var_lympho <- as.integer(data_lympho_2$IgGpostDoV) # since IgGpostDoV has been converted to character to exclude it from Imputation
lassologcv_lympho <- cv.glmnet(x_vars_lympho, y_var_lympho, family="binomial", alpha = 1,   nfolds = 5, type.measure = "auc", set.seed(42))
plot(lassologcv_lympho)

### 2.2.1 LASSO Min Model
lasso_min_lympho <- glmnet(x_vars_lympho, y_var_lympho, alpha = 1,  family="binomial", lambda = lassologcv_lympho$lambda.min)
coef(lasso_min_lympho)
# Save LASSO-Min Model
save(lasso_min_lympho, file="E:\\Impfungen_modell\\Final\\LASSO_Min_Revision.RData")
rm(lasso_min_lympho)
load("E:\\Impfungen_modell\\Final\\LASSO_Min_Revision.RData")
# Predict Scores to perform ROC Analysis and Calculate Thresholds
lasso_min_score_lympho <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = x_vars_lympho)
lasso_min_roc_lympho <- roc(y_var_lympho ~ as.numeric(lasso_min_score_lympho), plot = TRUE, print.auc = TRUE)
coords_min <- coords(lasso_min_roc_lympho, "best", best.method="closest.topleft", transpose = T)
cutoff_min <- as.numeric(coords_min[1])

### 2.2.2 LASSO 1SE Model
lasso_1se_lympho <- glmnet(x_vars_lympho, y_var_lympho, alpha = 1,  family="binomial", lambda = lassologcv_lympho$lambda.1se)
coef(lasso_1se_lympho)
# Save LASSO-1SE Model
save(lasso_1se_lympho, file="E:\\Impfungen_modell\\Final\\LASSO_1SE_Revision.RData")
rm(lasso_1se_lympho)
load("E:\\Impfungen_modell\\Final\\LASSO_1SE_Revision.RData")
# Predict Scores to perform ROC Analysis and Calculate Thresholds
lasso_1se_score_lympho <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = x_vars_lympho)
lasso_1se_roc_lympho <- roc(y_var_lympho ~ as.numeric(lasso_1se_score_lympho), plot = TRUE, print.auc = TRUE)
coords_1se <- coords(lasso_1se_roc_lympho, "best", best.method="closest.topleft", transpose = T)
cutoff_1se <- as.numeric(coords_1se[1])

# Risk Probability Cutoff for Online Calculator
prob_cutoff_min <- (1/(1+exp(-cutoff_min)))
prob_cutoff_1se <- (1/(1+exp(-cutoff_1se)))

#### 2.3 RF fitted on the entire development dataset
data_rf_lympho <- data_lympho_2
data_rf_lympho$IgGpostDoV <- as.factor(data_rf_lympho$IgGpostDoV)
control <- trainControl(method="repeatedcv", number=5, repeats=2, search="random")
rf_random_lympho <- train(IgGpostDoV~., data=data_rf_lympho, method="rf", metric="Accuracy", tuneLength=15, ntree=1000, trControl=control, importance=FALSE, proximity=FALSE)
rf_lympho <- randomForest(IgGpostDoV ~., data = data_rf_lympho, ntree=1000, mtry=rf_random_lympho$bestTune$mtry, importance=TRUE, proximity=FALSE)
# Cutoff for RF = 0.5
save(rf_lympho, file="E:\\Impfungen_modell\\Final\\RF_Lympho_Revision.RData")
rm(rf_lympho)
load("E:\\Impfungen_modell\\Final\\RF_Lympho_Revision.RData")


#### 2.4 GBRT fitted on the entire development dataset
grid <- expand.grid(.n.trees=seq(300,900,by=200),.interaction.depth=seq(2,16,by=2),.shrinkage=c(0.001, 0.01, 0.1), .n.minobsinnode=10)
control_gbm <- trainControl(method="repeatedcv", number=5, repeats=2)
gbm_train_lympho<-train(IgGpostDoV~.,data=data_rf_lympho, method='gbm', tuneGrid = grid)
gbm_lympho<-gbm(IgGpostDoV~.,data=data_rf_lympho,n.trees = gbm_train_lympho$bestTune$n.trees,interaction.depth = gbm_train_lympho$bestTune$interaction.depth,
                        shrinkage = gbm_train_lympho$bestTune$shrinkage, n.minobsinnode=gbm_train_lympho$bestTune$n.minobsinnode, distribution = 'gaussian')
gbm_score_lympho <- predict(gbm_lympho,newdata = data_rf_lympho)
gbm_roc_lympho <- roc(data_rf_lympho$IgGpostDoV ~ gbm_score_lympho, plot = TRUE, print.auc = TRUE)
coords_gbm <- coords(gbm_roc_lympho, "best", best.method="closest.topleft", transpose = T)
cutoff_gbm <- as.numeric(coords_gbm[1])
save(gbm_lympho, file="E:\\Impfungen_modell\\Final\\GBM_Lympho_Revision.RData")
rm(gbm_lympho)
load("E:\\Impfungen_modell\\Final\\GBM_Lympho_Revision.RData")

##### 3. Internal Validation

#  Define empty dataframes for Performance Metrics during Resampling
df_lr <- data.frame()
df_lr_imp <- data.frame()
df_lr_lympho <- data.frame() # Logistische Regression mit Lymphos
df_lr_lasso_min_lympho <- data.frame() # Logistische Regression mit Lymphos mit Lasso-Regularization mit minimum lambda in der CV
df_lr_lasso_1se_lympho <- data.frame() # Logistische Regression mit Lymphos mit Lasso-Regularization mit 1 SE lambda in der CV
df_rf_lympho <- data.frame() # Random Forest Model mit Lympho
df_gbm_lympho <- data.frame() # Gradient Boosted Model mit Lympho

# Define variables and dataframes to store coefficients during resampling to assess model stability.
df_coef_lasso_min <- data.frame(matrix(NA, nrow = 0, ncol = 21))
colnames(df_coef_lasso_min)<- c('X.Intercept.','IgGpreDoV', 'Sex','Age','BMI','mRNA','Retransplantation','TransplantAge', 'DialysisYears', 'Diabetes',
                                'Steroid','Belatacept','MPADose','CNI','MoreThan2IS','GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','VaccinationNo_3')

df_coef_lasso_1se <- data.frame(matrix(NA, nrow = 0, ncol = 21))
colnames(df_coef_lasso_1se)<- c('X.Intercept.','IgGpreDoV', 'Sex','Age','BMI','mRNA','Retransplantation','TransplantAge', 'DialysisYears', 'Diabetes',
                                'Steroid','Belatacept', 'MPADose','CNI','MoreThan2IS','GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','VaccinationNo_3')


# Loop with 100 resampling steps
# Each loop produces train/test splits of 70:30 for final dataset, single imputed dataset and dataset without lymphocyte count
# LR is fitted on all three datasets to test, if imputation or using no leukocyte count could be of benefit
# LASSO LR, RF, and GBRT are fitted on the data_lympho_2 dataset only, since this was chosen as the final development dataset 
for (i in 1:100){
  set.seed(i) #Random Seed changes with every loop, so 100 different train/test splits are achieved
  print(i)
  
  # Define indices for Train/Test in imputed Development dataset (70:30 Split)
  train_ind_imp <- sample(seq_len(nrow(data_model_3_imp)), size = nrow(data_model_3_imp)*0.7)
  train_imp <- data_model_3_imp[train_ind_imp, ]
  test_imp <- data_model_3_imp[-train_ind_imp, ]

  # Define indices for Train/Test in development dataset w/o lymphocyte count (70:30 Split)
  train_ind <- sample(seq_len(nrow(data_model_3)), size = nrow(data_model_3)*0.7)
  train <- data_model_3[train_ind, ]
  test <- data_model_3[-train_ind, ]
  
  # Define indices for Train/Test in final development dataset with lymphocyte count (70:30 Split)
  train_ind_lympho <- sample(seq_len(nrow(data_lympho_2)), size = nrow(data_lympho_2)*0.7)
  train_lympho <- data_lympho_2[train_ind_lympho,]
  test_lympho <- data_lympho_2[-train_ind_lympho,]
  
  ### Logistic Regression in imputed Development dataset
  lr_imp <- glm(IgGpostDoV ~ ., data = train_imp, family = "binomial")
  coef(lr_imp)
  lr_score_imp <- predict(lr_imp, newdata = test_imp, type = "response") # Predict on Testdata
  lr_roc_imp <- roc(test_imp$IgGpostDoV ~ lr_score_imp, plot = TRUE, print.auc = TRUE) # ROC-Analysis 
  lr_metric_imp <- coords(lr_roc_imp, "best", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), best.method="closest.topleft", transpose = T)
  lr_auc_imp <- as.numeric(lr_roc_imp$auc) # ROC-AUC as numeric
  lr_metric_imp[6]<-lr_auc_imp
  df_lr_imp <- rbind(df_lr_imp, lr_metric_imp)
  
  ### Logistic Regression in development dataset w/o lymphocyte count
  lr <- glm(IgGpostDoV ~ ., data = train, family = "binomial")
  lr_score <- predict(lr, newdata = test, type = "response") # Predict on Testdata
  lr_roc <- roc(test$IgGpostDoV ~ lr_score, plot = TRUE, print.auc = TRUE) # ROC-Analysis 
  lr_metric<- coords(lr_roc, "best", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), best.method="closest.topleft", transpose = T)
  lr_auc <- as.numeric(lr_roc$auc) # ROC-AUC as numeric
  lr_metric[6]<-lr_auc
  df_lr <- rbind(df_lr, lr_metric)
  
  ### Logistic Regression in final development dataset with lymphocyte count
  lr_lympho <- glm(IgGpostDoV ~., data = train_lympho, family = "binomial")
  lr_score_lympho <- predict(lr_lympho, newdata = test_lympho, type = "response")
  lr_roc_lympho <- roc(test_lympho$IgGpostDoV ~ lr_score_lympho, plot = TRUE, print.auc = TRUE)
  lr_metric_lympho <- coords(lr_roc_lympho, "best", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), best.method="closest.topleft", transpose = T)
  lr_auc_lympho <- as.numeric(lr_roc_lympho$auc)
  lr_metric_lympho[6]<-lr_auc_lympho
  df_lr_lympho <- rbind(df_lr_lympho, lr_metric_lympho)

  ### Lasso Logistic Regression in final development dataset with lymphocyte count
  x_vars_lympho <- model.matrix(IgGpostDoV~. , data_lympho_2)[,-1]
  y_var_lympho <- as.integer(data_lympho_2$IgGpostDoV)
  # Defining Hyperparameter Lambda
  cv_output_lympho <- cv.glmnet(x_vars_lympho[train_ind_lympho, ], y_var_lympho[train_ind_lympho], family="binomial", alpha = 1, nfolds = 5, type.measure = "auc", set.seed(i))
  
  # LASSO-Min Model
  lambda_min_lympho <- cv_output_lympho$lambda.min
  lasso_min_lympho <- glmnet(x_vars_lympho[train_ind_lympho,], y_var_lympho[train_ind_lympho], alpha = 1, lambda = lambda_min_lympho, family="binomial")
  lr_lasso_min_score_lympho <- predict(lasso_min_lympho, s = lambda_min_lympho, newx = x_vars_lympho[-train_ind_lympho,])
  lr_lasso_min_roc_lympho <- roc(y_var_lympho[-train_ind_lympho] ~ as.numeric(lr_lasso_min_score_lympho), plot = TRUE, print.auc = TRUE)
  lr_metric_lasso_min_lympho <- coords(lr_lasso_min_roc_lympho, "best", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), best.method="closest.topleft", transpose = T)
  lr_lasso_min_auc_lympho <- as.numeric(lr_lasso_min_roc_lympho$auc)
  lr_metric_lasso_min_lympho[6]<-lr_lasso_min_auc_lympho
  df_lr_lasso_min_lympho <- rbind(df_lr_lasso_min_lympho, lr_metric_lasso_min_lympho)
  # Save Coefficients from LASSO-Min Model
  coefs_lasso_min <- coef(lasso_min_lympho, s = "lambda.min")
  coefs_lasso_min_2<- t(data.frame(name = coefs_lasso_min@Dimnames[[1]][coefs_lasso_min@i + 1], coefficient = coefs_lasso_min@x))
  coefs_lasso_min_2 <- data.frame(janitor::row_to_names(coefs_lasso_min_2, row_number = 1))
  df_coef_lasso_min <- rbind.fill(df_coef_lasso_min, coefs_lasso_min_2)
  
  # LASSO-1SE Model
  lambda_1se_lympho <- cv_output_lympho$lambda.1se
  lasso_1se_lympho <- glmnet(x_vars_lympho[train_ind_lympho,], y_var_lympho[train_ind_lympho], alpha = 1, lambda = lambda_1se_lympho, family="binomial")
  lr_lasso_1se_score_lympho <- predict(lasso_1se_lympho, s = lambda_1se_lympho, newx = x_vars_lympho[-train_ind_lympho,])
  lr_lasso_1se_roc_lympho <- roc(y_var_lympho[-train_ind_lympho] ~ as.numeric(lr_lasso_1se_score_lympho), plot = TRUE, print.auc = TRUE)
  lr_metric_lasso_1se_lympho <- coords(lr_lasso_1se_roc_lympho, "best", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), best.method="closest.topleft", transpose = T)
  lr_lasso_1se_auc_lympho <- as.numeric(lr_lasso_1se_roc_lympho$auc)
  lr_metric_lasso_1se_lympho[6]<-lr_lasso_1se_auc_lympho
  df_lr_lasso_1se_lympho <- rbind(df_lr_lasso_1se_lympho, lr_metric_lasso_1se_lympho)
  # Save Coefficients from LASSO-1SE Model
  coefs_lasso_1se <- coef(lasso_1se_lympho, s = "lambda.1se")
  coefs_lasso_1se_2<- t(data.frame(name = coefs_lasso_1se@Dimnames[[1]][coefs_lasso_1se@i + 1], coefficient = coefs_lasso_1se@x))
  coefs_lasso_1se_2 <- data.frame(janitor::row_to_names(coefs_lasso_1se_2, row_number = 1))
  df_coef_lasso_1se <- rbind.fill(df_coef_lasso_1se, coefs_lasso_1se_2)
  
  ### RandomForest in final development dataset with lymphocyte count
  train_rf_lympho <- train_lympho
  train_rf_lympho$IgGpostDoV <- as.factor(train_rf_lympho$IgGpostDoV)
  control <- trainControl(method="repeatedcv", number=5, repeats=2, search="random")
  mtry <- sqrt(ncol(train_rf_lympho)-1)
  rf_random_lympho <- train(IgGpostDoV~., data=train_rf_lympho, method="rf", metric="Accuracy", tuneLength=15, ntree=1000, trControl=control, importance=FALSE, proximity=FALSE)
  rf_lympho <- randomForest(IgGpostDoV ~., data = train_lympho, ntree=1000, mtry=rf_random_lympho$bestTune$mtry, importance=FALSE, proximity=FALSE, type="class")
  rf_score_lympho <- predict(rf_lympho, newdata=test_lympho, type="class")
  rf_roc_lympho <- roc(test_lympho$IgGpostDoV ~ rf_score_lympho, plot = TRUE, print.auc = TRUE)
  rf_metric_lympho <- coords(rf_roc_lympho, 0.5, ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), input="threshold", transpose = T)
  rf_auc_lympho <- as.numeric(rf_roc_lympho$auc)
  rf_metric_lympho[6] <- rf_auc_lympho
  df_rf_lympho <- rbind(df_rf_lympho, rf_metric_lympho)
  
  ### GBM in final development dataset with lymphocyte count
  grid <- expand.grid(.n.trees=seq(300,900,by=200),.interaction.depth=seq(2,16,by=2),.shrinkage=c(0.001, 0.01, 0.1), .n.minobsinnode=10)
  control <- trainControl(method="repeatedcv", number=5, repeats=2)
  gbm.train_lympho<-train(IgGpostDoV~.,data=train_rf_lympho, method='gbm', tuneGrid = grid)
  gbm.postDoV_lympho<-gbm(IgGpostDoV~.,data=train_lympho,n.trees = gbm.train_lympho$bestTune$n.trees,interaction.depth = gbm.train_lympho$bestTune$interaction.depth,
                          shrinkage = gbm.train_lympho$bestTune$shrinkage, n.minobsinnode=gbm.train_lympho$bestTune$n.minobsinnode, distribution = 'gaussian')
  gbm.test_lympho<-predict(gbm.postDoV_lympho,newdata = test_lympho)
  gbm_roc_lympho <- roc(test_lympho$IgGpostDoV ~ gbm.test_lympho, plot = TRUE, print.auc = TRUE)
  gbm_metric_lympho <- coords(gbm_roc_lympho, "best", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), best.method="closest.topleft", transpose = T)
  gbm_auc_lympho <- as.numeric(gbm_roc_lympho$auc)
  gbm_metric_lympho[6]=gbm_auc_lympho
  df_gbm_lympho <- rbind(df_gbm_lympho, gbm_metric_lympho)
}

# Rename Columns
df_lr_imp <- setNames(df_lr_imp, c("sens","spec","acc","ppv","npv","auc"))
df_lr <- setNames(df_lr, c("sens","spec","acc","ppv","npv","auc"))
df_lr_lympho <- setNames(df_lr_lympho, c("sens","spec","acc","ppv","npv","auc"))
df_lr_lasso_min_lympho <- setNames(df_lr_lasso_min_lympho, c("sens","spec","acc","ppv","npv","auc"))
df_lr_lasso_1se_lympho <- setNames(df_lr_lasso_1se_lympho, c("sens","spec","acc","ppv","npv","auc"))
df_rf_lympho <- setNames(df_rf_lympho, c("sens","spec","acc","ppv","npv","auc"))
df_gbm_lympho <- setNames(df_gbm_lympho, c("sens","spec","acc","ppv","npv","auc"))

# Add Column "class" with model name
df_lr_imp$class <- "LR (imputed)"
df_lr$class <- "LR (no lympho)"
df_lr_lympho$class <- "LR (lympho)"
df_lr_lasso_min_lympho$class <- "LASSO-Min LR (lympho)"
df_lr_lasso_1se_lympho$class <- "LASSO-1SE LR (lympho)"
df_rf_lympho$class <- "RF (lympho)"
df_gbm_lympho$class <- "GBRT (lympho)"

#Save metrics
write.csv(df_lr_imp,"E:\\Impfungen_modell\\Final\\df_lr_imp_rev.csv", row.names = FALSE)
write.csv(df_lr,"E:\\Impfungen_modell\\Final\\df_lr_rev.csv", row.names = FALSE)
write.csv(df_lr_lympho,"E:\\Impfungen_modell\\Final\\df_lr_lympho_rev.csv", row.names = FALSE)
write.csv(df_lr_lasso_min_lympho,"E:\\Impfungen_modell\\Final\\df_lr_lasso_min_lympho_rev.csv", row.names = FALSE)
write.csv(df_lr_lasso_1se_lympho,"E:\\Impfungen_modell\\Final\\df_lr_lasso_1se_lympho_rev.csv", row.names = FALSE)
write.csv(df_rf_lympho,"E:\\Impfungen_modell\\Final\\df_rf_lympho_rev.csv", row.names = FALSE)
write.csv(df_gbm_lympho,"E:\\Impfungen_modell\\Final\\df_gbm_lympho_rev.csv", row.names = FALSE)
write.csv(df_coef_lasso_1se,"E:\\Impfungen_modell\\Final\\df_coef_lasso_1se_rev.csv", row.names = FALSE)
write.csv(df_coef_lasso_min,"E:\\Impfungen_modell\\Final\\df_coef_lasso_min_rev.csv", row.names = FALSE)


# Merge dataframes for internal validation Table/Figure
df_intern <- rbind(df_lr_lympho, df_lr_lasso_min_lympho)
df_intern <- rbind(df_intern, df_lr_lasso_1se_lympho)
df_intern <- rbind(df_intern, df_rf_lympho)
df_intern <- rbind(df_intern, df_gbm_lympho)

# Merge dataframes for supplement- comparison imputation vs. no lymphocyte count vs. final development dataset
df_supp <- rbind(df_lr_imp, df_lr)
df_supp <- rbind(df_supp, df_lr_lympho)


# Extract Mean, Median und 95%CI of Performance during Internal Validation for Tables
#LASSO-1SE
quantile(df_lr_lasso_1se_lympho$auc, prob=c(.025,.5,.975))
mean(df_lr_lasso_1se_lympho$auc)
quantile(df_lr_lasso_1se_lympho$sens, prob=c(.025,.5,.975))
mean(df_lr_lasso_1se_lympho$sens)
quantile(df_lr_lasso_1se_lympho$spec, prob=c(.025,.5,.975))
mean(df_lr_lasso_1se_lympho$spec)
quantile(df_lr_lasso_1se_lympho$acc, prob=c(.025,.5,.975))
mean(df_lr_lasso_1se_lympho$acc)
quantile(df_lr_lasso_1se_lympho$ppv, prob=c(.025,.5,.975))
mean(df_lr_lasso_1se_lympho$ppv)
quantile(df_lr_lasso_1se_lympho$npv, prob=c(.025,.5,.975))
mean(df_lr_lasso_1se_lympho$npv)

#LASSO-Min
quantile(df_lr_lasso_min_lympho$auc, prob=c(.025,.5,.975))
mean(df_lr_lasso_min_lympho$auc)
quantile(df_lr_lasso_min_lympho$sens, prob=c(.025,.5,.975))
mean(df_lr_lasso_min_lympho$sens)
quantile(df_lr_lasso_min_lympho$spec, prob=c(.025,.5,.975))
mean(df_lr_lasso_min_lympho$spec)
quantile(df_lr_lasso_min_lympho$acc, prob=c(.025,.5,.975))
mean(df_lr_lasso_min_lympho$acc)
quantile(df_lr_lasso_min_lympho$ppv, prob=c(.025,.5,.975))
mean(df_lr_lasso_min_lympho$ppv)
quantile(df_lr_lasso_min_lympho$npv, prob=c(.025,.5,.975))
mean(df_lr_lasso_min_lympho$npv)

#LR
quantile(df_lr_lympho$auc, prob=c(.025,.5,.975))
mean(df_lr_lympho$auc)
quantile(df_lr_lympho$sens, prob=c(.025,.5,.975))
mean(df_lr_lympho$sens)
quantile(df_lr_lympho$spec, prob=c(.025,.5,.975))
mean(df_lr_lympho$spec)
quantile(df_lr_lympho$acc, prob=c(.025,.5,.975))
mean(df_lr_lympho$acc)
quantile(df_lr_lympho$ppv, prob=c(.025,.5,.975))
mean(df_lr_lympho$ppv)
quantile(df_lr_lympho$npv, prob=c(.025,.5,.975))
mean(df_lr_lympho$npv)

#LR without lympho
quantile(df_lr$auc, prob=c(.025,.5,.975))
mean(df_lr$auc)
quantile(df_lr$sens, prob=c(.025,.5,.975))
mean(df_lr$sens)
quantile(df_lr$spec, prob=c(.025,.5,.975))
mean(df_lr$spec)
quantile(df_lr$acc, prob=c(.025,.5,.975))
mean(df_lr$acc)
quantile(df_lr$ppv, prob=c(.025,.5,.975))
mean(df_lr$ppv)
quantile(df_lr$npv, prob=c(.025,.5,.975))
mean(df_lr$npv)

#LR imputed
quantile(df_lr_imp$auc, prob=c(.025,.5,.975))
mean(df_lr_imp$auc)
quantile(df_lr_imp$sens, prob=c(.025,.5,.975))
mean(df_lr_imp$sens)
quantile(df_lr_imp$spec, prob=c(.025,.5,.975))
mean(df_lr_imp$spec)
quantile(df_lr_imp$acc, prob=c(.025,.5,.975))
mean(df_lr_imp$acc)
quantile(df_lr_imp$ppv, prob=c(.025,.5,.975))
mean(df_lr_imp$ppv)
quantile(df_lr_imp$npv, prob=c(.025,.5,.975))
mean(df_lr_imp$npv)


#RF
quantile(df_rf_lympho$auc, prob=c(.025,.5,.975))
mean(df_rf_lympho$auc)
quantile(df_rf_lympho$sens, prob=c(.025,.5,.975))
mean(df_rf_lympho$sens)
quantile(df_rf_lympho$spec, prob=c(.025,.5,.975))
mean(df_rf_lympho$spec)
quantile(df_rf_lympho$acc, prob=c(.025,.5,.975))
mean(df_rf_lympho$acc)
quantile(df_rf_lympho$ppv, prob=c(.025,.5,.975))
mean(df_rf_lympho$ppv)
quantile(df_rf_lympho$npv, prob=c(.025,.5,.975))
mean(df_rf_lympho$npv)

#GBRT
quantile(df_gbm_lympho$auc, prob=c(.025,.5,.975))
mean(df_gbm_lympho$auc)
quantile(df_gbm_lympho$sens, prob=c(.025,.5,.975))
mean(df_gbm_lympho$sens)
quantile(df_gbm_lympho$spec, prob=c(.025,.5,.975))
mean(df_gbm_lympho$spec)
quantile(df_gbm_lympho$acc, prob=c(.025,.5,.975))
mean(df_gbm_lympho$acc)
quantile(df_gbm_lympho$ppv, prob=c(.025,.5,.975))
mean(df_gbm_lympho$ppv)
quantile(df_gbm_lympho$npv, prob=c(.025,.5,.975))
mean(df_gbm_lympho$npv)

#Save data from internal validation
write.csv(df_intern,"E:\\Impfungen_modell\\Final\\Metrics_Internal_Val_Revision.csv", row.names = FALSE)

#Order for Plotting of internal validation
plot_order_int <- c("LR (lympho)", "LASSO-Min LR (lympho)", "LASSO-1SE LR (lympho)", "GBRT (lympho)", "RF (lympho)")
df_intern$class <- factor(df_intern$class, levels = plot_order_int)
#Boxplots of AUC
ggplot(data = df_intern, mapping = aes(x=class, y=auc)) + 
  geom_point(aes(),alpha=0.8) +
  geom_boxplot(fill="grey",color="black",alpha=0.2) + 
  theme_minimal()

# Order of Plotting for supplementary data
plot_order_supp <- c("LR (lympho)", "LR (imputed)", "LR (no lympho)")
df_supp$class <- factor(df_supp$class, levels = plot_order_supp)
# Boxplots of AUC
ggplot(data = df_supp, mapping = aes(x=class, y=auc)) + 
  geom_point(aes(),alpha=0.8) +
  geom_boxplot(fill="grey",color="black",alpha=0.2) + 
  theme_minimal()


### Analyze coefficients for LASSO LR

df_coef_lasso_1se = as.data.frame(sapply(df_coef_lasso_1se, as.numeric))
df_coef_lasso_min = as.data.frame(sapply(df_coef_lasso_min, as.numeric))

# Show how often how many variables were selected -> report most frequent variable count
hist(21-rowSums(is.na(df_coef_lasso_min))) # Result: 20
hist(21-rowSums(is.na(df_coef_lasso_1se))) # Result: 10

#Prepare Coefficient Data for Plotting LASSO-Min
lasso_min_intercept <- data.frame(df_coef_lasso_min[,"X.Intercept."])
lasso_min_intercept$Variable <- "Intercept"
colnames(lasso_min_intercept) <- c("Coefficient", "Variable")
lasso_min_IgGpreDoV <- data.frame(df_coef_lasso_min[,"IgGpreDoV"])
lasso_min_IgGpreDoV$Variable <- "IgGpreDoV"
colnames(lasso_min_IgGpreDoV) <- c("Coefficient", "Variable")
lasso_min_Sex <- data.frame(df_coef_lasso_min[,"Sex"])
lasso_min_Sex$Variable <- "Sex"
colnames(lasso_min_Sex) <- c("Coefficient", "Variable")
lasso_min_Age <- data.frame(df_coef_lasso_min[,"Age"])
lasso_min_Age$Variable <- "Age"
colnames(lasso_min_Age) <- c("Coefficient", "Variable")
lasso_min_BMI <- data.frame(df_coef_lasso_min[,"BMI"])
lasso_min_BMI$Variable <- "BMI"
colnames(lasso_min_BMI) <- c("Coefficient", "Variable")
lasso_min_mRNA <- data.frame(df_coef_lasso_min[,"mRNA"])
lasso_min_mRNA$Variable <- "mRNA"
colnames(lasso_min_mRNA) <- c("Coefficient", "Variable")
lasso_min_Retransplantation <- data.frame(df_coef_lasso_min[,"Retransplantation"])
lasso_min_Retransplantation$Variable <- "Retransplantation"
colnames(lasso_min_Retransplantation) <- c("Coefficient", "Variable")
lasso_min_TransplantAge <- data.frame(df_coef_lasso_min[,"TransplantAge"])
lasso_min_TransplantAge$Variable <- "TransplantAge"
colnames(lasso_min_TransplantAge) <- c("Coefficient", "Variable")
lasso_min_DialysisYears <- data.frame(df_coef_lasso_min[,"DialysisYears"])
lasso_min_DialysisYears$Variable <- "DialysisYears"
colnames(lasso_min_DialysisYears) <- c("Coefficient", "Variable")
lasso_min_Diabetes <- data.frame(df_coef_lasso_min[,"Diabetes"])
lasso_min_Diabetes$Variable <- "Diabetes"
colnames(lasso_min_Diabetes) <- c("Coefficient", "Variable")
lasso_min_Steroid <- data.frame(df_coef_lasso_min[,"Steroid"])
lasso_min_Steroid$Variable <- "Steroid"
colnames(lasso_min_Steroid) <- c("Coefficient", "Variable")
lasso_min_Belatacept <- data.frame(df_coef_lasso_min[,"Belatacept"])
lasso_min_Belatacept$Variable <- "Belatacept"
colnames(lasso_min_Belatacept) <- c("Coefficient", "Variable")
lasso_min_MPADose <- data.frame(df_coef_lasso_min[,"MPADose"])
lasso_min_MPADose$Variable <- "MPADose"
colnames(lasso_min_MPADose) <- c("Coefficient", "Variable")
lasso_min_CNI <- data.frame(df_coef_lasso_min[,"CNI"])
lasso_min_CNI$Variable <- "CNI"
colnames(lasso_min_CNI) <- c("Coefficient", "Variable")
lasso_min_MoreThan2IS <- data.frame(df_coef_lasso_min[,"MoreThan2IS"])
lasso_min_MoreThan2IS$Variable <- "MoreThan2IS"
colnames(lasso_min_MoreThan2IS) <- c("Coefficient", "Variable")
lasso_min_GFR4vHPSCH <- data.frame(df_coef_lasso_min[,"GFR4vHPSCH"])
lasso_min_GFR4vHPSCH$Variable <- "GFR4vHPSCH"
colnames(lasso_min_GFR4vHPSCH) <- c("Coefficient", "Variable")
lasso_min_LymphoEB <- data.frame(df_coef_lasso_min[,"LymphoEB"])
lasso_min_LymphoEB$Variable <- "LymphoEB"
colnames(lasso_min_LymphoEB) <- c("Coefficient", "Variable")
lasso_min_HbEB <- data.frame(df_coef_lasso_min[,"HbEB"])
lasso_min_HbEB$Variable <- "HbEB"
colnames(lasso_min_HbEB) <- c("Coefficient", "Variable")
lasso_min_Albuminuria <- data.frame(df_coef_lasso_min[,"Albuminuria"])
lasso_min_Albuminuria$Variable <- "Albuminuria"
colnames(lasso_min_Albuminuria) <- c("Coefficient", "Variable")
lasso_min_TimeSinceLastVacc <- data.frame(df_coef_lasso_min[,"TimeSinceLastVacc"])
lasso_min_TimeSinceLastVacc$Variable <- "TimeSinceLastVacc"
colnames(lasso_min_TimeSinceLastVacc) <- c("Coefficient", "Variable")
lasso_min_VaccinationNo_3 <- data.frame(df_coef_lasso_min[,"VaccinationNo_3"])
lasso_min_VaccinationNo_3$Variable <- "VaccinationNo_3"
colnames(lasso_min_VaccinationNo_3) <- c("Coefficient", "Variable")
#Combine Coefficient Data
df_coef_lasso_min_ggplot <- rbind(lasso_min_IgGpreDoV, lasso_min_Sex)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_Age)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_BMI)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_mRNA)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_Retransplantation)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_TransplantAge)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_DialysisYears)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_Diabetes)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_Steroid)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_Belatacept)
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_MPADose )
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_CNI )
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_MoreThan2IS )
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_GFR4vHPSCH )
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_LymphoEB )
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_HbEB )
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_Albuminuria )
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_TimeSinceLastVacc )
df_coef_lasso_min_ggplot <- rbind(df_coef_lasso_min_ggplot, lasso_min_VaccinationNo_3 )
summary(df_coef_lasso_min)

#Save Coefficient Data
write.csv(df_coef_lasso_min_ggplot,"S:\\C13\\MNP\\Digitale-Nephrologie\\vALID\\Impfungen\\Final\\Coefficients_LASSO_Min_Revision.csv", row.names = FALSE)

#Plot Coefficients
plot_order_coef_min <- c("IgGpreDoV", "MPADose", "GFR4vHPSCH", "VaccinationNo_3", "Belatacept", "LymphoEB",  "TransplantAge", "DialysisYears", "BMI", "HbEB", "Retransplantation", "MoreThan2IS", "Sex", "Albuminuria", "mRNA", "Age", "TimeSinceLastVacc", "Diabetes",  "Steroid", "CNI")
df_coef_lasso_min_ggplot$Variable <- factor(df_coef_lasso_min_ggplot$Variable, levels = plot_order_coef_min)
ggplot(data = df_coef_lasso_min_ggplot, mapping = aes(x=Variable, y=Coefficient)) + 
  geom_point(aes(),alpha=0.8) +
  geom_boxplot(fill="grey",color="black",alpha=0.2) + 
  theme_minimal() +
  coord_flip() +
  scale_x_discrete(limits = rev(levels(df_coef_lasso_min_ggplot$Variable)))

#Prepare Coefficient Data for Plotting LASSO-1SE
lasso_1se_intercept <- data.frame(df_coef_lasso_1se[,"X.Intercept."])
lasso_1se_intercept$Variable <- "Intercept"
colnames(lasso_1se_intercept) <- c("Coefficient", "Variable")
lasso_1se_IgGpreDoV <- data.frame(df_coef_lasso_1se[,"IgGpreDoV"])
lasso_1se_IgGpreDoV$Variable <- "IgGpreDoV"
colnames(lasso_1se_IgGpreDoV) <- c("Coefficient", "Variable")
lasso_1se_Sex <- data.frame(df_coef_lasso_1se[,"Sex"])
lasso_1se_Sex$Variable <- "Sex"
colnames(lasso_1se_Sex) <- c("Coefficient", "Variable")
lasso_1se_Age <- data.frame(df_coef_lasso_1se[,"Age"])
lasso_1se_Age$Variable <- "Age"
colnames(lasso_1se_Age) <- c("Coefficient", "Variable")
lasso_1se_BMI <- data.frame(df_coef_lasso_1se[,"BMI"])
lasso_1se_BMI$Variable <- "BMI"
colnames(lasso_1se_BMI) <- c("Coefficient", "Variable")
lasso_1se_mRNA <- data.frame(df_coef_lasso_1se[,"mRNA"])
lasso_1se_mRNA$Variable <- "mRNA"
colnames(lasso_1se_mRNA) <- c("Coefficient", "Variable")
lasso_1se_Retransplantation <- data.frame(df_coef_lasso_1se[,"Retransplantation"])
lasso_1se_Retransplantation$Variable <- "Retransplantation"
colnames(lasso_1se_Retransplantation) <- c("Coefficient", "Variable")
lasso_1se_TransplantAge <- data.frame(df_coef_lasso_1se[,"TransplantAge"])
lasso_1se_TransplantAge$Variable <- "TransplantAge"
colnames(lasso_1se_TransplantAge) <- c("Coefficient", "Variable")
lasso_1se_DialysisYears <- data.frame(df_coef_lasso_1se[,"DialysisYears"])
lasso_1se_DialysisYears$Variable <- "DialysisYears"
colnames(lasso_1se_DialysisYears) <- c("Coefficient", "Variable")
lasso_1se_Diabetes <- data.frame(df_coef_lasso_1se[,"Diabetes"])
lasso_1se_Diabetes$Variable <- "Diabetes"
colnames(lasso_1se_Diabetes) <- c("Coefficient", "Variable")
lasso_1se_Steroid <- data.frame(df_coef_lasso_1se[,"Steroid"])
lasso_1se_Steroid$Variable <- "Steroid"
colnames(lasso_1se_Steroid) <- c("Coefficient", "Variable")
lasso_1se_Belatacept <- data.frame(df_coef_lasso_1se[,"Belatacept"])
lasso_1se_Belatacept$Variable <- "Belatacept"
colnames(lasso_1se_Belatacept) <- c("Coefficient", "Variable")
lasso_1se_MPADose <- data.frame(df_coef_lasso_1se[,"MPADose"])
lasso_1se_MPADose$Variable <- "MPADose"
colnames(lasso_1se_MPADose) <- c("Coefficient", "Variable")
lasso_1se_CNI <- data.frame(df_coef_lasso_1se[,"CNI"])
lasso_1se_CNI$Variable <- "CNI"
colnames(lasso_1se_CNI) <- c("Coefficient", "Variable")
lasso_1se_MoreThan2IS <- data.frame(df_coef_lasso_1se[,"MoreThan2IS"])
lasso_1se_MoreThan2IS$Variable <- "MoreThan2IS"
colnames(lasso_1se_MoreThan2IS) <- c("Coefficient", "Variable")
lasso_1se_GFR4vHPSCH <- data.frame(df_coef_lasso_1se[,"GFR4vHPSCH"])
lasso_1se_GFR4vHPSCH$Variable <- "GFR4vHPSCH"
colnames(lasso_1se_GFR4vHPSCH) <- c("Coefficient", "Variable")
lasso_1se_LymphoEB <- data.frame(df_coef_lasso_1se[,"LymphoEB"])
lasso_1se_LymphoEB$Variable <- "LymphoEB"
colnames(lasso_1se_LymphoEB) <- c("Coefficient", "Variable")
lasso_1se_HbEB <- data.frame(df_coef_lasso_1se[,"HbEB"])
lasso_1se_HbEB$Variable <- "HbEB"
colnames(lasso_1se_HbEB) <- c("Coefficient", "Variable")
lasso_1se_Albuminuria <- data.frame(df_coef_lasso_1se[,"Albuminuria"])
lasso_1se_Albuminuria$Variable <- "Albuminuria"
colnames(lasso_1se_Albuminuria) <- c("Coefficient", "Variable")
lasso_1se_TimeSinceLastVacc <- data.frame(df_coef_lasso_1se[,"TimeSinceLastVacc"])
lasso_1se_TimeSinceLastVacc$Variable <- "TimeSinceLastVacc"
colnames(lasso_1se_TimeSinceLastVacc) <- c("Coefficient", "Variable")
lasso_1se_VaccinationNo_3 <- data.frame(df_coef_lasso_1se[,"VaccinationNo_3"])
lasso_1se_VaccinationNo_3$Variable <- "VaccinationNo_3"
colnames(lasso_1se_VaccinationNo_3) <- c("Coefficient", "Variable")

df_coef_lasso_1se_ggplot <- rbind(lasso_1se_IgGpreDoV, lasso_1se_Sex)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_Age)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_BMI)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_mRNA)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_Retransplantation)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_TransplantAge)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_DialysisYears)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_Diabetes)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_Steroid)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_Belatacept)
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_MPADose )
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_CNI )
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_MoreThan2IS )
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_GFR4vHPSCH )
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_LymphoEB )
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_HbEB )
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_Albuminuria )
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_TimeSinceLastVacc )
df_coef_lasso_1se_ggplot <- rbind(df_coef_lasso_1se_ggplot, lasso_1se_VaccinationNo_3 )

View(df_coef_lasso_1se_ggplot)

#Save Coefficients
write.csv(df_coef_lasso_1se_ggplot,"S:\\C13\\MNP\\Digitale-Nephrologie\\vALID\\Impfungen\\Final\\Coefficients_LASSO_1SE_Revision.csv", row.names = FALSE)

#Plot Coefficients
plot_order_coef_1se <- c("IgGpreDoV", "MPADose", "GFR4vHPSCH", "VaccinationNo_3", "Belatacept", "LymphoEB",  "MoreThan2IS", "TransplantAge", "HbEB", "DialysisYears", "BMI", "Retransplantation",  "Sex",  "Age", "Albuminuria", "mRNA", "Diabetes", "TimeSinceLastVacc", "CNI", "Steroid")
df_coef_lasso_1se_ggplot$Variable <- factor(df_coef_lasso_1se_ggplot$Variable, levels = plot_order_coef_1se)
ggplot(data = df_coef_lasso_1se_ggplot, mapping = aes(x=Variable, y=Coefficient)) + 
  geom_point(aes(),alpha=0.8) +
  geom_boxplot(fill="grey",color="black",alpha=0.2) + 
  theme_minimal() +
  coord_flip()+
  scale_x_discrete(limits = rev(levels(df_coef_lasso_1se_ggplot$Variable)))


##### 4. External Validation

#### 4.1 Validation Set 1
#Data Preparation
val_1 <- read.table(file = "E:\\Impfungen_modell\\Final\\Validation_UKD_BO.csv", sep=";", strip.white=TRUE, header=TRUE)

val_1[val_1 == "neg."] <- 0
val_1[val_1=="k.A."] <- NA
val_1$IgGpreDoV[val_1$IgGpreDoV == 0] <- 0
val_1$IgGpreDoV[val_1$IgGpreDoV > 0 & val_1$IgGpreDoV < 35.2] <- 1
val_1$IgGpostDoV[val_1$IgGpostDoV < 35.2] <- 0
val_1$IgGpostDoV[val_1$IgGpostDoV >= 35.2] <- 1
val_1$Albuminuria[val_1$Albuminuria == "<30"] <- 0
val_1$GFR4vHPSCH[val_1$GFR4vHPSCH == "<10"] <- 10
val_1$Sex <- as.numeric(val_1$Sex)
val_1$BMI <- as.numeric(val_1$BMI)
val_1$LymphoEB <- as.numeric(val_1$LymphoEB)
val_1$HbEB <- as.numeric(val_1$HbEB)
val_1$Albuminuria <- as.numeric(val_1$Albuminuria)
val_1$GFR4vHPSCH <- as.numeric(val_1$GFR4vHPSCH)
val_1$DialysisYears <- as.numeric(val_1$DialysisYears)
val_1$IgGpreDoV <- as.numeric(val_1$IgGpreDoV)
val_1$IgGpostDoV <- as.numeric(val_1$IgGpostDoV)
val_1$Albuminuria <- val_1$Albuminuria/1000
val_1$VaccinationNo[val_1$VaccinationNo == 3] <- 1
val_1$VaccinationNo[val_1$VaccinationNo == 4] <- 0
names(val_1)[names(val_1) == "VaccinationNo"] <- "VaccinationNo_3"

# Save Dataset for Multiple Imputation Validation Set 1
val_1_mi <- val_1[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                    'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                    'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_1_mi$IgGpostDoV <- as.character(val_1_mi$IgGpostDoV)

# Complete Case Analysis for LR, LASSO-Min LR, RF, GBRT
val_1_cc_min <- val_1[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                     'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                     'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_1_cc_min <- na.omit(val_1_cc_min)
length(unique(val_1_cc_min$PatientID)) # Show number of unique Patients
val_1_cc_min <- val_1_cc_min[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
summary(val_1_cc_min)
describe(val_1_cc_min)

# Complete Case Analysis for LASSO-1SE LR
val_1_cc_1se <- val_1[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','BMI','Retransplantation','TransplantAge',
                            'DialysisYears', 'Belatacept', 'MPADose', 'MoreThan2IS',
                            'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'IgGpostDoV')]
val_1_cc_1se <- na.omit(val_1_cc_1se)
length(unique(val_1_cc_1se$PatientID)) # Show number of unique Patients
val_1_cc_1se <- merge(val_1[, c('PatientID', 'VaccinationNo_3','Sex', 'Age', 'mRNA', 'Diabetes', 'Steroid', 'CNI', 'Albuminuria', 'TimeSinceLastVacc')], val_1_cc_1se, by=c("PatientID", "VaccinationNo_3"))
val_1_cc_1se <- val_1_cc_1se[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_1_cc_1se[is.na(val_1_cc_1se)] <- 0
summary(val_1_cc_1se)
describe(val_1_cc_1se)

# Mean/Median Imputation
val_1$mRNA[is.na(val_1$mRNA)] <- 1 #die 3 fehlenden waren 4.Impfung -> vermutlich mRNA
val_1$DialysisYears[is.na(val_1$DialysisYears)] <- 3.0 # bei missing dialysis years wird 3.14 imputiert
val_1$BMI[is.na(val_1$BMI)] <- 25.2
val_1$Albuminuria[is.na(val_1$Albuminuria)] <- 0.03
val_1$TimeSinceLastVacc[is.na(val_1$TimeSinceLastVacc)] <- 65
val_1$TransplantAge[is.na(val_1$TransplantAge)] <- 7.8
val_1$Diabetes[is.na(val_1$Diabetes)] <- 0

#Data Preparation
val_1_1 <- val_1[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                   'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                   'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_1_2 <- na.omit(val_1_1)

length(unique(val_1_2$PatientID)) # Show number of unique Patients
summary(val_1_2)
describe(val_1_2)
val_1_2 <- val_1_2[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                    'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                    'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

### Point estimated LASSO-Min and LASSO-1SE
val_1_x_vars <- model.matrix(IgGpostDoV~. , val_1_2)[,-1]
val_1_y_var <- val_1_2$IgGpostDoV
val_1_lasso_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = val_1_x_vars)
val_1_lasso_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = val_1_x_vars)
### ROC-Analysis
val_1_lasso_min_roc <- roc(val_1_y_var ~ as.numeric(val_1_lasso_min_score), plot = TRUE, print.auc = TRUE)
val_1_lasso_1se_roc <- roc(val_1_y_var ~ as.numeric(val_1_lasso_1se_score), plot = TRUE, print.auc = TRUE)
coords(val_1_lasso_min_roc, "best", best.method="closest.topleft", transpose = T)
coords(val_1_lasso_min_roc, cutoff_min, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)
coords(val_1_lasso_1se_roc, cutoff_1se, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)


#### 4.2 Validation Set 2
# Data Preparation
val_2 <- read.table(file = "E:\\Impfungen_modell\\Final\\Validation_Wien_Versuch_2_BO.csv", sep=";", strip.white=TRUE, header=TRUE)
val_2_15 <- val_2 # Save to generate dataset with cutoff 15 U/mL

val_2$BMI <- NA
val_2$Diabetes <- NA
val_2$DialysisYears <- NA
val_2$GFR4vHPSCH <-  NA
val_2$Sex[val_2$Sex == "m"] <- 0
val_2$Sex[val_2$Sex == "w"] <- 1
#val_2$IgGpreDoV[val_2$IgGpreDoV == "<0.4"] <- 0
val_2$IgGpostDoV[val_2$IgGpostDoV < 0.8] <- 0
val_2$IgGpostDoV[val_2$IgGpostDoV >= 0.8] <- 1
val_2$antiHBsSE[val_2$antiHBsSE == ">1000"] <- 1000
val_2$antiHBsSE[val_2$antiHBsSE == "<10"] <- 0
val_2$antiHBsSE[val_2$antiHBsSE == ""] <- NA
val_2$Sex <- as.numeric(val_2$Sex)
val_2$BMI <- as.numeric(val_2$BMI)
val_2$KreatininHP <- as.numeric(val_2$KreatininHP)
val_2$LymphoEB <- as.numeric(val_2$LymphoEB)
val_2$HbEB <- as.numeric(val_2$HbEB)
val_2$Albuminuria <- as.numeric(val_2$Albuminuria)
val_2$DialysisYears <- as.numeric(val_2$DialysisYears)
val_2$IgGpreDoV <- as.numeric(val_2$IgGpreDoV)
val_2$IgGpostDoV <- as.numeric(val_2$IgGpostDoV)
val_2$MPADose <- val_2$MPADose/1000
val_2$Albuminuria <- val_2$Albuminuria/1000
val_2$VaccinationNo[val_2$VaccinationNo == 3] <- 1
val_2$VaccinationNo[val_2$VaccinationNo == 4] <- 0
names(val_2)[names(val_2) == "VaccinationNo"] <- "VaccinationNo_3"

# eGFR from Kreatinin following CKD-EPI Formula
for (i in 1:nrow(val_2)){
  if (val_2$Sex[i]==0){
    if (val_2$KreatininHP[i]/0.9<1){
      val_2$GFR4vHPSCH[i] <- 142*((val_2$KreatininHP[i]/0.9)^-0.302)*0.9938^val_2$Age[i]
    }
    else{
      val_2$GFR4vHPSCH[i] <- 142*((val_2$KreatininHP[i]/0.9)^-1.2)*0.9938^val_2$Age[i]
    }
  }
  else {
    if (val_2$KreatininHP[i]/0.7<1){
      val_2$GFR4vHPSCH[i] <- 142*((val_2$KreatininHP[i]/0.7)^-0.241)*(0.9938^val_2$Age[i])*1.012
    }
    else{
      val_2$GFR4vHPSCH[i] <- 142*((val_2$KreatininHP[i]/0.7)^-1.2)*(0.9938^val_2$Age[i])*1.012
    }
  }
}

# No imputation method since BMI, DialysisYears, Diabetes are completely missing

# Complete Case Analysis for LR, LASSO-Min LR, RF, GBRT
val_2_cc_min <- val_2[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                         'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_2_cc_min <- na.omit(val_2_cc_min)
length(unique(val_2_cc_min$PatientID)) # Show number of unique Patients
val_2_cc_min <- val_2_cc_min[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                         'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
summary(val_2_cc_min)

# Complete Case Analysis for LASSO-1SE LR
val_2_cc_1se <- val_2[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','BMI','Retransplantation','TransplantAge',
                         'DialysisYears', 'Belatacept', 'MPADose', 'MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'IgGpostDoV')]
val_2_cc_1se <- na.omit(val_2_cc_1se)
length(unique(val_2_cc_1se$PatientID)) # Show number of unique Patients
val_2_cc_1se <- merge(val_2[, c('PatientID', 'VaccinationNo_3','Sex', 'Age', 'mRNA', 'Diabetes', 'Steroid', 'CNI', 'Albuminuria', 'TimeSinceLastVacc')], val_2_cc_1se, by=c("PatientID", "VaccinationNo_3"))
val_2_cc_1se <- val_2_cc_1se[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
summary(val_2_cc_1se)


### For both 20 and 12 variable models -> complete case analysis does not work for validation set 2

### Mean/Median Imputation
val_2$DialysisYears[is.na(val_2$DialysisYears)] <- 3.0 
val_2$BMI[is.na(val_2$BMI)] <- 25.2
val_2$Albuminuria[is.na(val_2$Albuminuria)] <- 0.03
val_2$TimeSinceLastVacc[is.na(val_2$TimeSinceLastVacc)] <- 65
val_2$TransplantAge[is.na(val_2$TransplantAge)] <- 7.8
val_2$Diabetes[is.na(val_2$Diabetes)] <- 0


### Select Dataset with 20 predictor variables
val_2_1 <- val_2[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                    'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                    'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_2_2 <- na.omit(val_2_1)
summary(val_2_2)
describe(val_2_2)

### Test LASSO Min und 1SE Model
val_2_x_vars <- model.matrix(IgGpostDoV~. , val_2_2)[,-1]
val_2_y_var <- val_2_2$IgGpostDoV
val_2_lasso_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = val_2_x_vars)
val_2_lasso_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = val_2_x_vars)
# ROC Analysis
val_2_lasso_min_roc <- roc(val_2_y_var ~ as.numeric(val_2_lasso_min_score), plot = TRUE, print.auc = TRUE)
val_2_lasso_1se_roc <- roc(val_2_y_var ~ as.numeric(val_2_lasso_1se_score), plot = TRUE, print.auc = TRUE)
coords(val_2_lasso_min_roc, "best", best.method="closest.topleft", transpose = T)
coords(val_2_lasso_min_roc, cutoff_min, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)
coords(val_2_lasso_1se_roc, cutoff_1se, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)

### Validation Set 2 with Cutoff 15 U/mL for Elecsys (Roche)
val_2_15$BMI <- NA
val_2_15$Diabetes <- NA
val_2_15$DialysisYears <- NA
val_2_15$GFR4vHPSCH <-  NA
val_2_15$Sex[val_2_15$Sex == "m"] <- 0
val_2_15$Sex[val_2_15$Sex == "w"] <- 1
val_2_15$IgGpostDoV[val_2_15$IgGpostDoV < 15] <- 0
val_2_15$IgGpostDoV[val_2_15$IgGpostDoV >= 15] <- 1
val_2_15$Sex <- as.numeric(val_2_15$Sex)
val_2_15$BMI <- as.numeric(val_2_15$BMI)
val_2_15$KreatininHP <- as.numeric(val_2_15$KreatininHP)
val_2_15$LymphoEB <- as.numeric(val_2_15$LymphoEB)
val_2_15$HbEB <- as.numeric(val_2_15$HbEB)
val_2_15$Albuminuria <- as.numeric(val_2_15$Albuminuria)
val_2_15$DialysisYears <- as.numeric(val_2_15$DialysisYears)
val_2_15$IgGpreDoV <- as.numeric(val_2_15$IgGpreDoV)
val_2_15$IgGpostDoV <- as.numeric(val_2_15$IgGpostDoV)
val_2_15$MPADose <- val_2_15$MPADose/1000
val_2_15$Albuminuria <- val_2_15$Albuminuria/1000
val_2_15$VaccinationNo[val_2_15$VaccinationNo == 3] <- 1
val_2_15$VaccinationNo[val_2_15$VaccinationNo == 4] <- 0
names(val_2_15)[names(val_2_15) == "VaccinationNo"] <- "VaccinationNo_3"

for (i in 1:nrow(val_2_15)){
  if (val_2_15$Sex[i]==0){
    if (val_2_15$KreatininHP[i]/0.9<1){
      val_2_15$GFR4vHPSCH[i] <- 142*((val_2_15$KreatininHP[i]/0.9)^-0.302)*0.9938^val_2_15$Age[i]
    }
    else{
      val_2_15$GFR4vHPSCH[i] <- 142*((val_2_15$KreatininHP[i]/0.9)^-1.2)*0.9938^val_2_15$Age[i]
    }
  }
  else {
    if (val_2_15$KreatininHP[i]/0.7<1){
      val_2_15$GFR4vHPSCH[i] <- 142*((val_2_15$KreatininHP[i]/0.7)^-0.241)*(0.9938^val_2_15$Age[i])*1.012
    }
    else{
      val_2_15$GFR4vHPSCH[i] <- 142*((val_2_15$KreatininHP[i]/0.7)^-1.2)*(0.9938^val_2_15$Age[i])*1.012
    }
  }
}

# Complete Case Analysis for LR, LASSO-Min LR, RF, GBRT
val_2_15_cc_min <- val_2_15[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                         'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_2_15_cc_min <- na.omit(val_2_15_cc_min)
length(unique(val_2_15_cc_min$PatientID)) # Show number of unique Patients
val_2_15_cc_min <- val_2_15_cc_min[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

summary(val_2_15_cc_min)

# Complete Case Analysis for LASSO-1SE LR
val_2_15_cc_1se <- val_2_15[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','BMI','Retransplantation','TransplantAge',
                         'DialysisYears', 'Belatacept', 'MPADose', 'MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'IgGpostDoV')]
val_2_15_cc_1se <- na.omit(val_2_15_cc_1se)
length(unique(val_2_15_cc_1se$PatientID)) # Show number of unique Patients
val_2_15_cc_1se <- merge(val_2_15[, c('PatientID', 'VaccinationNo_3','Sex', 'Age', 'mRNA', 'Diabetes', 'Steroid', 'CNI', 'Albuminuria', 'TimeSinceLastVacc')], val_2_15_cc_1se, by=c("PatientID", "VaccinationNo_3"))
val_2_15_cc_1se <- val_2_15_cc_1se[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
summary(val_2_15_cc_1se)

### For both 20 and 12 variable models -> complete case analysis does not work for validation set 2 with cutoff 15 U/mL

### Mean/Median Imputation Validation Set 2 Cutoff 15 U/mL
val_2_15$DialysisYears[is.na(val_2_15$DialysisYears)] <- 3.0 # bei missing dialysis years wird 3.14 imputiert
val_2_15$BMI[is.na(val_2_15$BMI)] <- 25.2
val_2_15$Albuminuria[is.na(val_2_15$Albuminuria)] <- 0.03
val_2_15$TimeSinceLastVacc[is.na(val_2_15$TimeSinceLastVacc)] <- 65
val_2_15$TransplantAge[is.na(val_2_15$TransplantAge)] <- 7.8
val_2_15$Diabetes[is.na(val_2_15$Diabetes)] <- 0

### Select Dataset with 20 predictor variables
val_2_15_1 <- val_2_15[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                       'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                       'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_2_15_2 <- na.omit(val_2_15_1)

View(val_2_15_2)
summary(val_2_15_2)
describe(val_2_15_2)


### Test LASSO Min und 1SE Model
val_2_15_x_vars <- model.matrix(IgGpostDoV~. , val_2_15_2)[,-1]
val_2_15_y_var <- val_2_15_2$IgGpostDoV
val_2_15_lasso_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = val_2_15_x_vars)
val_2_15_lasso_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = val_2_15_x_vars)
### ROC-Analysis
val_2_15_lasso_min_roc <- roc(val_2_15_y_var ~ as.numeric(val_2_15_lasso_min_score), plot = TRUE, print.auc = TRUE)
val_2_15_lasso_1se_roc <- roc(val_2_15_y_var ~ as.numeric(val_2_15_lasso_1se_score), plot = TRUE, print.auc = TRUE)
coords(val_2_15_lasso_min_roc, cutoff_min, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)
coords(val_2_15_lasso_1se_roc, cutoff_1se, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)

#### Validation Set 3
val_3 <- read.table(file = "E:\\Impfungen_modell\\Final\\Validation_Strasbourg_BO_2.csv", sep=";", strip.white=TRUE, header=TRUE)

#Data Preparation
val_3$Sex[val_3$Sex == "m"] <- 0
val_3$Sex[val_3$Sex == "w"] <- 1
val_3$IgGpreDoV[val_3$IgGpreDoV != 0] <- 1
val_3$IgGpostDoV[val_3$IgGpostDoV < 7] <- 0
val_3$IgGpostDoV[val_3$IgGpostDoV >= 7] <- 1
val_3$Sex <- as.numeric(val_3$Sex)
val_3$BMI <- as.numeric(val_3$BMI)
val_3$GFR4vHPSCH <- as.numeric(val_3$GFR4vHPSCH)
val_3$LymphoEB <- as.numeric(val_3$LymphoEB)
val_3$HbEB <- as.numeric(val_3$HbEB)
val_3$DialysisYears <- as.numeric(val_3$DialysisYears)
val_3$IgGpreDoV <- as.numeric(val_3$IgGpreDoV)
val_3$IgGpostDoV <- as.numeric(val_3$IgGpostDoV)
val_3$MPADose <- val_3$MPADose/1000
val_3$UACR <- val_3$UACR*1000/113
val_3$UPCR <- val_3$UPCR*1000/113
val_3$VaccinationNo[val_3$VaccinationNo == 3] <- 1
val_3$VaccinationNo[val_3$VaccinationNo == 4] <- 0
names(val_3)[names(val_3) == "VaccinationNo"] <- "VaccinationNo_3"
# Convert UPCR in UACR
for (i in 1:nrow(val_3)){
  if (is.na(val_3$UACR[i]) & !is.na(val_3$UPCR[i])){
    val_3$UACR[i] <- exp(0.2445*log(min(c(val_3$UPCR[i]/50, 1))) + 1.5531*log(max(c(min(c(val_3$UPCR[i]/500, 1)), 0.1))) + 1.1057*log(max(c(val_3$UPCR[i]/500, 1))) + 5.2562 - 0.0793*max(c(val_3$Sex[i], 0)) + 0.0802 * max(c(val_3$Diabetes[i], 0)) + 0.1339)
  }
}
val_3$Albuminuria <- as.numeric(val_3$UACR/1000)

nrow(val_3)

### Save Validation Set 3 for Multiple Imputation
val_3_mi <- val_3[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                     'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                     'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

# Complete Case Analysis for LR, LASSO-Min LR, RF, GBRT
val_3_cc_min <- val_3[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                         'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_3_cc_min <- na.omit(val_3_cc_min)
length(unique(val_3_cc_min$PatientID)) # Show number of unique Patients
val_3_cc_min <- val_3_cc_min[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
summary(val_3_cc_min)
describe(val_3_cc_min)

# Complete Case Analysis for LASSO-1SE LR
val_3_cc_1se <- val_3[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','BMI','Retransplantation','TransplantAge',
                         'DialysisYears', 'Belatacept', 'MPADose', 'MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'IgGpostDoV')]
val_3_cc_1se <- na.omit(val_3_cc_1se)
nrow(val_3_cc_1se)
length(unique(val_3_cc_1se$PatientID)) # Show number of unique Patients

### Same Patients as for complete variable set
val_3_cc_1se <- val_3_cc_min


## Median/Mean Imputation
val_3$DialysisYears[is.na(val_3$DialysisYears)] <- 3.0 # for missing dialysis years 3.0 years
val_3$BMI[is.na(val_3$BMI)] <- 25.2
val_3$Albuminuria[is.na(val_3$Albuminuria)] <- 0.03
val_3$TimeSinceLastVacc[is.na(val_3$TimeSinceLastVacc)] <- 65
val_3$TransplantAge[is.na(val_3$TransplantAge)] <- 7.8
val_3$Diabetes[is.na(val_3$Diabetes)] <- 0

### Select Dataset with 20 predictor variables
val_3_1 <- val_3[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                    'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                    'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_3_2 <- na.omit(val_3_1)

length(unique(val_3_2$PatientID))
summary(val_3_2)
describe(val_3_2)

val_3_2 <- val_3_2[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                    'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                    'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]


### Test LASSO Min und 1SE Model
val_3_x_vars <- model.matrix(IgGpostDoV~. , val_3_2)[,-1]
val_3_y_var <- val_3_2$IgGpostDoV
val_3_lasso_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = val_3_x_vars)
val_3_lasso_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = val_3_x_vars)
### ROC-Analysis
val_3_lasso_min_roc <- roc(val_3_y_var ~ as.numeric(val_3_lasso_min_score), plot = TRUE, print.auc = TRUE)
val_3_lasso_1se_roc <- roc(val_3_y_var ~ as.numeric(val_3_lasso_1se_score), plot = TRUE, print.auc = TRUE)
coords(val_3_lasso_min_roc, cutoff_min, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)
coords(val_3_lasso_1se_roc, cutoff_1se, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)

#### Validation Set 4
val_4 <- read.table(file = "E:\\Impfungen_modell\\Final\\Validation_Nantes_BO_final_manual.csv", sep=";", strip.white=TRUE, header=TRUE)
val_4$IgGpreDoV[val_4$IgGpreDoV == "<0,4"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoV == "<0,5"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoV == "<1"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoV == "<4"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoV == "<4,81"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoV == "<4,9"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoV == "<7,8"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoV == "<21"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoV == "<50"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoV == ">10"] <- 10
val_4$IgGpreDoV[val_4$IgGpreDoV == ">150"] <- 264
val_4$IgGpreDoV[val_4$IgGpreDoV == ">250"] <- 264
val_4$IgGpreDoV[val_4$IgGpreDoV == ">400"] <- 400
val_4$IgGpreDoV <- as.numeric(val_4$IgGpreDoV)
val_4$IgGpostDoV[val_4$IgGpostDoV == "<0.4"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoV == "<0.5"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoV == "<1"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoV == "<4"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoV == "<4.81"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoV == "<4.9"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoV == "<6.8"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoV == "<7.8"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoV == "<21"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoV == ">10"] <- 10
val_4$IgGpostDoV[val_4$IgGpostDoV == ">150"] <- 150
val_4$IgGpostDoV[val_4$IgGpostDoV == ">250"] <- 264
val_4$IgGpostDoV[val_4$IgGpostDoV == ">750"] <- 750
val_4$IgGpostDoV[val_4$IgGpostDoV == ">2080"] <- 2080
val_4$IgGpostDoV[val_4$IgGpostDoV == ">2500"] <- 2500
val_4$IgGpostDoV[val_4$IgGpostDoV == ">3270"] <- 3270
val_4$IgGpostDoV[val_4$IgGpostDoV == ">12500"] <- 12500
val_4$IgGpostDoV <- as.numeric(val_4$IgGpostDoV)

val_4_15 <- val_4


for (i in 1:nrow(val_4_15)){
  if (!is.na(val_4_15$IgGpreDoV[i]) & !is.na(val_4_15$CutoffpreDoV[i])){
    if (val_4_15$AssaypreDoV[i] == "ECLIA (Roche)" & val_4_15$IgGpreDoV[i] >= val_4_15$CutoffpreDoV[i] & val_4_15$IgGpreDoV[i] < 15) {
      val_4_15$IgGpreDoVint[i] <- "low pos"
    }
    else if (val_4_15$AssaypreDoV[i] == "ECLIA (Roche)" & val_4_15$IgGpreDoV[i] >= 15) {
      val_4_15$IgGpreDoVint[i] <- "pos"
    }
  }
  if (!is.na(val_4_15$IgGpostDoV[i]) & !is.na(val_4_15$CutoffpostDoV[i])){
    if (val_4_15$Assaypost[i] == "ECLIA (Roche)" & val_4_15$IgGpostDoV[i] >= val_4_15$CutoffpostDoV[i] & val_4_15$IgGpostDoV[i] < 15) {
      val_4_15$IgGpostDoVint[i] <- "low pos"
    }
    else if (val_4_15$Assaypost[i] == "ECLIA (Roche)" & val_4_15$IgGpostDoV[i] >= 15) {
      val_4_15$IgGpostDoVint[i] <- "pos"
    }
  }
}

#### Applying Exclusion Criteria and counting how often which exclusion criterion applies
nrow(val_4)
val_4_3 <- subset (val_4, val_4$VaccinationNo==3)
nrow(val_4_3)
val_4_4 <- subset (val_4, val_4$VaccinationNo==4)
nrow(val_4_4)
val_4_3_pos <- subset(val_4, val_4$IgGpreDoVint=="pos" & val_4$VaccinationNo==3)
nrow(val_4_3_pos)
val_4_4_pos <- subset(val_4, val_4$IgGpreDoVint=="pos" & val_4$VaccinationNo==4)
nrow(val_4_4_pos)
val_4_3_noassay <- subset(val_4, (is.na(val_4$CutoffpreDoV) | is.na(val_4$CutoffpostDoV)) & val_4$IgGpreDoVint!="pos" & val_4$VaccinationNo==3)
nrow(val_4_3_noassay)
val_4_4_noassay <- subset(val_4, (is.na(val_4$CutoffpreDoV) | is.na(val_4$CutoffpostDoV)) & val_4$IgGpreDoVint!="pos" & val_4$VaccinationNo==4)
nrow(val_4_4_noassay)

nrow(val_4_15)
val_4_15_3 <- subset (val_4_15, val_4_15$VaccinationNo==3)
nrow(val_4_15_3)
val_4_15_4 <- subset (val_4_15, val_4_15$VaccinationNo==4)
nrow(val_4_15_4)
val_4_15_3_pos <- subset(val_4_15, val_4_15$IgGpreDoVint=="pos" & val_4_15$VaccinationNo==3)
nrow(val_4_15_3_pos)
val_4_15_4_pos <- subset(val_4_15, val_4_15$IgGpreDoVint=="pos" & val_4_15$VaccinationNo==4)
nrow(val_4_15_4_pos)
val_4_15_3_noassay <- subset(val_4_15, (is.na(val_4_15$CutoffpreDoV) | is.na(val_4_15$CutoffpostDoV)) & (val_4_15$IgGpreDoVint=="neg" | val_4_15$IgGpreDoVint=="low pos") & val_4_15$VaccinationNo==3)
nrow(val_4_15_3_noassay)
val_4_15_4_noassay <- subset(val_4_15, (is.na(val_4_15$CutoffpreDoV) | is.na(val_4_15$CutoffpostDoV)) & val_4_15$IgGpreDoVint!="pos" & val_4_15$VaccinationNo==4)
nrow(val_4_15_4_noassay)

val_4 <- subset (val_4, !is.na(val_4$CutoffpreDoV) & !is.na(val_4$CutoffpostDoV))
nrow(val_4)
val_4 <- subset (val_4, (val_4$IgGpreDoVint=="neg" | val_4$IgGpreDoVint=="low pos")  )
nrow(val_4)
val_4_15 <- subset (val_4_15, (val_4_15$IgGpreDoVint=="neg" | val_4_15$IgGpreDoVint=="low pos") & !is.na(val_4_15$CutoffpreDoV) & !is.na(val_4_15$CutoffpostDoV))
nrow(val_4_15)


### Data Preparation
val_4$Sex[val_4$Sex == "M"] <- 0
val_4$Sex[val_4$Sex == "F"] <- 1
val_4$IgGpreDoV[val_4$IgGpreDoVint == "neg"] <- 0
val_4$IgGpreDoV[val_4$IgGpreDoVint == "low pos"] <- 1
val_4$IgGpostDoV[val_4$IgGpostDoVint == "neg"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoVint == "low pos"] <- 0
val_4$IgGpostDoV[val_4$IgGpostDoVint == "pos"] <- 1
val_4$Sex <- as.numeric(val_4$Sex)
val_4$BMI <- as.numeric(val_4$BMI)
val_4$GFR4vHPSCH <- as.numeric(val_4$GFR4vHPSCH)
val_4$LymphoEB <- as.numeric(val_4$LymphoEB)
val_4$HbEB <- as.numeric(val_4$HbEB)
val_4$DialysisYears <- as.numeric(val_4$DialysisYears)
val_4$IgGpreDoV <- as.numeric(val_4$IgGpreDoV)
val_4$IgGpostDoV <- as.numeric(val_4$IgGpostDoV)
val_4$MPADose <- val_4$MPADose/1000
val_4$Proteinuria <- val_4$Proteinuria*1000
val_4$VaccinationNo[val_4$VaccinationNo == 3] <- 1
val_4$VaccinationNo[val_4$VaccinationNo == 4] <- 0
names(val_4)[names(val_4) == "VaccinationNo"] <- "VaccinationNo_3"
# Convert Proteinuria to Albuminuria
for (i in 1:nrow(val_4)){
  if (is.na(val_4$Albuminuria[i]) & !is.na(val_4$Proteinuria[i])){
    val_4$Albuminuria[i] <- exp(0.2445*log(min(c(val_4$Proteinuria[i]/50, 1))) + 1.5531*log(max(c(min(c(val_4$Proteinuria[i]/500, 1)), 0.1))) + 1.1057*log(max(c(val_4$Proteinuria[i]/500, 1))) + 5.2562 - 0.0793*max(c(val_4$Sex[i], 0)) + 0.0802 * max(c(val_4$Diabetes[i], 0)) + 0.1339)
  }
}
val_4$Albuminuria <- as.numeric(val_4$Albuminuria)/1000

# Apply Additional Exclusion Criteria
nrow(val_4)
val_4 <- subset(val_4, (!is.na(val_4$GFR4vHPSCH) & val_4$VaccinationNo_3==1) | val_4$VaccinationNo_3==0)
nrow(val_4)
val_4 <- subset(val_4, (!is.na(val_4$LymphoEB) & val_4$VaccinationNo_3==1) | val_4$VaccinationNo_3==0)
nrow(val_4)
val_4 <- subset(val_4, (!is.na(val_4$HbEB) & val_4$VaccinationNo_3==1) | val_4$VaccinationNo_3==0)
nrow(val_4)

val_4 <- subset(val_4, (!is.na(val_4$GFR4vHPSCH) & val_4$VaccinationNo_3==0) | val_4$VaccinationNo_3==1)
nrow(val_4)
val_4 <- subset(val_4, (!is.na(val_4$LymphoEB) & val_4$VaccinationNo_3==0) | val_4$VaccinationNo_3==1)
nrow(val_4)
val_4 <- subset(val_4, (!is.na(val_4$HbEB) & val_4$VaccinationNo_3==0) | val_4$VaccinationNo_3==1)
nrow(val_4)

# Save Dataset for multiple imputation
val_4_mi <- val_4[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                     'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                     'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

# Complete Case Analysis for LR, LASSO-Min LR, RF, GBRT
val_4_cc_min <- val_4[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                         'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_4_cc_min <- na.omit(val_4_cc_min)
nrow(val_4_cc_min)
length(unique(val_4_cc_min$PatientID)) # Show number of unique Patients

val_4_cc_min <- val_4_cc_min[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
summary(val_4_cc_min)
describe(val_4_cc_min)

# Complete Case Analysis for LASSO-1SE LR
val_4_cc_1se <- val_4[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','BMI','Retransplantation','TransplantAge',
                         'DialysisYears', 'Belatacept', 'MPADose', 'MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'IgGpostDoV')]
val_4_cc_1se <- na.omit(val_4_cc_1se)
nrow(val_4_cc_1se)
length(unique(val_4_cc_1se$PatientID)) # Show number of unique Patients
val_4_cc_1se <- merge(val_4[, c('PatientID', 'VaccinationNo_3','Sex', 'Age', 'mRNA', 'Diabetes', 'Steroid', 'CNI', 'Albuminuria', 'TimeSinceLastVacc')], val_4_cc_1se, by=c("PatientID", "VaccinationNo_3"))
val_4_cc_1se <- val_4_cc_1se[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_4_cc_1se[is.na(val_4_cc_1se)] <- 0
summary(val_4_cc_1se)
describe(val_4_cc_1se)


### Mean Median Imputation Validation Set 4 
val_4$DialysisYears[is.na(val_4$DialysisYears)] <- 3.0
val_4$BMI[is.na(val_4$BMI)] <- 25.2
val_4$Albuminuria[is.na(val_4$Albuminuria)] <- 0.03
val_4$TimeSinceLastVacc[is.na(val_4$TimeSinceLastVacc)] <- 65
val_4$TransplantAge[is.na(val_4$TransplantAge)] <- 7.8
val_4$Diabetes[is.na(val_4$Diabetes)] <- 0

### Select Dataset with 20 predictor variables
val_4_1 <- val_4[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                    'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                    'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_4_2 <- na.omit(val_4_1)
length(unique(val_4_2$PatientID))
summary(val_4_2)
describe(val_4_2)

## Count Events for Figure S4
val_4_3_pos <- subset(val_4, val_4$VaccinationNo_3==1 & val_4$IgGpostDoV==1)
nrow(val_4_3_pos)
val_4_3_neg <- subset(val_4, val_4$VaccinationNo_3==1 & val_4$IgGpostDoV==0)
nrow(val_4_3_neg)
val_4_4_pos <- subset(val_4, val_4$VaccinationNo_3==0 & val_4$IgGpostDoV==1)
nrow(val_4_4_pos)
val_4_4_neg <- subset(val_4, val_4$VaccinationNo_3==0 & val_4$IgGpostDoV==0)
nrow(val_4_4_neg)

val_4_2 <- val_4_2[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                    'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                    'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

### LASSO Min und 1SE Model testen
val_4_x_vars <- model.matrix(IgGpostDoV~. , val_4_2)[,-1]
val_4_y_var <- val_4_2$IgGpostDoV
val_4_lasso_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = val_4_x_vars)
val_4_lasso_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = val_4_x_vars)
### ROC-Analyse auf dem kompletten Datensatz -> überschätzt Performance, ist aber für Cutoffbestimmung geeignet
val_4_lasso_min_roc <- roc(val_4_y_var ~ as.numeric(val_4_lasso_min_score), plot = TRUE, print.auc = TRUE)
val_4_lasso_1se_roc <- roc(val_4_y_var ~ as.numeric(val_4_lasso_1se_score), plot = TRUE, print.auc = TRUE)
coords(val_4_lasso_min_roc, cutoff_min, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)
coords(val_4_lasso_1se_roc, cutoff_1se, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)


### Data Preparation Validation Set 4 Cutoff 15 U/mL for Elecsys assay
val_4_15$Sex[val_4_15$Sex == "M"] <- 0
val_4_15$Sex[val_4_15$Sex == "F"] <- 1
val_4_15$IgGpreDoV[val_4_15$IgGpreDoVint == "neg"] <- 0
val_4_15$IgGpreDoV[val_4_15$IgGpreDoVint == "low pos"] <- 1
val_4_15$IgGpostDoV[val_4_15$IgGpostDoVint == "neg"] <- 0
val_4_15$IgGpostDoV[val_4_15$IgGpostDoVint == "low pos"] <- 0
val_4_15$IgGpostDoV[val_4_15$IgGpostDoVint == "pos"] <- 1
val_4_15$Sex <- as.numeric(val_4_15$Sex)
val_4_15$BMI <- as.numeric(val_4_15$BMI)
val_4_15$GFR4vHPSCH <- as.numeric(val_4_15$GFR4vHPSCH)
val_4_15$LymphoEB <- as.numeric(val_4_15$LymphoEB)
val_4_15$HbEB <- as.numeric(val_4_15$HbEB)
val_4_15$DialysisYears <- as.numeric(val_4_15$DialysisYears)
val_4_15$IgGpreDoV <- as.numeric(val_4_15$IgGpreDoV)
val_4_15$IgGpostDoV <- as.numeric(val_4_15$IgGpostDoV)
val_4_15$MPADose <- val_4_15$MPADose/1000
val_4_15$Proteinuria <- val_4_15$Proteinuria*1000
val_4_15$VaccinationNo[val_4_15$VaccinationNo == 3] <- 1
val_4_15$VaccinationNo[val_4_15$VaccinationNo == 4] <- 0
names(val_4_15)[names(val_4_15) == "VaccinationNo"] <- "VaccinationNo_3"
# Albuminuria from Proteinuria
for (i in 1:nrow(val_4_15)){
  if (is.na(val_4_15$Albuminuria[i]) & !is.na(val_4_15$Proteinuria[i])){
    val_4_15$Albuminuria[i] <- exp(0.2445*log(min(c(val_4_15$Proteinuria[i]/50, 1))) + 1.5531*log(max(c(min(c(val_4_15$Proteinuria[i]/500, 1)), 0.1))) + 1.1057*log(max(c(val_4_15$Proteinuria[i]/500, 1))) + 5.2562 - 0.0793*max(c(val_4_15$Sex[i], 0)) + 0.0802 * max(c(val_4_15$Diabetes[i], 0)) + 0.1339)
  }
}
val_4_15$Albuminuria <- as.numeric(val_4_15$Albuminuria)/1000

nrow(val_4_15)
val_4_15 <- subset(val_4_15, (!is.na(val_4_15$GFR4vHPSCH) & val_4_15$VaccinationNo_3==1) | val_4_15$VaccinationNo_3==0)
nrow(val_4_15)
val_4_15 <- subset(val_4_15, (!is.na(val_4_15$LymphoEB) & val_4_15$VaccinationNo_3==1) | val_4_15$VaccinationNo_3==0)
nrow(val_4_15)
val_4_15 <- subset(val_4_15, (!is.na(val_4_15$HbEB) & val_4_15$VaccinationNo_3==1) | val_4_15$VaccinationNo_3==0)
nrow(val_4_15)

val_4_15 <- subset(val_4_15, (!is.na(val_4_15$GFR4vHPSCH) & val_4_15$VaccinationNo_3==0) | val_4_15$VaccinationNo_3==1)
nrow(val_4_15)
val_4_15 <- subset(val_4_15, (!is.na(val_4_15$LymphoEB) & val_4_15$VaccinationNo_3==0) | val_4_15$VaccinationNo_3==1)
nrow(val_4_15)
val_4_15 <- subset(val_4_15, (!is.na(val_4_15$HbEB) & val_4_15$VaccinationNo_3==0) | val_4_15$VaccinationNo_3==1)
nrow(val_4_15)

# Save dataset for multiple imputation
val_4_15_mi <- val_4_15[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                           'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                           'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

# Complete Case Analysis for LR, LASSO-Min LR, RF, GBRT
val_4_15_cc_min <- val_4_15[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                         'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_4_15_cc_min <- na.omit(val_4_15_cc_min)
length(unique(val_4_15_cc_min$PatientID)) # Show number of unique Patients

val_4_15_cc_min <- val_4_15_cc_min[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
summary(val_4_15_cc_min)
describe(val_4_15_cc_min)

# Complete Case Analysis for LASSO-1SE LR
val_4_15_cc_1se <- val_4_15[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','BMI','Retransplantation','TransplantAge',
                         'DialysisYears', 'Belatacept', 'MPADose', 'MoreThan2IS',
                         'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'IgGpostDoV')]
val_4_15_cc_1se <- na.omit(val_4_15_cc_1se)
nrow(val_4_15_cc_1se)
length(unique(val_4_15_cc_1se$PatientID)) # Show number of unique Patients
val_4_15_cc_1se <- merge(val_4_15[, c('PatientID', 'VaccinationNo_3','Sex', 'Age', 'mRNA', 'Diabetes', 'Steroid', 'CNI', 'Albuminuria', 'TimeSinceLastVacc')], val_4_15_cc_1se, by=c("PatientID", "VaccinationNo_3"))
val_4_15_cc_1se <- val_4_15_cc_1se[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                                'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                                'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

val_4_15_cc_1se[is.na(val_4_15_cc_1se)] <- 0
summary(val_4_15_cc_1se)
describe(val_4_15_cc_1se)


### Median/Mean Imputation Validation Set 4 Cutoff 15 U/mL
val_4_15$DialysisYears[is.na(val_4_15$DialysisYears)] <- 3.0 # Median
val_4_15$BMI[is.na(val_4_15$BMI)] <- 25.2 # Mean
val_4_15$Albuminuria[is.na(val_4_15$Albuminuria)] <- 0.03 # Median
val_4_15$TimeSinceLastVacc[is.na(val_4_15$TimeSinceLastVacc)] <- 65 # Median
val_4_15$TransplantAge[is.na(val_4_15$TransplantAge)] <- 7.8 # Median
val_4_15$Diabetes[is.na(val_4_15$Diabetes)] <- 0 # Median

### Select dataset with 20 predictor variables
val_4_15_1 <- val_4_15[,c('PatientID','IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                       'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                       'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]
val_4_15_2 <- na.omit(val_4_15_1)
length(unique(val_4_15_2$PatientID))
summary(val_4_15_2)
describe(val_4_15_2)

val_4_15_3_pos <- subset(val_4_15, val_4_15$VaccinationNo_3==1 & val_4_15$IgGpostDoV==1)
nrow(val_4_15_3_pos)
val_4_15_3_neg <- subset(val_4_15, val_4_15$VaccinationNo_3==1 & val_4_15$IgGpostDoV==0)
nrow(val_4_15_3_neg)
val_4_15_4_pos <- subset(val_4_15, val_4_15$VaccinationNo_3==0 & val_4_15$IgGpostDoV==1)
nrow(val_4_15_4_pos)
val_4_15_4_neg <- subset(val_4_15, val_4_15$VaccinationNo_3==0 & val_4_15$IgGpostDoV==0)
nrow(val_4_15_4_neg)

val_4_15_2 <- val_4_15_2[,c('IgGpreDoV', 'VaccinationNo_3','Sex','Age','BMI','mRNA','Retransplantation','TransplantAge',
                       'DialysisYears', 'Diabetes', 'Steroid','Belatacept', 'MPADose', 'CNI','MoreThan2IS',
                       'GFR4vHPSCH', 'LymphoEB', 'HbEB', 'Albuminuria','TimeSinceLastVacc','IgGpostDoV')]

### LASSO Min und 1SE Model testen
val_4_15_x_vars <- model.matrix(IgGpostDoV~. , val_4_15_2)[,-1]
val_4_15_y_var <- val_4_15_2$IgGpostDoV
val_4_15_lasso_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = val_4_15_x_vars)
val_4_15_lasso_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = val_4_15_x_vars)
### ROC-Analyse auf dem kompletten Datensatz -> überschätzt Performance, ist aber für Cutoffbestimmung geeignet
val_4_15_lasso_min_roc <- roc(val_4_15_y_var ~ as.numeric(val_4_15_lasso_min_score), plot = TRUE, print.auc = TRUE)
val_4_15_lasso_1se_roc <- roc(val_4_15_y_var ~ as.numeric(val_4_15_lasso_1se_score), plot = TRUE, print.auc = TRUE)
coords(val_4_15_lasso_min_roc, cutoff_min, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)
coords(val_4_15_lasso_1se_roc, cutoff_1se, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)

##### Pooled Validation Sets with Cutoff 0.8
val_all <- rbind(val_1_2, val_2_2)
val_all <- rbind(val_all, val_3_2)
val_all <- rbind(val_all, val_4_2)

val_all_x_vars <- model.matrix(IgGpostDoV~. , val_all)[,-1]
val_all_y_var <- val_all$IgGpostDoV
val_all_lasso_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = val_all_x_vars)
val_all_lasso_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = val_all_x_vars)
### ROC-Analyis
val_all_lasso_min_roc <- roc(val_all_y_var ~ as.numeric(val_all_lasso_min_score), plot = TRUE, print.auc = TRUE)
val_all_lasso_1se_roc <- roc(val_all_y_var ~ as.numeric(val_all_lasso_1se_score), plot = TRUE, print.auc = TRUE)
coords(val_all_lasso_min_roc, cutoff_min, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)
coords(val_all_lasso_1se_roc, cutoff_1se, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)

##### Multiple Imputation Validation with Cutoff 0.8
val_all_mi <- rbind(val_1_mi, val_2_mi)
val_all_mi <- rbind(val_all_mi, val_3_mi)
val_all_mi <- rbind(val_all_mi, val_4_mi)

val_all_mi$IgGpostDoV <- as.character(val_all_mi$IgGpostDoV)

imp_val_all <- mice(val_all_mi, m=5, maxit = 50, method = 'pmm', seed = 0)
summary(imp_val_all)
val_all_mi_1 <- complete(imp_val_all,1)
val_all_mi_2 <- complete(imp_val_all,2)
val_all_mi_3 <- complete(imp_val_all,3)
val_all_mi_4 <- complete(imp_val_all,4)
val_all_mi_5 <- complete(imp_val_all,5)

##### Complete Validation with Cutoff 15
val_all_15 <- rbind(val_1_2, val_2_15_2)
val_all_15 <- rbind(val_all_15, val_3_2)
val_all_15 <- rbind(val_all_15, val_4_15_2)

val_all_15_x_vars <- model.matrix(IgGpostDoV~. , val_all_15)[,-1]
val_all_15_y_var <- val_all_15$IgGpostDoV
val_all_15_lasso_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = val_all_15_x_vars)
val_all_15_lasso_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = val_all_15_x_vars)
### ROC-Analysis
val_all_15_lasso_min_roc <- roc(val_all_15_y_var ~ as.numeric(val_all_15_lasso_min_score), plot = TRUE, print.auc = TRUE)
val_all_15_lasso_1se_roc <- roc(val_all_15_y_var ~ as.numeric(val_all_15_lasso_1se_score), plot = TRUE, print.auc = TRUE)
coords(val_all_15_lasso_min_roc, "best", best.method="closest.topleft", transpose = T)
coords(val_all_15_lasso_min_roc, cutoff_min, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)
coords(val_all_15_lasso_1se_roc, cutoff_1se, input="threshold",  ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = F)

##### Multiple Imputation Validation with Cutoff 15
val_all_15_mi <- rbind(val_1_mi, val_2_15_mi)
val_all_15_mi <- rbind(val_all_15_mi, val_3_mi)
val_all_15_mi <- rbind(val_all_15_mi, val_4_15_mi)

val_all_15_mi$IgGpostDoV <- as.character(val_all_15_mi$IgGpostDoV)

imp_val_all_15 <- mice(val_all_15_mi, m=5, maxit = 50, method = 'pmm', seed = 0)
summary(imp_val_all_15)
val_all_15_mi_1 <- complete(imp_val_all_15,1)
val_all_15_mi_2 <- complete(imp_val_all_15,2)
val_all_15_mi_3 <- complete(imp_val_all_15,3)
val_all_15_mi_4 <- complete(imp_val_all_15,4)
val_all_15_mi_5 <- complete(imp_val_all_15,5)


#### Bootstrapped External Validation

df_sample_val_1_min <- data.frame()
df_sample_val_1_1se <- data.frame()
df_sample_val_1_cc_min <- data.frame()
df_sample_val_1_cc_1se <- data.frame()
df_sample_val_2_min <- data.frame()
df_sample_val_2_1se <- data.frame()
df_sample_val_2_15_min <- data.frame()
df_sample_val_2_15_1se <- data.frame()
df_sample_val_3_min <- data.frame()
df_sample_val_3_1se <- data.frame()
df_sample_val_3_cc_min <- data.frame()
df_sample_val_3_cc_1se <- data.frame()
df_sample_val_4_min <- data.frame()
df_sample_val_4_1se <- data.frame()
df_sample_val_4_cc_min <- data.frame()
df_sample_val_4_cc_1se <- data.frame()
df_sample_val_4_15_min <- data.frame()
df_sample_val_4_15_1se <- data.frame()
df_sample_val_4_15_cc_min <- data.frame()
df_sample_val_4_15_cc_1se <- data.frame()
df_sample_val_all_min <- data.frame()
df_sample_val_all_1se <- data.frame()
df_sample_val_all_mi_min <- data.frame()
df_sample_val_all_mi_1se <- data.frame()
df_sample_val_all_15_min <- data.frame()
df_sample_val_all_15_1se <- data.frame()
df_sample_val_all_15_mi_min <- data.frame()
df_sample_val_all_15_mi_1se <- data.frame()
df_sample_val_all_lr <- data.frame()
df_sample_val_all_rf <- data.frame()
df_sample_val_all_gbm <- data.frame()
df_sample_val_all_15_lr <- data.frame()
df_sample_val_all_15_rf <- data.frame()
df_sample_val_all_15_gbm <- data.frame()


### 1000-fold bootstrapping
for (i in 1:1000) {
  set.seed(i)
  
  # Validation Set 1 Mean/Median imputed
  sample_val_1 <- val_1_2[sample(1:nrow(val_1_2), nrow(val_1_2), replace = TRUE), ]
  sample_val_1_x_vars <- model.matrix(IgGpostDoV~. , sample_val_1)[,-1]
  sample_val_1_y_var <- sample_val_1$IgGpostDoV
  
  sample_val_1_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_1_x_vars)
  sample_val_1_min_roc <- roc(sample_val_1_y_var ~ as.numeric(sample_val_1_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_1_metric_min <- coords(sample_val_1_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_1_min_auc <- as.numeric(sample_val_1_min_roc$auc)
  sample_val_1_metric_min[6]<-sample_val_1_min_auc
  df_sample_val_1_min <- rbind(df_sample_val_1_min, sample_val_1_metric_min)
  
  sample_val_1_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_1_x_vars)
  sample_val_1_1se_roc <- roc(sample_val_1_y_var ~ as.numeric(sample_val_1_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_1_metric_1se <- coords(sample_val_1_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_1_1se_auc <- as.numeric(sample_val_1_1se_roc$auc)
  sample_val_1_metric_1se[6]<-sample_val_1_1se_auc
  df_sample_val_1_1se <- rbind(df_sample_val_1_1se, sample_val_1_metric_1se)
  
  # Validation Set 1 Complete Case analysis LASSO Min
  sample_val_1_cc_min <- val_1_cc_min[sample(1:nrow(val_1_cc_min), nrow(val_1_cc_min), replace = TRUE), ]
  sample_val_1_cc_min_x_vars <- model.matrix(IgGpostDoV~. , sample_val_1_cc_min)[,-1]
  sample_val_1_cc_min_y_var <- sample_val_1_cc_min$IgGpostDoV
  
  sample_val_1_cc_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_1_cc_min_x_vars)
  sample_val_1_cc_min_roc <- roc(sample_val_1_cc_min_y_var ~ as.numeric(sample_val_1_cc_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_1_cc_min_metric <- coords(sample_val_1_cc_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_1_cc_min_auc <- as.numeric(sample_val_1_cc_min_roc$auc)
  sample_val_1_cc_min_metric[6]<-sample_val_1_cc_min_auc
  df_sample_val_1_cc_min <- rbind(df_sample_val_1_cc_min, sample_val_1_cc_min_metric)
  
  # Validation Set 1 Complete Case analysis LASSO 1SE
  sample_val_1_cc_1se <- val_1_cc_1se[sample(1:nrow(val_1_cc_1se), nrow(val_1_cc_1se), replace = TRUE), ]
  sample_val_1_cc_1se_x_vars <- model.matrix(IgGpostDoV~. , sample_val_1_cc_1se)[,-1]
  sample_val_1_cc_1se_y_var <- sample_val_1_cc_1se$IgGpostDoV
  
  sample_val_1_cc_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_1_cc_1se_x_vars)
  sample_val_1_cc_1se_roc <- roc(sample_val_1_cc_1se_y_var ~ as.numeric(sample_val_1_cc_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_1_cc_1se_metric <- coords(sample_val_1_cc_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_1_cc_1se_auc <- as.numeric(sample_val_1_cc_1se_roc$auc)
  sample_val_1_cc_1se_metric[6]<-sample_val_1_cc_1se_auc
  df_sample_val_1_cc_1se <- rbind(df_sample_val_1_cc_1se, sample_val_1_cc_1se_metric)
 
  # Validation Set 2 Mean/Median imputed Cutoff 0.8 U/mL
  sample_val_2 <- val_2_2[sample(1:nrow(val_2_2), nrow(val_2_2), replace = TRUE), ]
  sample_val_2_x_vars <- model.matrix(IgGpostDoV~. , sample_val_2)[,-1]
  sample_val_2_y_var <- sample_val_2$IgGpostDoV
  
  sample_val_2_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_2_x_vars)
  sample_val_2_min_roc <- roc(sample_val_2_y_var ~ as.numeric(sample_val_2_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_2_metric_min <- coords(sample_val_2_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_2_min_auc <- as.numeric(sample_val_2_min_roc$auc)
  sample_val_2_metric_min[6]<-sample_val_2_min_auc
  df_sample_val_2_min <- rbind(df_sample_val_2_min, sample_val_2_metric_min)
  
  sample_val_2_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_2_x_vars)
  sample_val_2_1se_roc <- roc(sample_val_2_y_var ~ as.numeric(sample_val_2_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_2_metric_1se <- coords(sample_val_2_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_2_1se_auc <- as.numeric(sample_val_2_1se_roc$auc)
  sample_val_2_metric_1se[6]<-sample_val_2_1se_auc
  df_sample_val_2_1se <- rbind(df_sample_val_2_1se, sample_val_2_metric_1se)
  
  # Validation Set 2 Mean/Median imputed Cutoff 15 U/mL
  sample_val_2_15 <- val_2_15_2[sample(1:nrow(val_2_15_2), nrow(val_2_15_2), replace = TRUE), ]
  sample_val_2_15_x_vars <- model.matrix(IgGpostDoV~. , sample_val_2_15)[,-1]
  sample_val_2_15_y_var <- sample_val_2_15$IgGpostDoV
  
  sample_val_2_15_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_2_15_x_vars)
  sample_val_2_15_min_roc <- roc(sample_val_2_15_y_var ~ as.numeric(sample_val_2_15_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_2_15_metric_min <- coords(sample_val_2_15_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_2_15_min_auc <- as.numeric(sample_val_2_15_min_roc$auc)
  sample_val_2_15_metric_min[6]<-sample_val_2_15_min_auc
  df_sample_val_2_15_min <- rbind(df_sample_val_2_15_min, sample_val_2_15_metric_min)
  
  sample_val_2_15_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_2_15_x_vars)
  sample_val_2_15_1se_roc <- roc(sample_val_2_15_y_var ~ as.numeric(sample_val_2_15_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_2_15_metric_1se <- coords(sample_val_2_15_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_2_15_1se_auc <- as.numeric(sample_val_2_15_1se_roc$auc)
  sample_val_2_15_metric_1se[6]<-sample_val_2_15_1se_auc
  df_sample_val_2_15_1se <- rbind(df_sample_val_2_15_1se, sample_val_2_15_metric_1se)
  
  # Validation Set 3 Mean/Median imputed
  sample_val_3 <- val_3_2[sample(1:nrow(val_3_2), nrow(val_3_2), replace = TRUE), ]
  sample_val_3_x_vars <- model.matrix(IgGpostDoV~. , sample_val_3)[,-1]
  sample_val_3_y_var <- sample_val_3$IgGpostDoV
  
  sample_val_3_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_3_x_vars)
  sample_val_3_min_roc <- roc(sample_val_3_y_var ~ as.numeric(sample_val_3_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_3_metric_min <- coords(sample_val_3_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_3_min_auc <- as.numeric(sample_val_3_min_roc$auc)
  sample_val_3_metric_min[6]<-sample_val_3_min_auc
  df_sample_val_3_min <- rbind(df_sample_val_3_min, sample_val_3_metric_min)
  
  sample_val_3_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_3_x_vars)
  sample_val_3_1se_roc <- roc(sample_val_3_y_var ~ as.numeric(sample_val_3_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_3_metric_1se <- coords(sample_val_3_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_3_1se_auc <- as.numeric(sample_val_3_1se_roc$auc)
  sample_val_3_metric_1se[6]<-sample_val_3_1se_auc
  df_sample_val_3_1se <- rbind(df_sample_val_3_1se, sample_val_3_metric_1se)
  
  # Validation Set 3 Complete Case analysis LASSO Min
  sample_val_3_cc_min <- val_3_cc_min[sample(1:nrow(val_3_cc_min), nrow(val_3_cc_min), replace = TRUE), ]
  sample_val_3_cc_min_x_vars <- model.matrix(IgGpostDoV~. , sample_val_3_cc_min)[,-1]
  sample_val_3_cc_min_y_var <- sample_val_3_cc_min$IgGpostDoV
  
  sample_val_3_cc_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_3_cc_min_x_vars)
  sample_val_3_cc_min_roc <- roc(sample_val_3_cc_min_y_var ~ as.numeric(sample_val_3_cc_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_3_cc_min_metric <- coords(sample_val_3_cc_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_3_cc_min_auc <- as.numeric(sample_val_3_cc_min_roc$auc)
  sample_val_3_cc_min_metric[6]<-sample_val_3_cc_min_auc
  df_sample_val_3_cc_min <- rbind(df_sample_val_3_cc_min, sample_val_3_cc_min_metric)
  
  # Validation Set 3 Complete Case analysis LASSO 1SE
  sample_val_3_cc_1se <- val_3_cc_1se[sample(1:nrow(val_3_cc_1se), nrow(val_3_cc_1se), replace = TRUE), ]
  sample_val_3_cc_1se_x_vars <- model.matrix(IgGpostDoV~. , sample_val_3_cc_1se)[,-1]
  sample_val_3_cc_1se_y_var <- sample_val_3_cc_1se$IgGpostDoV
  
  sample_val_3_cc_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_3_cc_1se_x_vars)
  sample_val_3_cc_1se_roc <- roc(sample_val_3_cc_1se_y_var ~ as.numeric(sample_val_3_cc_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_3_cc_1se_metric <- coords(sample_val_3_cc_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_3_cc_1se_auc <- as.numeric(sample_val_3_cc_1se_roc$auc)
  sample_val_3_cc_1se_metric[6]<-sample_val_3_cc_1se_auc
  df_sample_val_3_cc_1se <- rbind(df_sample_val_3_cc_1se, sample_val_3_cc_1se_metric)
  
  # Validation Set 4 Mean/Median imputed Cutoff 0.8 U/mL
  sample_val_4 <- val_4_2[sample(1:nrow(val_4_2), nrow(val_4_2), replace = TRUE), ]
  sample_val_4_x_vars <- model.matrix(IgGpostDoV~. , sample_val_4)[,-1]
  sample_val_4_y_var <- sample_val_4$IgGpostDoV
  
  sample_val_4_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_4_x_vars)
  sample_val_4_min_roc <- roc(sample_val_4_y_var ~ as.numeric(sample_val_4_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_4_metric_min <- coords(sample_val_4_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_4_min_auc <- as.numeric(sample_val_4_min_roc$auc)
  sample_val_4_metric_min[6]<-sample_val_4_min_auc
  df_sample_val_4_min <- rbind(df_sample_val_4_min, sample_val_4_metric_min)
  
  sample_val_4_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_4_x_vars)
  sample_val_4_1se_roc <- roc(sample_val_4_y_var ~ as.numeric(sample_val_4_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_4_metric_1se <- coords(sample_val_4_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_4_1se_auc <- as.numeric(sample_val_4_1se_roc$auc)
  sample_val_4_metric_1se[6]<-sample_val_4_1se_auc
  df_sample_val_4_1se <- rbind(df_sample_val_4_1se, sample_val_4_metric_1se)
  
  # Validation Set 4 Complete Case analysis LASSO Min Cutoff 0.8 U/mL
  sample_val_4_cc_min <- val_4_cc_min[sample(1:nrow(val_4_cc_min), nrow(val_4_cc_min), replace = TRUE), ]
  sample_val_4_cc_min_x_vars <- model.matrix(IgGpostDoV~. , sample_val_4_cc_min)[,-1]
  sample_val_4_cc_min_y_var <- sample_val_4_cc_min$IgGpostDoV
  
  sample_val_4_cc_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_4_cc_min_x_vars)
  sample_val_4_cc_min_roc <- roc(sample_val_4_cc_min_y_var ~ as.numeric(sample_val_4_cc_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_4_cc_min_metric <- coords(sample_val_4_cc_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_4_cc_min_auc <- as.numeric(sample_val_4_cc_min_roc$auc)
  sample_val_4_cc_min_metric[6]<-sample_val_4_cc_min_auc
  df_sample_val_4_cc_min <- rbind(df_sample_val_4_cc_min, sample_val_4_cc_min_metric)
  
  # Validation Set 4 Complete Case analysis LASSO 1SE Cutoff 0.8 U/mL
  sample_val_4_cc_1se <- val_4_cc_1se[sample(1:nrow(val_4_cc_1se), nrow(val_4_cc_1se), replace = TRUE), ]
  sample_val_4_cc_1se_x_vars <- model.matrix(IgGpostDoV~. , sample_val_4_cc_1se)[,-1]
  sample_val_4_cc_1se_y_var <- sample_val_4_cc_1se$IgGpostDoV
  
  sample_val_4_cc_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_4_cc_1se_x_vars)
  sample_val_4_cc_1se_roc <- roc(sample_val_4_cc_1se_y_var ~ as.numeric(sample_val_4_cc_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_4_cc_1se_metric <- coords(sample_val_4_cc_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_4_cc_1se_auc <- as.numeric(sample_val_4_cc_1se_roc$auc)
  sample_val_4_cc_1se_metric[6]<-sample_val_4_cc_1se_auc
  df_sample_val_4_cc_1se <- rbind(df_sample_val_4_cc_1se, sample_val_4_cc_1se_metric)
  
  # Validation Set 4 Mean/Median imputed Cutoff 15 U/mL
  sample_val_4_15 <- val_4_15_2[sample(1:nrow(val_4_15_2), nrow(val_4_15_2), replace = TRUE), ]
  sample_val_4_15_x_vars <- model.matrix(IgGpostDoV~. , sample_val_4_15)[,-1]
  sample_val_4_15_y_var <- sample_val_4_15$IgGpostDoV
  
  sample_val_4_15_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_4_15_x_vars)
  sample_val_4_15_min_roc <- roc(sample_val_4_15_y_var ~ as.numeric(sample_val_4_15_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_4_15_metric_min <- coords(sample_val_4_15_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_4_15_min_auc <- as.numeric(sample_val_4_15_min_roc$auc)
  sample_val_4_15_metric_min[6]<-sample_val_4_15_min_auc
  df_sample_val_4_15_min <- rbind(df_sample_val_4_15_min, sample_val_4_15_metric_min)
  
  sample_val_4_15_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_4_15_x_vars)
  sample_val_4_15_1se_roc <- roc(sample_val_4_15_y_var ~ as.numeric(sample_val_4_15_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_4_15_metric_1se <- coords(sample_val_4_15_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_4_15_1se_auc <- as.numeric(sample_val_4_15_1se_roc$auc)
  sample_val_4_15_metric_1se[6]<-sample_val_4_15_1se_auc
  df_sample_val_4_15_1se <- rbind(df_sample_val_4_15_1se, sample_val_4_15_metric_1se)
  
  # Validation Set 4 Complete Case analysis LASSO Min Cutoff 15 U/mL
  sample_val_4_15_cc_min <- val_4_15_cc_min[sample(1:nrow(val_4_15_cc_min), nrow(val_4_15_cc_min), replace = TRUE), ]
  sample_val_4_15_cc_min_x_vars <- model.matrix(IgGpostDoV~. , sample_val_4_15_cc_min)[,-1]
  sample_val_4_15_cc_min_y_var <- sample_val_4_15_cc_min$IgGpostDoV
  
  sample_val_4_15_cc_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_4_15_cc_min_x_vars)
  sample_val_4_15_cc_min_roc <- roc(sample_val_4_15_cc_min_y_var ~ as.numeric(sample_val_4_15_cc_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_4_15_cc_min_metric <- coords(sample_val_4_15_cc_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_4_15_cc_min_auc <- as.numeric(sample_val_4_15_cc_min_roc$auc)
  sample_val_4_15_cc_min_metric[6]<-sample_val_4_15_cc_min_auc
  df_sample_val_4_15_cc_min <- rbind(df_sample_val_4_15_cc_min, sample_val_4_15_cc_min_metric)
  
  # Validation Set 4 Complete Case analysis LASSO 1SE Cutoff 15 U/mL
  sample_val_4_15_cc_1se <- val_4_15_cc_1se[sample(1:nrow(val_4_15_cc_1se), nrow(val_4_15_cc_1se), replace = TRUE), ]
  sample_val_4_15_cc_1se_x_vars <- model.matrix(IgGpostDoV~. , sample_val_4_15_cc_1se)[,-1]
  sample_val_4_15_cc_1se_y_var <- sample_val_4_15_cc_1se$IgGpostDoV
  
  sample_val_4_15_cc_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_4_15_cc_1se_x_vars)
  sample_val_4_15_cc_1se_roc <- roc(sample_val_4_15_cc_1se_y_var ~ as.numeric(sample_val_4_15_cc_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_4_15_cc_1se_metric <- coords(sample_val_4_15_cc_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_4_15_cc_1se_auc <- as.numeric(sample_val_4_15_cc_1se_roc$auc)
  sample_val_4_15_cc_1se_metric[6]<-sample_val_4_15_cc_1se_auc
  df_sample_val_4_15_cc_1se <- rbind(df_sample_val_4_15_cc_1se, sample_val_4_15_cc_1se_metric)
  
  sample_val_all <- val_all[sample(1:nrow(val_all), nrow(val_all), replace = TRUE), ]
  sample_val_all_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all)[,-1]
  sample_val_all_y_var <- sample_val_all$IgGpostDoV
  
  sample_val_all_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_x_vars)
  sample_val_all_min_roc <- roc(sample_val_all_y_var ~ as.numeric(sample_val_all_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_metric_min <- coords(sample_val_all_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_min_auc <- as.numeric(sample_val_all_min_roc$auc)
  sample_val_all_metric_min[6]<-sample_val_all_min_auc
  df_sample_val_all_min <- rbind(df_sample_val_all_min, sample_val_all_metric_min)
  
  sample_val_all_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_x_vars)
  sample_val_all_1se_roc <- roc(sample_val_all_y_var ~ as.numeric(sample_val_all_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_metric_1se <- coords(sample_val_all_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_1se_auc <- as.numeric(sample_val_all_1se_roc$auc)
  sample_val_all_metric_1se[6]<-sample_val_all_1se_auc
  df_sample_val_all_1se <- rbind(df_sample_val_all_1se, sample_val_all_metric_1se)
  
  # Random Forest
  sample_val_all_rf <- sample_val_all
  sample_val_all_rf$IgGpostDoV <- as.factor(sample_val_all_rf$IgGpostDoV)
  
  sample_val_all_rf_score <- predict(rf_lympho, newdata=sample_val_all_rf, type="class")
  sample_val_all_rf_roc <- roc(sample_val_all_rf$IgGpostDoV ~ sample_val_all_rf_score, plot = FALSE, print.auc = FALSE)
  sample_val_all_metric_rf <- coords(sample_val_all_rf_roc, 0.5, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_rf_auc <- as.numeric(sample_val_all_rf_roc$auc)
  sample_val_all_metric_rf[6]<-sample_val_all_rf_auc
  df_sample_val_all_rf <- rbind(df_sample_val_all_rf, sample_val_all_metric_rf)  
  
  # GBM
  sample_val_all_gbm_score<-predict(gbm_lympho, newdata = sample_val_all_rf)
  sample_val_all_gbm_roc <- roc(sample_val_all_rf$IgGpostDoV ~ sample_val_all_gbm_score, plot = FALSE, print.auc = FALSE)
  sample_val_all_metric_gbm <- coords(sample_val_all_gbm_roc, cutoff_gbm, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_gbm_auc <- as.numeric(sample_val_all_gbm_roc$auc)
  sample_val_all_metric_gbm[6]<-sample_val_all_gbm_auc
  df_sample_val_all_gbm <- rbind(df_sample_val_all_gbm, sample_val_all_metric_gbm)  
  
  # LR (full)
  sample_val_all_lr_score <- predict(lr_lympho, newdata = sample_val_all_rf)
  sample_val_all_lr_roc <- roc(sample_val_all_rf$IgGpostDoV ~ as.numeric(sample_val_all_lr_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_metric_lr <- coords(sample_val_all_lr_roc, cutoff_lr, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_lr_auc <- as.numeric(sample_val_all_lr_roc$auc)
  sample_val_all_metric_lr[6]<-sample_val_all_lr_auc
  df_sample_val_all_lr <- rbind(df_sample_val_all_lr, sample_val_all_metric_lr)
  
  ## Imp 1 Cutoff 0.8
  sample_val_all_mi_1 <- val_all_mi_1[sample(1:nrow(val_all_mi_1), nrow(val_all_mi_1), replace = TRUE), ]
  sample_val_all_mi_1_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_mi_1)[,-1]
  sample_val_all_mi_1_y_var <- as.integer(sample_val_all_mi_1$IgGpostDoV)
  
  sample_val_all_mi_1_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_mi_1_x_vars)
  sample_val_all_mi_1_min_roc <- roc(sample_val_all_mi_1_y_var ~ as.numeric(sample_val_all_mi_1_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_1_metric_min <- coords(sample_val_all_mi_1_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_1_min_auc <- as.numeric(sample_val_all_mi_1_min_roc$auc)
  sample_val_all_mi_1_metric_min[6]<-sample_val_all_mi_1_min_auc
  
  
  sample_val_all_mi_1_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_mi_1_x_vars)
  sample_val_all_mi_1_1se_roc <- roc(sample_val_all_mi_1_y_var ~ as.numeric(sample_val_all_mi_1_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_1_metric_1se <- coords(sample_val_all_mi_1_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_1_1se_auc <- as.numeric(sample_val_all_mi_1_1se_roc$auc)
  sample_val_all_mi_1_metric_1se[6]<-sample_val_all_mi_1_1se_auc

  ## Imp 2 Cutoff 0.8
  sample_val_all_mi_2 <- val_all_mi_2[sample(1:nrow(val_all_mi_2), nrow(val_all_mi_2), replace = TRUE), ]
  sample_val_all_mi_2_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_mi_2)[,-1]
  sample_val_all_mi_2_y_var <- as.integer(sample_val_all_mi_2$IgGpostDoV)
  
  sample_val_all_mi_2_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_mi_2_x_vars)
  sample_val_all_mi_2_min_roc <- roc(sample_val_all_mi_2_y_var ~ as.numeric(sample_val_all_mi_2_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_2_metric_min <- coords(sample_val_all_mi_2_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_2_min_auc <- as.numeric(sample_val_all_mi_2_min_roc$auc)
  sample_val_all_mi_2_metric_min[6]<-sample_val_all_mi_2_min_auc
  
  
  sample_val_all_mi_2_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_mi_2_x_vars)
  sample_val_all_mi_2_1se_roc <- roc(sample_val_all_mi_2_y_var ~ as.numeric(sample_val_all_mi_2_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_2_metric_1se <- coords(sample_val_all_mi_2_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_2_1se_auc <- as.numeric(sample_val_all_mi_2_1se_roc$auc)
  sample_val_all_mi_2_metric_1se[6]<-sample_val_all_mi_2_1se_auc
  
  ## Imp 3 Cutoff 0.8
  sample_val_all_mi_3 <- val_all_mi_3[sample(1:nrow(val_all_mi_3), nrow(val_all_mi_3), replace = TRUE), ]
  sample_val_all_mi_3_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_mi_3)[,-1]
  sample_val_all_mi_3_y_var <- as.integer(sample_val_all_mi_3$IgGpostDoV)
  
  sample_val_all_mi_3_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_mi_3_x_vars)
  sample_val_all_mi_3_min_roc <- roc(sample_val_all_mi_3_y_var ~ as.numeric(sample_val_all_mi_3_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_3_metric_min <- coords(sample_val_all_mi_3_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_3_min_auc <- as.numeric(sample_val_all_mi_3_min_roc$auc)
  sample_val_all_mi_3_metric_min[6]<-sample_val_all_mi_3_min_auc
  
  
  sample_val_all_mi_3_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_mi_3_x_vars)
  sample_val_all_mi_3_1se_roc <- roc(sample_val_all_mi_3_y_var ~ as.numeric(sample_val_all_mi_3_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_3_metric_1se <- coords(sample_val_all_mi_3_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_3_1se_auc <- as.numeric(sample_val_all_mi_3_1se_roc$auc)
  sample_val_all_mi_3_metric_1se[6]<-sample_val_all_mi_3_1se_auc
  
  ## Imp 4 Cutoff 0.8
  sample_val_all_mi_4 <- val_all_mi_4[sample(1:nrow(val_all_mi_4), nrow(val_all_mi_4), replace = TRUE), ]
  sample_val_all_mi_4_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_mi_4)[,-1]
  sample_val_all_mi_4_y_var <- as.integer(sample_val_all_mi_4$IgGpostDoV)
  
  sample_val_all_mi_4_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_mi_4_x_vars)
  sample_val_all_mi_4_min_roc <- roc(sample_val_all_mi_4_y_var ~ as.numeric(sample_val_all_mi_4_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_4_metric_min <- coords(sample_val_all_mi_4_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_4_min_auc <- as.numeric(sample_val_all_mi_4_min_roc$auc)
  sample_val_all_mi_4_metric_min[6]<-sample_val_all_mi_4_min_auc
  
  
  sample_val_all_mi_4_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_mi_4_x_vars)
  sample_val_all_mi_4_1se_roc <- roc(sample_val_all_mi_4_y_var ~ as.numeric(sample_val_all_mi_4_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_4_metric_1se <- coords(sample_val_all_mi_4_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_4_1se_auc <- as.numeric(sample_val_all_mi_4_1se_roc$auc)
  sample_val_all_mi_4_metric_1se[6]<-sample_val_all_mi_4_1se_auc
  
  ## Imp 5 Cutoff 0.8
  sample_val_all_mi_5 <- val_all_mi_5[sample(1:nrow(val_all_mi_5), nrow(val_all_mi_5), replace = TRUE), ]
  sample_val_all_mi_5_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_mi_5)[,-1]
  sample_val_all_mi_5_y_var <- as.integer(sample_val_all_mi_5$IgGpostDoV)
  
  sample_val_all_mi_5_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_mi_5_x_vars)
  sample_val_all_mi_5_min_roc <- roc(sample_val_all_mi_5_y_var ~ as.numeric(sample_val_all_mi_5_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_5_metric_min <- coords(sample_val_all_mi_5_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_5_min_auc <- as.numeric(sample_val_all_mi_5_min_roc$auc)
  sample_val_all_mi_5_metric_min[6]<-sample_val_all_mi_5_min_auc
  
  
  sample_val_all_mi_5_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_mi_5_x_vars)
  sample_val_all_mi_5_1se_roc <- roc(sample_val_all_mi_5_y_var ~ as.numeric(sample_val_all_mi_5_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_mi_5_metric_1se <- coords(sample_val_all_mi_5_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_mi_5_1se_auc <- as.numeric(sample_val_all_mi_5_1se_roc$auc)
  sample_val_all_mi_5_metric_1se[6]<-sample_val_all_mi_5_1se_auc
  
  sample_val_all_mi_metric_1se <- rbind(sample_val_all_mi_1_metric_1se, sample_val_all_mi_2_metric_1se, sample_val_all_mi_3_metric_1se, sample_val_all_mi_4_metric_1se, sample_val_all_mi_5_metric_1se)
  sample_val_all_mi_metric_min <- rbind(sample_val_all_mi_1_metric_min, sample_val_all_mi_2_metric_min, sample_val_all_mi_3_metric_min, sample_val_all_mi_4_metric_min, sample_val_all_mi_5_metric_min)

  df_sample_val_all_mi_min <- rbind(df_sample_val_all_mi_min,   colMeans(sample_val_all_mi_metric_min))
  df_sample_val_all_mi_1se <- rbind(df_sample_val_all_mi_1se, colMeans(sample_val_all_mi_metric_1se))
  
  
  # All 15 U/mL
  sample_val_all_15 <- val_all_15[sample(1:nrow(val_all_15), nrow(val_all_15), replace = TRUE), ]
  sample_val_all_15_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_15)[,-1]
  sample_val_all_15_y_var <- sample_val_all_15$IgGpostDoV
  
  sample_val_all_15_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_15_x_vars)
  sample_val_all_15_min_roc <- roc(sample_val_all_15_y_var ~ as.numeric(sample_val_all_15_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_metric_min <- coords(sample_val_all_15_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_min_auc <- as.numeric(sample_val_all_15_min_roc$auc)
  sample_val_all_15_metric_min[6]<-sample_val_all_15_min_auc
  df_sample_val_all_15_min <- rbind(df_sample_val_all_15_min, sample_val_all_15_metric_min)
  
  sample_val_all_15_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_15_x_vars)
  sample_val_all_15_1se_roc <- roc(sample_val_all_15_y_var ~ as.numeric(sample_val_all_15_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_metric_1se <- coords(sample_val_all_15_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_1se_auc <- as.numeric(sample_val_all_15_1se_roc$auc)
  sample_val_all_15_metric_1se[6]<-sample_val_all_15_1se_auc
  df_sample_val_all_15_1se <- rbind(df_sample_val_all_15_1se, sample_val_all_15_metric_1se)
  
  # Random Forest
  sample_val_all_15_rf <- sample_val_all_15
  sample_val_all_15_rf$IgGpostDoV <- as.factor(sample_val_all_15_rf$IgGpostDoV)
  
  sample_val_all_15_rf_score <- predict(rf_lympho, newdata=sample_val_all_15_rf, type="class")
  sample_val_all_15_rf_roc <- roc(sample_val_all_15_rf$IgGpostDoV ~ sample_val_all_15_rf_score, plot = FALSE, print.auc = FALSE)
  sample_val_all_15_metric_rf <- coords(sample_val_all_15_rf_roc, 0.5, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_rf_auc <- as.numeric(sample_val_all_15_rf_roc$auc)
  sample_val_all_15_metric_rf[6]<-sample_val_all_15_rf_auc
  df_sample_val_all_15_rf <- rbind(df_sample_val_all_15_rf, sample_val_all_15_metric_rf)  
  
  # GBM
  sample_val_all_15_gbm_score<-predict(gbm_lympho, newdata = sample_val_all_15_rf)
  sample_val_all_15_gbm_roc <- roc(sample_val_all_15_rf$IgGpostDoV ~ sample_val_all_15_gbm_score, plot = FALSE, print.auc = FALSE)
  sample_val_all_15_metric_gbm <- coords(sample_val_all_15_gbm_roc, cutoff_gbm, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_gbm_auc <- as.numeric(sample_val_all_15_gbm_roc$auc)
  sample_val_all_15_metric_gbm[6]<-sample_val_all_15_gbm_auc
  df_sample_val_all_15_gbm <- rbind(df_sample_val_all_15_gbm, sample_val_all_15_metric_gbm)  
  
  # LR (full)
  sample_val_all_15_lr_score <- predict(lr_lympho, newdata = sample_val_all_15_rf)
  sample_val_all_15_lr_roc <- roc(sample_val_all_15_rf$IgGpostDoV ~ as.numeric(sample_val_all_15_lr_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_metric_lr <- coords(sample_val_all_15_lr_roc, cutoff_lr, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_lr_auc <- as.numeric(sample_val_all_15_lr_roc$auc)
  sample_val_all_15_metric_lr[6]<-sample_val_all_15_lr_auc
  df_sample_val_all_15_lr <- rbind(df_sample_val_all_15_lr, sample_val_all_15_metric_lr)
  
  ## Imp 1 Cutoff 15
  sample_val_all_15_mi_1 <- val_all_15_mi_1[sample(1:nrow(val_all_15_mi_1), nrow(val_all_15_mi_1), replace = TRUE), ]
  sample_val_all_15_mi_1_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_15_mi_1)[,-1]
  sample_val_all_15_mi_1_y_var <- as.integer(sample_val_all_15_mi_1$IgGpostDoV)
  
  sample_val_all_15_mi_1_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_15_mi_1_x_vars)
  sample_val_all_15_mi_1_min_roc <- roc(sample_val_all_15_mi_1_y_var ~ as.numeric(sample_val_all_15_mi_1_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_1_metric_min <- coords(sample_val_all_15_mi_1_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_1_min_auc <- as.numeric(sample_val_all_15_mi_1_min_roc$auc)
  sample_val_all_15_mi_1_metric_min[6]<-sample_val_all_15_mi_1_min_auc
  
  
  sample_val_all_15_mi_1_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_15_mi_1_x_vars)
  sample_val_all_15_mi_1_1se_roc <- roc(sample_val_all_15_mi_1_y_var ~ as.numeric(sample_val_all_15_mi_1_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_1_metric_1se <- coords(sample_val_all_15_mi_1_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_1_1se_auc <- as.numeric(sample_val_all_15_mi_1_1se_roc$auc)
  sample_val_all_15_mi_1_metric_1se[6]<-sample_val_all_15_mi_1_1se_auc
  
  ## Imp 2 Cutoff 15
  sample_val_all_15_mi_2 <- val_all_15_mi_2[sample(1:nrow(val_all_15_mi_2), nrow(val_all_15_mi_2), replace = TRUE), ]
  sample_val_all_15_mi_2_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_15_mi_2)[,-1]
  sample_val_all_15_mi_2_y_var <- as.integer(sample_val_all_15_mi_2$IgGpostDoV)
  
  sample_val_all_15_mi_2_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_15_mi_2_x_vars)
  sample_val_all_15_mi_2_min_roc <- roc(sample_val_all_15_mi_2_y_var ~ as.numeric(sample_val_all_15_mi_2_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_2_metric_min <- coords(sample_val_all_15_mi_2_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_2_min_auc <- as.numeric(sample_val_all_15_mi_2_min_roc$auc)
  sample_val_all_15_mi_2_metric_min[6]<-sample_val_all_15_mi_2_min_auc
  
  
  sample_val_all_15_mi_2_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_15_mi_2_x_vars)
  sample_val_all_15_mi_2_1se_roc <- roc(sample_val_all_15_mi_2_y_var ~ as.numeric(sample_val_all_15_mi_2_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_2_metric_1se <- coords(sample_val_all_15_mi_2_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_2_1se_auc <- as.numeric(sample_val_all_15_mi_2_1se_roc$auc)
  sample_val_all_15_mi_2_metric_1se[6]<-sample_val_all_15_mi_2_1se_auc
  
  ## Imp 3 Cutoff 15
  sample_val_all_15_mi_3 <- val_all_15_mi_3[sample(1:nrow(val_all_15_mi_3), nrow(val_all_15_mi_3), replace = TRUE), ]
  sample_val_all_15_mi_3_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_15_mi_3)[,-1]
  sample_val_all_15_mi_3_y_var <- as.integer(sample_val_all_15_mi_3$IgGpostDoV)
  
  sample_val_all_15_mi_3_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_15_mi_3_x_vars)
  sample_val_all_15_mi_3_min_roc <- roc(sample_val_all_15_mi_3_y_var ~ as.numeric(sample_val_all_15_mi_3_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_3_metric_min <- coords(sample_val_all_15_mi_3_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_3_min_auc <- as.numeric(sample_val_all_15_mi_3_min_roc$auc)
  sample_val_all_15_mi_3_metric_min[6]<-sample_val_all_15_mi_3_min_auc
  
  
  sample_val_all_15_mi_3_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_15_mi_3_x_vars)
  sample_val_all_15_mi_3_1se_roc <- roc(sample_val_all_15_mi_3_y_var ~ as.numeric(sample_val_all_15_mi_3_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_3_metric_1se <- coords(sample_val_all_15_mi_3_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_3_1se_auc <- as.numeric(sample_val_all_15_mi_3_1se_roc$auc)
  sample_val_all_15_mi_3_metric_1se[6]<-sample_val_all_15_mi_3_1se_auc
  
  ## Imp 4 Cutoff 15
  sample_val_all_15_mi_4 <- val_all_15_mi_4[sample(1:nrow(val_all_15_mi_4), nrow(val_all_15_mi_4), replace = TRUE), ]
  sample_val_all_15_mi_4_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_15_mi_4)[,-1]
  sample_val_all_15_mi_4_y_var <- as.integer(sample_val_all_15_mi_4$IgGpostDoV)
  
  sample_val_all_15_mi_4_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_15_mi_4_x_vars)
  sample_val_all_15_mi_4_min_roc <- roc(sample_val_all_15_mi_4_y_var ~ as.numeric(sample_val_all_15_mi_4_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_4_metric_min <- coords(sample_val_all_15_mi_4_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_4_min_auc <- as.numeric(sample_val_all_15_mi_4_min_roc$auc)
  sample_val_all_15_mi_4_metric_min[6]<-sample_val_all_15_mi_4_min_auc
  
  
  sample_val_all_15_mi_4_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_15_mi_4_x_vars)
  sample_val_all_15_mi_4_1se_roc <- roc(sample_val_all_15_mi_4_y_var ~ as.numeric(sample_val_all_15_mi_4_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_4_metric_1se <- coords(sample_val_all_15_mi_4_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_4_1se_auc <- as.numeric(sample_val_all_15_mi_4_1se_roc$auc)
  sample_val_all_15_mi_4_metric_1se[6]<-sample_val_all_15_mi_4_1se_auc
  
  ## Imp 5 Cutoff 15
  sample_val_all_15_mi_5 <- val_all_15_mi_5[sample(1:nrow(val_all_15_mi_5), nrow(val_all_15_mi_5), replace = TRUE), ]
  sample_val_all_15_mi_5_x_vars <- model.matrix(IgGpostDoV~. , sample_val_all_15_mi_5)[,-1]
  sample_val_all_15_mi_5_y_var <- as.integer(sample_val_all_15_mi_5$IgGpostDoV)
  
  sample_val_all_15_mi_5_min_score <- predict(lasso_min_lympho, s = lassologcv_lympho$lambda.min, newx = sample_val_all_15_mi_5_x_vars)
  sample_val_all_15_mi_5_min_roc <- roc(sample_val_all_15_mi_5_y_var ~ as.numeric(sample_val_all_15_mi_5_min_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_5_metric_min <- coords(sample_val_all_15_mi_5_min_roc, cutoff_min, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_5_min_auc <- as.numeric(sample_val_all_15_mi_5_min_roc$auc)
  sample_val_all_15_mi_5_metric_min[6]<-sample_val_all_15_mi_5_min_auc
  
  
  sample_val_all_15_mi_5_1se_score <- predict(lasso_1se_lympho, s = lassologcv_lympho$lambda.1se, newx = sample_val_all_15_mi_5_x_vars)
  sample_val_all_15_mi_5_1se_roc <- roc(sample_val_all_15_mi_5_y_var ~ as.numeric(sample_val_all_15_mi_5_1se_score), plot = FALSE, print.auc = FALSE) # ROC-Analyse Vergleich mit echten Outcomes 
  sample_val_all_15_mi_5_metric_1se <- coords(sample_val_all_15_mi_5_1se_roc, cutoff_1se, input="threshold", ret=c("sensitivity", "specificity", "accuracy", "ppv", "npv"), transpose = T)
  sample_val_all_15_mi_5_1se_auc <- as.numeric(sample_val_all_15_mi_5_1se_roc$auc)
  sample_val_all_15_mi_5_metric_1se[6]<-sample_val_all_15_mi_5_1se_auc
  
  sample_val_all_15_mi_metric_1se <- rbind(sample_val_all_15_mi_1_metric_1se, sample_val_all_15_mi_2_metric_1se, sample_val_all_15_mi_3_metric_1se, sample_val_all_15_mi_4_metric_1se, sample_val_all_15_mi_5_metric_1se)
  sample_val_all_15_mi_metric_min <- rbind(sample_val_all_15_mi_1_metric_min, sample_val_all_15_mi_2_metric_min, sample_val_all_15_mi_3_metric_min, sample_val_all_15_mi_4_metric_min, sample_val_all_15_mi_5_metric_min)
  
  df_sample_val_all_15_mi_min <- rbind(df_sample_val_all_15_mi_min,   colMeans(sample_val_all_15_mi_metric_min))
  df_sample_val_all_15_mi_1se <- rbind(df_sample_val_all_15_mi_1se, colMeans(sample_val_all_15_mi_metric_1se))
}

### Rename Column Names correctly
df_sample_val_1_min <- setNames(df_sample_val_1_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_1_1se <- setNames(df_sample_val_1_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_1_cc_min <- setNames(df_sample_val_1_cc_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_1_cc_1se <- setNames(df_sample_val_1_cc_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_2_min <- setNames(df_sample_val_2_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_2_1se <- setNames(df_sample_val_2_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_2_15_min <- setNames(df_sample_val_2_15_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_2_15_1se <- setNames(df_sample_val_2_15_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_3_min <- setNames(df_sample_val_3_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_3_1se <- setNames(df_sample_val_3_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_3_cc_min <- setNames(df_sample_val_3_cc_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_3_cc_1se <- setNames(df_sample_val_3_cc_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_4_min <- setNames(df_sample_val_4_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_4_1se <- setNames(df_sample_val_4_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_4_cc_min <- setNames(df_sample_val_4_cc_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_4_cc_1se <- setNames(df_sample_val_4_cc_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_4_15_min <- setNames(df_sample_val_4_15_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_4_15_1se <- setNames(df_sample_val_4_15_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_4_15_cc_min <- setNames(df_sample_val_4_15_cc_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_4_15_cc_1se <- setNames(df_sample_val_4_15_cc_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_min <- setNames(df_sample_val_all_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_1se <- setNames(df_sample_val_all_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_lr <- setNames(df_sample_val_all_lr, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_rf <- setNames(df_sample_val_all_rf, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_gbm <- setNames(df_sample_val_all_gbm, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_mi_min <- setNames(df_sample_val_all_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_mi_1se <- setNames(df_sample_val_all_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_15_min <- setNames(df_sample_val_all_15_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_15_1se <- setNames(df_sample_val_all_15_1se, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_15_lr <- setNames(df_sample_val_all_15_lr, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_15_rf <- setNames(df_sample_val_all_15_rf, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_15_gbm <- setNames(df_sample_val_all_15_gbm, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_15_mi_min <- setNames(df_sample_val_all_15_min, c("sens","spec","acc","ppv","npv","auc"))
df_sample_val_all_15_mi_1se <- setNames(df_sample_val_all_15_1se, c("sens","spec","acc","ppv","npv","auc"))


### Add Column "class" with Model name
df_sample_val_1_min$class <- "Val-1 20-var"
df_sample_val_1_1se$class <- "Val-1 10-var"
df_sample_val_1_cc_min$class <- "Val-1 20-var (cc)"
df_sample_val_1_cc_1se$class <- "Val-1 12-var (cc)"
df_sample_val_2_min$class <- "Val-2 20-var"
df_sample_val_2_1se$class <- "Val-2 10-var"
df_sample_val_2_15_min$class <- "Val-2 20-var 15U/mL"
df_sample_val_2_15_1se$class <- "Val-2 10-var 15U/mL"
df_sample_val_3_min$class <- "Val-3 20-var"
df_sample_val_3_1se$class <- "Val-3 10-var"
df_sample_val_3_cc_min$class <- "Val-3 20-var (cc)"
df_sample_val_3_cc_1se$class <- "Val-3 10-var (cc)"
df_sample_val_4_min$class <- "Val-4 20-var"
df_sample_val_4_1se$class <- "Val-4 10-var"
df_sample_val_4_cc_min$class <- "Val-4 20-var (cc)"
df_sample_val_4_cc_1se$class <- "Val-4 10-var (cc)"
df_sample_val_4_15_min$class <- "Val-4 20-var 15U/mL"
df_sample_val_4_15_1se$class <- "Val-4 10-var 15U/mL"
df_sample_val_4_15_cc_min$class <- "Val-4 20-var 15U/mL (cc)"
df_sample_val_4_15_cc_1se$class <- "Val-4 10-var 15U/mL (cc)"
df_sample_val_all_min$class <- "all 20-var"
df_sample_val_all_1se$class <- "all 10-var"
df_sample_val_all_lr$class <- "all LR"
df_sample_val_all_rf$class <- "all RF"
df_sample_val_all_gbm$class <- "all GBRT"
df_sample_val_all_mi_min$class <- "all 20-var (MI)"
df_sample_val_all_mi_1se$class <- "all 10-var (MI)"
df_sample_val_all_15_min$class <- "all 20-var 15U/mL"
df_sample_val_all_15_1se$class <- "all 10-var 15U/mL"
df_sample_val_all_15_lr$class <- "all LR 15U/mL"
df_sample_val_all_15_rf$class <- "all RF 15U/mL"
df_sample_val_all_15_gbm$class <- "all GBRT 15U/mL"
df_sample_val_all_15_mi_min$class <- "all 20-var 15U/mL (MI)"
df_sample_val_all_15_mi_1se$class <- "all 10-var 15U/mL (MI)"

#Save Model Performance Data
write.csv(df_sample_val_1_min,"E:\\Impfungen_modell\\Final\\df_sample_val_1_min.csv", row.names = FALSE)
write.csv(df_sample_val_1_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_1_1se.csv", row.names = FALSE)
write.csv(df_sample_val_1_cc_min,"E:\\Impfungen_modell\\Final\\df_sample_val_1_cc_min.csv", row.names = FALSE)
write.csv(df_sample_val_1_cc_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_1_cc_1se.csv", row.names = FALSE)
write.csv(df_sample_val_2_min,"E:\\Impfungen_modell\\Final\\df_sample_val_2_min.csv", row.names = FALSE)
write.csv(df_sample_val_2_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_2_1se.csv", row.names = FALSE)
write.csv(df_sample_val_2_15_min,"E:\\Impfungen_modell\\Final\\df_sample_val_2_15_min.csv", row.names = FALSE)
write.csv(df_sample_val_2_15_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_2_15_1se.csv", row.names = FALSE)
write.csv(df_sample_val_3_min,"E:\\Impfungen_modell\\Final\\df_sample_val_3_min.csv", row.names = FALSE)
write.csv(df_sample_val_3_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_3_1se.csv", row.names = FALSE)
write.csv(df_sample_val_3_cc_min,"E:\\Impfungen_modell\\Final\\df_sample_val_3_cc_min.csv", row.names = FALSE)
write.csv(df_sample_val_3_cc_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_3_cc_1se.csv", row.names = FALSE)
write.csv(df_sample_val_4_min,"E:\\Impfungen_modell\\Final\\df_sample_val_4_min.csv", row.names = FALSE)
write.csv(df_sample_val_4_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_4_1se.csv", row.names = FALSE)
write.csv(df_sample_val_4_cc_min,"E:\\Impfungen_modell\\Final\\df_sample_val_4_cc_min.csv", row.names = FALSE)
write.csv(df_sample_val_4_cc_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_4_cc_1se.csv", row.names = FALSE)
write.csv(df_sample_val_4_15_min,"E:\\Impfungen_modell\\Final\\df_sample_val_4_15_min.csv", row.names = FALSE)
write.csv(df_sample_val_4_15_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_4_15_1se.csv", row.names = FALSE)
write.csv(df_sample_val_4_15_cc_min,"E:\\Impfungen_modell\\Final\\df_sample_val_4_15_cc_min.csv", row.names = FALSE)
write.csv(df_sample_val_4_15_cc_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_4_15_cc_1se.csv", row.names = FALSE)
write.csv(df_sample_val_all_min,"E:\\Impfungen_modell\\Final\\df_sample_val_all_min.csv", row.names = FALSE)
write.csv(df_sample_val_all_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_all_1se.csv", row.names = FALSE)
write.csv(df_sample_val_all_lr,"E:\\Impfungen_modell\\Final\\df_sample_val_all_lr.csv", row.names = FALSE)
write.csv(df_sample_val_all_rf,"E:\\Impfungen_modell\\Final\\df_sample_val_all_rf.csv", row.names = FALSE)
write.csv(df_sample_val_all_gbm,"E:\\Impfungen_modell\\Final\\df_sample_val_all_gbm.csv", row.names = FALSE)
write.csv(df_sample_val_all_mi_min,"E:\\Impfungen_modell\\Final\\df_sample_val_all_mi_min.csv", row.names = FALSE)
write.csv(df_sample_val_all_mi_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_all_mi_1se.csv", row.names = FALSE)
write.csv(df_sample_val_all_15_min,"E:\\Impfungen_modell\\Final\\df_sample_val_all_15_min.csv", row.names = FALSE)
write.csv(df_sample_val_all_15_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_all_15_1se.csv", row.names = FALSE)
write.csv(df_sample_val_all_15_lr,"E:\\Impfungen_modell\\Final\\df_sample_val_all_15_lr.csv", row.names = FALSE)
write.csv(df_sample_val_all_15_rf,"E:\\Impfungen_modell\\Final\\df_sample_val_all_15_rf.csv", row.names = FALSE)
write.csv(df_sample_val_all_15_gbm,"E:\\Impfungen_modell\\Final\\df_sample_val_all_15_gbm.csv", row.names = FALSE)
write.csv(df_sample_val_all_15_mi_min,"E:\\Impfungen_modell\\Final\\df_sample_val_all_15_mi_min.csv", row.names = FALSE)
write.csv(df_sample_val_all_15_mi_1se,"E:\\Impfungen_modell\\Final\\df_sample_val_all_15_mi_1se.csv", row.names = FALSE)


### Prepare Plotting for 10-variable model
df_val_10 <- rbind(df_sample_val_1_1se, df_sample_val_1_cc_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_2_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_2_15_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_3_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_3_cc_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_4_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_4_cc_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_4_15_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_4_15_cc_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_all_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_all_mi_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_all_15_1se)
df_val_10 <- rbind(df_val_10, df_sample_val_all_15_mi_1se)

plot_order_10 <- c("Val-1 10-var", "Val-1 10-var (cc)", "Val-2 10-var", "Val-2 10-var 15U/mL", "Val-3 10-var",  "Val-3 10-var (cc)", "Val-4 10-var", "Val-4 10-var (cc)", "Val-4 10-var 15U/mL", "Val-4 10-var 15U/mL (cc)", "all 10-var", "all 10-var (MI)",  "all 10-var 15U/mL", "all 10-var 15U/mL (MI)")

df_val_10$class <- factor(df_val_10$class, levels = plot_order_10)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

ggplot(data = df_val_10, mapping = aes(x=class, y=auc)) + 
  geom_point(aes(),alpha=0.8) +
  geom_boxplot(fill="grey",color="black",alpha=0.2) + 
  scale_x_discrete(breaks=unique(df_val_10$class), 
                   labels=addline_format(plot_order_10))+
  theme_minimal()


### Dataframes 20-variable model
df_val_20 <- rbind(df_sample_val_1_min, df_sample_val_1_cc_min)
df_val_20 <- rbind(df_val_20, df_sample_val_2_min)
df_val_20 <- rbind(df_val_20, df_sample_val_2_15_min)
df_val_20 <- rbind(df_val_20, df_sample_val_3_min)
df_val_20 <- rbind(df_val_20, df_sample_val_3_cc_min)
df_val_20 <- rbind(df_val_20, df_sample_val_4_min)
df_val_20 <- rbind(df_val_20, df_sample_val_4_cc_min)
df_val_20 <- rbind(df_val_20, df_sample_val_4_15_min)
df_val_20 <- rbind(df_val_20, df_sample_val_4_15_cc_min)
df_val_20 <- rbind(df_val_20, df_sample_val_all_min)
df_val_20 <- rbind(df_val_20, df_sample_val_all_mi_min)
df_val_20 <- rbind(df_val_20, df_sample_val_all_15_min)
df_val_20 <- rbind(df_val_20, df_sample_val_all_15_mi_min)


plot_order_20 <- c("Val-1 20-var", "Val-1 20-var (cc)", "Val-2 20-var", "Val-2 20-var 15U/mL", "Val-3 20-var",  "Val-3 20-var (cc)", "Val-4 20-var", "Val-4 20-var (cc)", "Val-4 20-var 15U/mL", "Val-4 20-var 15U/mL (cc)", "all 20-var", "all 20-var (MI)",  "all 20-var 15U/mL", "all 20-var 15U/mL (MI)")

df_val_20$class <- factor(df_val_20$class, levels = plot_order_20)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

ggplot(data = df_val_20, mapping = aes(x=class, y=auc)) + 
  geom_point(aes(),alpha=0.8) +
  geom_boxplot(fill="grey",color="black",alpha=0.2) + 
  scale_x_discrete(breaks=unique(df_val_20$class), 
                   labels=addline_format(plot_order_20))+
  theme_minimal()

#### Performance on pooled datasets for all 5 models
df_val_all_models <- rbind(df_sample_val_all_lr, df_sample_val_all_min)
df_val_all_models <- rbind(df_val_all_models, df_sample_val_all_1se)
df_val_all_models <- rbind(df_val_all_models, df_sample_val_all_rf)
df_val_all_models <- rbind(df_val_all_models, df_sample_val_all_gbm)
df_val_all_models <- rbind(df_val_all_models, df_sample_val_all_15_lr)
df_val_all_models <- rbind(df_val_all_models, df_sample_val_all_15_min)
df_val_all_models <- rbind(df_val_all_models, df_sample_val_all_15_1se)
df_val_all_models <- rbind(df_val_all_models, df_sample_val_all_15_rf)
df_val_all_models <- rbind(df_val_all_models, df_sample_val_all_15_gbm)

View(df_val_all_models)

plot_order_all_models <- c("all LR", "all 20-var", "all 10-var", "all RF", "all GBRT", "all LR 15U/mL", "all 20-var 15U/mL", "all 10-var 15U/mL", "all RF 15U/mL", "all GBRT 15U/mL")

df_val_all_models$class <- factor(df_val_all_models$class, levels = plot_order_all_models)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}


ggplot(data = df_val_all_models, mapping = aes(x=class, y=auc)) + 
  geom_point(aes(),alpha=0.8) +
  geom_boxplot(fill="grey",color="black",alpha=0.2) + 
  scale_x_discrete(breaks=unique(df_val_all_models$class), 
                   labels=addline_format(plot_order_all_models))+
  theme_minimal()


### Plotting of complete performance data during external validation
df_val_complete <- rbind(df_sample_val_1_1se, df_sample_val_1_cc_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_1_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_1_cc_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_2_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_2_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_2_15_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_2_15_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_3_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_3_cc_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_3_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_3_cc_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_4_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_4_cc_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_4_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_4_cc_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_4_15_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_4_15_cc_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_4_15_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_4_15_cc_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_all_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_all_mi_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_all_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_all_mi_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_all_15_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_all_15_mi_1se)
df_val_complete <- rbind(df_val_complete, df_sample_val_all_15_min)
df_val_complete <- rbind(df_val_complete, df_sample_val_all_15_mi_min)



plot_order_complete <- c("Val-1 10-var", "Val-1 10-var (cc)", "Val-1 20-var", "Val-1 20-var (cc)", "Val-2 10-var",  "Val-2 20-var", "Val-2 10-var 15U/mL", "Val-2 20-var 15U/mL", "Val-3 10-var",  "Val-3 10-var (cc)", "Val-3 20-var",  "Val-3 20-var (cc)", "Val-4 10-var", "Val-4 10-var (cc)", "Val-4 20-var", "Val-4 20-var (cc)", "Val-4 10-var 15U/mL", "Val-4 10-var 15U/mL (cc)", "Val-4 20-var 15U/mL", "Val-4 20-var 15U/mL (cc)", "all 10-var", "all 10-var (MI)",  "all 20-var", "all 20-var (MI)",  "all 10-var 15U/mL", "all 10-var 15U/mL (MI)", "all 20-var 15U/mL", "all 20-var 15U/mL (MI)")
    
df_val_complete$class <- factor(df_val_complete$class, levels = plot_order_complete)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

ggplot(data = df_val_complete, mapping = aes(x=class, y=auc)) + 
  geom_point(aes(),alpha=0.8) +
  geom_boxplot(fill="grey",color="black",alpha=0.2) + 
  scale_x_discrete(breaks=unique(df_val_complete$class), 
                   labels=addline_format(plot_order_complete))+
  theme_minimal()


t.test(df_sample_val_1_min$auc, df_sample_val_1_1se$auc)

### Mean, Median und 95%CI during resampling
quantile(df_sample_val_1_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_1_min$auc)
quantile(df_sample_val_1_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_1_min$sens)
quantile(df_sample_val_1_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_1_min$spec)
quantile(df_sample_val_1_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_1_min$acc)
quantile(df_sample_val_1_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_1_min$ppv)
quantile(df_sample_val_1_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_1_min$npv)

quantile(df_sample_val_1_cc_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_min$auc)
quantile(df_sample_val_1_cc_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_min$sens)
quantile(df_sample_val_1_cc_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_min$spec)
quantile(df_sample_val_1_cc_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_min$acc)
quantile(df_sample_val_1_cc_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_min$ppv)
quantile(df_sample_val_1_cc_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_min$npv)

quantile(df_sample_val_1_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_1_1se$auc)
quantile(df_sample_val_1_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_1_1se$sens)
quantile(df_sample_val_1_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_1_1se$spec)
quantile(df_sample_val_1_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_1_1se$acc)
quantile(df_sample_val_1_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_1_1se$ppv)
quantile(df_sample_val_1_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_1_1se$npv)

quantile(df_sample_val_1_cc_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_1se$auc)
quantile(df_sample_val_1_cc_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_1se$sens)
quantile(df_sample_val_1_cc_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_1se$spec)
quantile(df_sample_val_1_cc_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_1se$acc)
quantile(df_sample_val_1_cc_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_1se$ppv)
quantile(df_sample_val_1_cc_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_1_cc_1se$npv)

quantile(df_sample_val_2_min$auc, prob=c(.025,.5,.975))
quantile(df_sample_val_2_min$sens, prob=c(.025,.5,.975))
quantile(df_sample_val_2_min$spec, prob=c(.025,.5,.975))
quantile(df_sample_val_2_min$acc, prob=c(.025,.5,.975))
quantile(df_sample_val_2_min$ppv, prob=c(.025,.5,.975))
quantile(df_sample_val_2_min$npv, prob=c(.025,.5,.975))

quantile(df_sample_val_2_1se$auc, prob=c(.025,.5,.975))
quantile(df_sample_val_2_1se$sens, prob=c(.025,.5,.975))
quantile(df_sample_val_2_1se$spec, prob=c(.025,.5,.975))
quantile(df_sample_val_2_1se$acc, prob=c(.025,.5,.975))
quantile(df_sample_val_2_1se$ppv, prob=c(.025,.5,.975))
quantile(df_sample_val_2_1se$npv, prob=c(.025,.5,.975))

quantile(df_sample_val_2_15_min$auc, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_min$sens, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_min$spec, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_min$acc, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_min$ppv, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_min$npv, prob=c(.025,.5,.975))

quantile(df_sample_val_2_15_1se$auc, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_1se$sens, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_1se$spec, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_1se$acc, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_1se$ppv, prob=c(.025,.5,.975))
quantile(df_sample_val_2_15_1se$npv, prob=c(.025,.5,.975))

quantile(df_sample_val_3_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_3_min$auc)
quantile(df_sample_val_3_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_3_min$sens)
quantile(df_sample_val_3_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_3_min$spec)
quantile(df_sample_val_3_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_3_min$acc)
quantile(df_sample_val_3_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_3_min$ppv)
quantile(df_sample_val_3_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_3_min$npv)

quantile(df_sample_val_3_cc_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_min$auc)
quantile(df_sample_val_3_cc_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_min$sens)
quantile(df_sample_val_3_cc_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_min$spec)
quantile(df_sample_val_3_cc_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_min$acc)
quantile(df_sample_val_3_cc_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_min$ppv)
quantile(df_sample_val_3_cc_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_min$npv)

quantile(df_sample_val_3_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_3_1se$auc)
quantile(df_sample_val_3_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_3_1se$sens)
quantile(df_sample_val_3_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_3_1se$spec)
quantile(df_sample_val_3_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_3_1se$acc)
quantile(df_sample_val_3_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_3_1se$ppv)
quantile(df_sample_val_3_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_3_1se$npv)

quantile(df_sample_val_3_cc_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_1se$auc)
quantile(df_sample_val_3_cc_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_1se$sens)
quantile(df_sample_val_3_cc_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_1se$spec)
quantile(df_sample_val_3_cc_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_1se$acc)
quantile(df_sample_val_3_cc_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_1se$ppv)
quantile(df_sample_val_3_cc_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_3_cc_1se$npv)

quantile(df_sample_val_4_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_4_min$auc)
quantile(df_sample_val_4_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_4_min$sens)
quantile(df_sample_val_4_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_4_min$spec)
quantile(df_sample_val_4_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_4_min$acc)
quantile(df_sample_val_4_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_4_min$ppv)
quantile(df_sample_val_4_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_4_min$npv)

quantile(df_sample_val_4_cc_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_min$auc)
quantile(df_sample_val_4_cc_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_min$sens)
quantile(df_sample_val_4_cc_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_min$spec)
quantile(df_sample_val_4_cc_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_min$acc)
quantile(df_sample_val_4_cc_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_min$ppv)
quantile(df_sample_val_4_cc_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_min$npv)

quantile(df_sample_val_4_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_4_1se$auc)
quantile(df_sample_val_4_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_4_1se$sens)
quantile(df_sample_val_4_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_4_1se$spec)
quantile(df_sample_val_4_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_4_1se$acc)
quantile(df_sample_val_4_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_4_1se$ppv)
quantile(df_sample_val_4_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_4_1se$npv)

quantile(df_sample_val_4_cc_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_1se$auc)
quantile(df_sample_val_4_cc_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_1se$sens)
quantile(df_sample_val_4_cc_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_1se$spec)
quantile(df_sample_val_4_cc_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_1se$acc)
quantile(df_sample_val_4_cc_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_1se$ppv)
quantile(df_sample_val_4_cc_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_4_cc_1se$npv)

quantile(df_sample_val_4_15_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_min$auc)
quantile(df_sample_val_4_15_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_min$sens)
quantile(df_sample_val_4_15_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_min$spec)
quantile(df_sample_val_4_15_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_min$acc)
quantile(df_sample_val_4_15_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_min$ppv)
quantile(df_sample_val_4_15_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_min$npv)

quantile(df_sample_val_4_15_cc_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_min$auc)
quantile(df_sample_val_4_15_cc_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_min$sens)
quantile(df_sample_val_4_15_cc_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_min$spec)
quantile(df_sample_val_4_15_cc_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_min$acc)
quantile(df_sample_val_4_15_cc_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_min$ppv)
quantile(df_sample_val_4_15_cc_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_min$npv)

quantile(df_sample_val_4_15_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_1se$auc)
quantile(df_sample_val_4_15_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_1se$sens)
quantile(df_sample_val_4_15_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_1se$spec)
quantile(df_sample_val_4_15_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_1se$acc)
quantile(df_sample_val_4_15_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_1se$ppv)
quantile(df_sample_val_4_15_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_1se$npv)

quantile(df_sample_val_4_15_cc_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_1se$auc)
quantile(df_sample_val_4_15_cc_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_1se$sens)
quantile(df_sample_val_4_15_cc_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_1se$spec)
quantile(df_sample_val_4_15_cc_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_1se$acc)
quantile(df_sample_val_4_15_cc_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_1se$ppv)
quantile(df_sample_val_4_15_cc_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_4_15_cc_1se$npv)

quantile(df_sample_val_all_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_min$auc)
quantile(df_sample_val_all_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_min$sens)
quantile(df_sample_val_all_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_min$spec)
quantile(df_sample_val_all_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_min$acc)
quantile(df_sample_val_all_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_min$ppv)
quantile(df_sample_val_all_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_min$npv)

quantile(df_sample_val_all_mi_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_min$auc)
quantile(df_sample_val_all_mi_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_min$sens)
quantile(df_sample_val_all_mi_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_min$spec)
quantile(df_sample_val_all_mi_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_min$acc)
quantile(df_sample_val_all_mi_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_min$ppv)
quantile(df_sample_val_all_mi_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_min$npv)

quantile(df_sample_val_all_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_1se$auc)
quantile(df_sample_val_all_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_1se$sens)
quantile(df_sample_val_all_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_1se$spec)
quantile(df_sample_val_all_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_1se$acc)
quantile(df_sample_val_all_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_1se$ppv)
quantile(df_sample_val_all_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_1se$npv)

quantile(df_sample_val_all_mi_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_1se$auc)
quantile(df_sample_val_all_mi_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_1se$sens)
quantile(df_sample_val_all_mi_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_1se$spec)
quantile(df_sample_val_all_mi_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_1se$acc)
quantile(df_sample_val_all_mi_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_1se$ppv)
quantile(df_sample_val_all_mi_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_mi_1se$npv)

quantile(df_sample_val_all_15_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_min$auc)
quantile(df_sample_val_all_15_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_min$sens)
quantile(df_sample_val_all_15_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_min$spec)
quantile(df_sample_val_all_15_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_min$acc)
quantile(df_sample_val_all_15_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_min$ppv)
quantile(df_sample_val_all_15_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_min$npv)

quantile(df_sample_val_all_15_mi_min$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_min$auc)
quantile(df_sample_val_all_15_mi_min$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_min$sens)
quantile(df_sample_val_all_15_mi_min$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_min$spec)
quantile(df_sample_val_all_15_mi_min$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_min$acc)
quantile(df_sample_val_all_15_mi_min$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_min$ppv)
quantile(df_sample_val_all_15_mi_min$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_min$npv)

quantile(df_sample_val_all_15_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_1se$auc)
quantile(df_sample_val_all_15_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_1se$sens)
quantile(df_sample_val_all_15_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_1se$spec)
quantile(df_sample_val_all_15_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_1se$acc)
quantile(df_sample_val_all_15_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_1se$ppv)
quantile(df_sample_val_all_15_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_1se$npv)

quantile(df_sample_val_all_15_mi_1se$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_1se$auc)
quantile(df_sample_val_all_15_mi_1se$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_1se$sens)
quantile(df_sample_val_all_15_mi_1se$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_1se$spec)
quantile(df_sample_val_all_15_mi_1se$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_1se$acc)
quantile(df_sample_val_all_15_mi_1se$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_1se$ppv)
quantile(df_sample_val_all_15_mi_1se$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_mi_1se$npv)

quantile(df_sample_val_all_lr$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_lr$auc)
quantile(df_sample_val_all_lr$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_lr$sens)
quantile(df_sample_val_all_lr$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_lr$spec)
quantile(df_sample_val_all_lr$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_lr$acc)
quantile(df_sample_val_all_lr$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_lr$ppv)
quantile(df_sample_val_all_lr$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_lr$npv)

quantile(df_sample_val_all_15_lr$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_lr$auc)
quantile(df_sample_val_all_15_lr$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_lr$sens)
quantile(df_sample_val_all_15_lr$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_lr$spec)
quantile(df_sample_val_all_15_lr$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_lr$acc)
quantile(df_sample_val_all_15_lr$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_lr$ppv)
quantile(df_sample_val_all_15_lr$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_lr$npv)

quantile(df_sample_val_all_rf$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_rf$auc)
quantile(df_sample_val_all_rf$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_rf$sens)
quantile(df_sample_val_all_rf$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_rf$spec)
quantile(df_sample_val_all_rf$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_rf$acc)
quantile(df_sample_val_all_rf$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_rf$ppv)
quantile(df_sample_val_all_rf$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_lr$npv)

quantile(df_sample_val_all_15_rf$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_rf$auc)
quantile(df_sample_val_all_15_rf$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_rf$sens)
quantile(df_sample_val_all_15_rf$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_rf$spec)
quantile(df_sample_val_all_15_rf$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_rf$acc)
quantile(df_sample_val_all_15_rf$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_rf$ppv)
quantile(df_sample_val_all_15_rf$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_rf$npv)

quantile(df_sample_val_all_gbm$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_gbm$auc)
quantile(df_sample_val_all_gbm$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_gbm$sens)
quantile(df_sample_val_all_gbm$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_gbm$spec)
quantile(df_sample_val_all_gbm$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_gbm$acc)
quantile(df_sample_val_all_gbm$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_gbm$ppv)
quantile(df_sample_val_all_gbm$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_gbm$npv)

quantile(df_sample_val_all_15_gbm$auc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_gbm$auc)
quantile(df_sample_val_all_15_gbm$sens, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_gbm$sens)
quantile(df_sample_val_all_15_gbm$spec, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_gbm$spec)
quantile(df_sample_val_all_15_gbm$acc, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_gbm$acc)
quantile(df_sample_val_all_15_gbm$ppv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_gbm$ppv)
quantile(df_sample_val_all_15_gbm$npv, prob=c(.025,.5,.975))
mean(df_sample_val_all_15_gbm$npv)

