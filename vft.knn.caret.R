library(tidyverse)
library(pROC)
library(caret)
library(nnet)

rm(list = ls())
######## prepare input data 
source("cv_nbmfpca.R")
demo_vft_data = readRDS('demoInfo.rds')
rej.ind.vft = readRDS("rej_ind_vft.rds")
vft_dat = readRDS('vft_data.rds')

ind.mat.vft = matrix(1:nrow(vft_dat), nrow = 52)
#dat.vft = vft_dat[-(as.vector(ind.mat.vft[,33])),]
dat.vft = vft_dat
#5-fold 
set.seed(1)
ind.train = createFolds(y = demo_vft_data$cam, k = 5, returnTrain = T) 

nsubj = 53; nchannels = 52;
ind.mat = matrix(1:(nsubj*nchannels), ncol = nsubj, nrow = nchannels)
rej.mat = matrix(1, nrow = nsubj, ncol = nchannels)
rej.mat[rej.ind.vft] = 0
lvl = 1:6# levels that included in the analysis

scores.vft = list()
cv.vft = list()
######### compute CV scores
for (nset in 1:length(ind.train)){
  scores.vft[[nset]] = cv_nbmfpca(rejection.index = rej.ind.vft, dat = dat.vft, nsubj = 53, nchannels = 52, trainSet = ind.train[[nset]]) 
  print(nset)
}

saveRDS(scores.vft, "CV_VFT.rds")
scores.vft = readRDS("CV_VFT.rds")

######### K-NN Impute
for (nset in 1:length(ind.train)){

  mfpca.fit.train = scores.vft[[nset]]$score.train
  mfpca.fit.test  = scores.vft[[nset]]$score.test
  K1.train = dim(mfpca.fit.train$scores1)[2]
  K1.test  = dim(mfpca.fit.test$scores1)[2]
  K2.train = dim(mfpca.fit.train$scores2)[3]
  K2.test  = dim(mfpca.fit.test$scores2)[3]
  
    # first level scores
  lvl.length.1 = min(K1.train, K1.test, length(lvl))
  lvl.length.2 = min(K2.train, K2.test, length(lvl))
  
  s1.vft.train = mfpca.fit.train$scores1[,1:lvl.length.1]; colnames(s1.vft.train) = paste0("s1.l",1:lvl.length.1)
  s1.vft.test  = mfpca.fit.test$scores1[,1:lvl.length.1];  colnames(s1.vft.test)  = paste0("s1.l",1:lvl.length.1)
  
  # second level scores with imputation
  s2.vft.train = s2.vft.test = list()
  
  for(l in 1:lvl.length.2){ 
    s2.vft.train.impute = mfpca.fit.train$scores2[,,l]; colnames(s2.vft.train.impute) = paste0("C",1:52,"_",l)
    s2.vft.test.impute = mfpca.fit.test$scores2[,,l]; colnames(s2.vft.test.impute) = paste0("C",1:52,"_",l)
    
    set.seed(123);  preProc.s2 = preProcess(s2.vft.train.impute, method = "knnImpute")
    set.seed(123);  s2.vft.train[[l]] = predict(preProc.s2, s2.vft.train.impute)
    set.seed(123);  s2.vft.test[[l]]  = predict(preProc.s2, s2.vft.test.impute)
    }
  
  dat.train = cbind(demo_vft_data[ind.train[[nset]],],s1.vft.train, do.call('cbind', s2.vft.train))
  dat.test  = cbind(demo_vft_data[-ind.train[[nset]],],s1.vft.test, do.call('cbind', s2.vft.test))
  
  cv.vft[[nset]] = list(train = dat.train, test = dat.test)
}


cv.vft.caret = list(index = ind.train, results = cv.vft)
saveRDS(cv.vft.caret, "CV_VFT_caret.rds")
# cv.vft.caret = readRDS("CV_VFT_caret.rds")


### Single hidden layer neural network
vft.result = data.frame(all = rep(NA, nsubj),
                            brain = rep(NA, nsubj), 
                            demo.s1 = rep(NA, nsubj),
                            demo = rep(NA, nsubj))

ctrl <- trainControl(method="repeatedcv", repeats=5, number=5,
                     summaryFunction=twoClassSummary,
                     classProbs=TRUE,savePredictions = "all")

method = "nnet"
ind.train = cv.vft.caret$index
cv.vft.caret.result = cv.vft.caret$results

for (nset in 1:length(ind.train)){
  
  dat.train = cv.vft.caret.result[[nset]]$train
  dat.train$cam = as.factor(dat.train$cam)
  dat.train$gender = as.factor(dat.train$gender)
  levels(dat.train$cam) <- c("Healthy", "Disease")
  
  dat.test = cv.vft.caret.result[[nset]]$test
  dat.test$gender = as.factor(dat.test$gender)
  dat.test$cam = as.factor(dat.test$cam)
  levels(dat.test$cam) <- c("Healthy", "Disease")
  
  dat.train.brain = dat.train[,-(2:10)]
  dat.train.demo = dat.train[,1:10] 
  dat.train.demo.s1 = dat.train[,1:12] 
  dat.test.brain = dat.test[,-(2:10)]
  dat.test.demo = dat.test[,1:10] 
  dat.test.demo.s1 = dat.test[,1:12]
  dat.train = dat.train
  dat.test = dat.test
  
  set.seed(123)
  fit.train  = train(cam ~., data=dat.train, method = method, MaxNWts = 1e4, trace = FALSE,
                        metric= "ROC", trControl=ctrl, preProcess = c("center", "scale"))
  test.pred = predict(fit.train , newdata = dat.test, type = "prob")[,2]
  
  ## train the model use scores data only. 
  set.seed(123)
  fit.train.brain  = train(cam ~., data=dat.train.brain, method = method, MaxNWts = 1e4,trace = FALSE,
                              metric= "ROC", trControl=ctrl, preProcess = c("center", "scale"))
  test.brain.pred = predict(fit.train.brain , dat.test.brain,type = "prob")[,2]
  
  
  ## train the model use demo and level 1 scores data only. 
  set.seed(123)
  fit.train.demo.s1  = train(cam ~., data=dat.train.demo.s1, method = method, MaxNWts = 1e4,trace = FALSE,
                                metric= "ROC", trControl=ctrl, preProcess = c("center", "scale"))
  test.demo.s1.pred = predict(fit.train.demo.s1 , dat.test.demo.s1,type = "prob")[,2]
  
  set.seed(123)
  fit.train.demo  = train(cam ~., data=dat.train.demo, method = method, MaxNWts = 1e4,trace = FALSE,
                                metric= "ROC", trControl=ctrl, preProcess = c("center", "scale"))
  test.demo.pred = predict(fit.train.demo , dat.test.demo,type = "prob")[,2]
  
  vft.result$all[-ind.train[[nset]]] = test.pred
  vft.result$brain[-ind.train[[nset]]] = test.brain.pred
  vft.result$demo.s1[-ind.train[[nset]]] = test.demo.s1.pred
  vft.result$demo[-ind.train[[nset]]] = test.demo.pred
  
  print(c(nset,
          53 - length(ind.train[[nset]]) - length(test.pred),
          53 - length(ind.train[[nset]]) - length(test.brain.pred), 
          53 - length(ind.train[[nset]]) - length(test.demo.s1.pred)
          ))
} 

 saveRDS(vft.result, 'vftCaretResults.rds')
 vft.result1 = vft.result
 vft.result = readRDS('vftCaretResults.rds')
 

par(mfrow = c(1,3))

# png(filename = "vft_scores12.png", width = 1000, height = 800)
rocobj.vft.brain = plot.roc(demo_vft_data$cam[1:53], vft.result$brain,  main="ROC curve (VFT fitted with both scores)", 
                            percent=TRUE,  ci=TRUE,  print.auc=T, print.auc.x = 25, print.auc.y = 1) 
ciobj <- ci.se(rocobj.vft.brain, specificities=seq(0, 100, 5)) # over a select set of specificities
plot(ciobj, type="shape", col="#1c61b6AA")
# dev.off()

# png(filename = "vft_demo_scores1.png", width = 1000, height = 800)
rocobj.vft.demo.s1 <- plot.roc(demo_vft_data$cam[1:53], vft.result$demo.s1,  main="ROC curve (VFT fitted with health info & subject-specific scores)",
                               percent=TRUE,  ci=TRUE,  print.auc=T, print.auc.x = 25, print.auc.y = 1) 
ciobj.vft.demo.s1 <- ci.se(rocobj.vft.demo.s1, specificities=seq(0, 100, 5)) # over a select set of specificities
plot(ciobj.vft.demo.s1, type="shape", col="#1c61b6AA") # plot as a blue shape
# dev.off()

# png(filename = "vft_demo_scores12.png", width = 1000, height = 800)
rocobj.vft = plot.roc(demo_vft_data$cam[1:53], vft.result$all,  main="ROC curve (VFT fitted with health info & both scores)", 
                      percent=TRUE,  ci=TRUE,  print.auc=T, print.auc.x = 25, print.auc.y = 1)
ciobj <- ci.se(rocobj.vft, specificities=seq(0, 100, 5)) # over a select set of specificities
plot(ciobj, type="shape", col="#1c61b6AA") # plot as a blue shape
# dev.off()

# png(filename = "vft_demo.png", width = 1000, height = 800)
rocobj.vft.s1 <- plot.roc(demo_vft_data$cam[1:53], vft.result$demo,  main="ROC curve (VFT data fitted with demographic info)", 
                          percent=TRUE,  ci=TRUE,  print.auc=T, print.auc.x = 15, print.auc.y = 1) 
ciobj.vft.s1 <- ci.se(rocobj.vft.s1, specificities=seq(0, 100, 5)) # over a select set of specificities
plot(ciobj.vft.s1, type="shape", col="#1c61b6AA") # plot as a blue shape
# dev.off()
