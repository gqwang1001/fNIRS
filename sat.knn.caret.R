library(tidyverse)
library(pROC)
library(caret)
library(nnet)
library(ICSNP)
library(RANN)

rm(list = ls())

source("cv_nbmfpca.R")
demo_sat_data = readRDS('demoInfo.rds')
rej.ind.sat= readRDS("rej_ind.rds")
sat_dat = readRDS('delbox_data.rds')

set.seed(1)
ind.train = createFolds(y = demo_sat_data$cam, k = 5, returnTrain = T) 

nsubj = 53; nchannels = 52;
ind.mat = matrix(1:(nsubj*nchannels), ncol = nsubj, nrow = nchannels)
rej.mat = matrix(1, nrow = nsubj, ncol = nchannels)
rej.mat[rej.ind.sat] = 0

sat.result = data.frame(all = rep(NA, nsubj), brain = rep(NA, nsubj))
sat.result$demo = rep(NA, nsubj)
sat.result$demo.s1 = rep(NA, nsubj)
cv.mfpca.results.sat = vector("list", length = length(ind.train))

scores.sat = cv.sat = list()
seed = 1

#### get cross-validated scores 
for (nset in 1:length(ind.train)){
  scores.sat[[nset]] = cv_nbmfpca(rejection.index = rej.ind.sat, dat = sat_dat, nsubj = 53, 
                                  nchannels = 52, trainSet = ind.train[[nset]]) 
  print(nset)
}

for (nset in 1:length(ind.train)){

  mfpca.fit.train = scores.sat[[nset]]$score.train
  mfpca.fit.test  = scores.sat[[nset]]$score.test
  
  s1.sat.train = mfpca.fit.train$scores1[,1:2]; colnames(s1.sat.train) = c("s1.l1","s1.l2")
  s1.sat.test  = mfpca.fit.test$scores1[,1:2]; colnames(s1.sat.test) = c("s1.l1","s1.l2")
  
  s2.sat.train = cbind(mfpca.fit.train$scores2[,,1], mfpca.fit.train$scores2[,,2])
  s2.sat.test = cbind(mfpca.fit.test$scores2[,,1], mfpca.fit.test$scores2[,,2])
  colnames(s2.sat.train) = c(paste0("C",1:nchannels,"_1"), paste0("C",1:nchannels,"_2"))
  colnames(s2.sat.test) = c(paste0("C",1:nchannels,"_1"), paste0("C",1:nchannels,"_2"))
  
  set.seed(seed);  preProc.s2l1 = preProcess(s2.sat.train[, 1:nchannels], method = "knnImpute")
  set.seed(seed);  preProc.s2l2 = preProcess(s2.sat.train[, -(1:nchannels)], method = "knnImpute")
  
  set.seed(seed);s2.l1.sat.train.imputed = predict(preProc.s2l1, s2.sat.train[, 1:nchannels])
  set.seed(seed);s2.l1.sat.test.imputed  = predict(preProc.s2l1, s2.sat.test[, 1:nchannels])
  set.seed(seed);s2.l2.sat.train.imputed = predict(preProc.s2l2, s2.sat.train[, -(1:nchannels)])
  set.seed(seed);s2.l2.sat.test.imputed  = predict(preProc.s2l2, s2.sat.test[, -(1:nchannels)])
  
  dat.train = cbind(demo_sat_data[ind.train[[nset]],], s1.sat.train, s2.l1.sat.train.imputed, s2.l2.sat.train.imputed)
  dat.test  = cbind(demo_sat_data[-ind.train[[nset]],], s1.sat.test, s2.l1.sat.test.imputed, s2.l2.sat.test.imputed)
  
  # dat.train = cbind(demo_sat_data[ind.train[[nset]],], s1.sat.train, s2.sat.train)
  # dat.test  = cbind(demo_sat_data[-ind.train[[nset]],], s1.sat.test, s2.sat.test)
  
  cv.sat[[nset]] = list(train = dat.train, test = dat.test)
}

cv.sat.caret = list(index = ind.train, results = cv.sat)
#saveRDS(cv.sat.caret, "CV_SAT_caret.rds")

# cv.sat.caret = readRDS("CV_SAT_caret.rds")

ctrl <- trainControl(method = "repeatedcv", repeats = 5, number = 5,
                     summaryFunction = twoClassSummary,
                     classProbs=TRUE, savePredictions = "all")
seed = 123
method = 'nnet'
metric= 'ROC'
preProc = c("center", "scale")


cv.mfpca.results.sat = cv.sat.caret$results
ind.train = cv.sat.caret$index


for (nset in 1:length(ind.train)){
    
  dat.train = cv.mfpca.results.sat[[nset]]$train; 
  dat.train$cam = as.factor(dat.train$cam); 
  dat.train$gender = as.factor(dat.train$gender); 
  levels(dat.train$cam) <- c("Healthy", "Disease")
  
  dat.test = cv.mfpca.results.sat[[nset]]$test; 
  dat.test$cam = as.factor(dat.test$cam); 
  dat.test$gender = as.factor(dat.test$gender);   
  levels(dat.test$cam) <- c("Healthy", "Disease")

  dat.train.brain = dat.train[,-(2:10)]
  dat.train.demo = dat.train[, 1:10] 
  dat.train.demo.s1 = dat.train[, 1:12] 
  dat.test.brain = dat.test[, -(2:10)]
  dat.test.demo = dat.test[, 1:10] 
  dat.test.demo.s1 = dat.test[, 1:12]
  
  set.seed(seed)
  fit.train  = train(cam ~., data=dat.train, method = method, metric= metric, trace = FALSE,
                        trControl=ctrl, preProcess = preProc)
  test.pred = predict(fit.train , newdata = dat.test, type = "prob")[,2]
  
  ## train the model use scores data only. 
  set.seed(seed)
  fit.train.brain  = train(cam ~., data=dat.train.brain, method = method, metric= metric,trace = FALSE,
                              trControl=ctrl, preProcess = preProc)
  test.brain.pred = predict(fit.train.brain , dat.test.brain, type = "prob")[,2]
  
  ## train the model use demo and level 1 scores data only. 

  set.seed(seed)
  fit.train.demo.s1  = train(cam ~., data=dat.train.demo.s1, method = method, metric= metric,trace = FALSE, 
                                trControl=ctrl, preProcess = preProc)
  test.demo.s1.pred = predict(fit.train.demo.s1 , dat.test.demo.s1, type = "prob")[,2]
  
  ## train the model use demo only. 
  
  set.seed(seed)
  fit.train.demo  = train(cam ~., data=dat.train.demo, method = method, metric= metric, trace = FALSE,
                             trControl=ctrl, preProcess = preProc)
  test.demo.pred = predict(fit.train.demo , dat.test.demo, type = "prob")[,2]
  
  sat.result$all[-ind.train[[nset]]] = test.pred
  sat.result$brain[-ind.train[[nset]]] = test.brain.pred
  sat.result$demo.s1[-ind.train[[nset]]] = test.demo.s1.pred
  sat.result$demo[-ind.train[[nset]]] = test.demo.pred
  
  print(nset)
}

saveRDS(sat.result, 'satCaretResults.rds')
sat.result = readRDS('satCaretResults.rds')


par(mfrow = c(1,3))

#png(filename = "sat_scores12.png", width = 1000, height = 800)
rocobj.sat.brain <- plot.roc(demo_sat_data$cam[1:53], sat.result$brain,  main="ROC curve (SAT data fitted with both scores)",
                      percent=TRUE,  ci=TRUE,  print.auc=T, print.auc.x = 30, print.auc.y = 1)
ciobj.sat.brain <- ci.se(rocobj.sat.brain, specificities=seq(0, 100, 5)) # over a select set of specificities  
plot(ciobj.sat.brain, type="shape", col="#1c61b6AA") # plot as a blue shape 
#dev.off()

#png(filename = "sat_demo_scores1.png", width = 1000, height = 800)
rocobj.sat.demo.s1 <- plot.roc(demo_sat_data$cam[1:53], sat.result$demo.s1,  main="ROC curve (SAT fitted with health info and subject-specific scores)",
                               percent=TRUE,  ci=TRUE,  print.auc=T, print.auc.x = 30, print.auc.y = 1)
ciobj.sat.demo.s1 <- ci.se(rocobj.sat.demo.s1, specificities=seq(0, 100, 5)) # over a select set of specificities  
plot(ciobj.sat.demo.s1, type="shape", col="#1c61b6AA") # plot as a blue shape 
# dev.off()

# png(filename = "sat_demo_scores12.png", width = 1000, height = 800)
rocobj.sat <- plot.roc(demo_sat_data$cam[1:53], sat.result$all ,  main="ROC curve (SAT fitted with health info and both scores)",
                       percent=TRUE,  ci=TRUE,  print.auc=T, print.auc.x = 30, print.auc.y = 1) 
ciobj.sat <- ci.se(rocobj.sat, specificities=seq(0, 100, 5)) # over a select set of specificities  
plot(ciobj.sat, type="shape", col="#1c61b6AA") # plot as a blue shape 
#dev.off()

# png(filename = "demo.png", width = 1000, height = 800)
rocobj.sat.demo <- plot.roc(demo_sat_data$cam[1:53], sat.result$demo ,  main="ROC curve (Data fitted with demographic info)",
                            percent=TRUE,  ci=TRUE,  print.auc=T, print.auc.x = 30, print.auc.y = 1)
ciobj.sat.demo <- ci.se(rocobj.sat.demo, specificities=seq(0, 100, 5)) # over a select set of specificities  
plot(ciobj.sat.demo, type="shape", col="#1c61b6AA") # plot as a blue shape 
# dev.off()



