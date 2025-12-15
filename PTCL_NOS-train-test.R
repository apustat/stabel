# Load necessary library
library(caret)
library(SIS)
library(stabel)
library(SuperLearner)
library(pROC)
library(glmnet)

df=read.csv("C:/Users/apust/OneDrive - University of Nebraska Medical Center/Desktop/R package paper/PTCL.GEP.DATA/combined.data/PTCL_NOS.csv")[,-1]
##GATA(0)=40, Tbet(1)=59

n <- nrow(df)
p <- ncol(df)
cutoff=0.70

set.seed(199055)
train_idx=sample(1:n, n*0.6, replace=FALSE)
test_idx= -train_idx

train_data <- df[train_idx, ]
test_data <- df[test_idx, ]

X=as.matrix(train_data[, -p])
Y=train_data[, p]

sis_mod<- SIS(X, Y, family='binomial', penalty = "lasso", tune="cv", nsis=30, iter.max =100, iter = FALSE) 
sig_var=sis_mod$sis.ix0 #significant variables by SIS

cv=cv.stabel(X=X[, sig_var],Y=Y,family = "binomial",fit.fun = c("Lasso","sparseSVM"),eval.metric.Lasso = "mse",
    eval.metric.sparseSVM = "me", nfolds = 4, seed=500)

vs=stabel.vs(X=X[, sig_var],Y=Y,cutoff = cutoff,B = 50, bestlam.Lasso = 0.0003, maxit.Lasso = 1e+05,
             bestlam.sparseSVM = 0.68, maxit.sparseSVM=1000, gamma=0.1,
             family = "binomial", dfmax = 10, ntree=300, mcAdj=TRUE, maxRuns = 300, pValue = 0.2,
             fit.fun = c("Lasso", "sparseSVM", "RF"),comb.method = "average", seed = 10)

vs.avg=names(which((rowSums(cbind(vs$selection.probabilities$Lasso, vs$selection.probabilities$sparseSVM, 
                                  vs$selection.probabilities$RF))/3) > cutoff)) #STABEL with average
vs.union= colnames(df)[union(vs$selected.variables$Lasso, union(vs$selected.variables$sparseSVM, vs$selected.variables$RF))] #STABEL with union

pred.avg=stabel.pred(X=X[, vs.avg, drop = FALSE],Y=Y, newX = as.matrix(test_data[, vs.avg, drop = FALSE]),
                     newY = test_data$Response, method = "method.NNLS", SL.library=c("SL.randomForest", "SL.svm","SL.lda"), 
                     family = "binomial",thr.prob=NULL, use.youden = TRUE, nfolds = 4, target.specificity=0.985, seed = 10)


pred.union=stabel.pred(X=X[, vs.union, drop = FALSE],Y=Y,newX = as.matrix(test_data[, vs.union, drop = FALSE]),
                       newY = test_data$Response,method = "method.NNLS", SL.library=c("SL.randomForest", "SL.svm","SL.lda"), 
                       family = "binomial", thr.prob=NULL, use.youden = TRUE, nfolds = 4, target.specificity=0.985, seed = 10)


pred.RF=stabel.pred(X=X[, vs.avg, drop = FALSE],Y=Y, newX = as.matrix(test_data[, vs.avg, drop = FALSE]),
                    newY = test_data$Response, method = "method.NNLS", SL.library="SL.randomForest", 
                    family = "binomial",thr.prob=NULL, use.youden = TRUE, nfolds = 4, target.specificity=0.985, seed = 10)

#######Lasso
set.seed(100)
cvfit.lasso <- cv.glmnet(x=X[, sig_var], y= Y, family = "binomial", alpha=1, type.measure="mse", nfolds=4) #finding the best value of lambda
bestlam.t.lasso <- cvfit.lasso$lambda.min
t.lasso=glmnet(x=X[, sig_var],y=Y, family="binomial", alpha=1, lambda=bestlam.t.lasso, maxit=100000, intercept=TRUE)
#lambda=24*bestlam.t.lasso selects exactly 5 genes same as STABEL
vs.t.lasso=colnames(X[, sig_var])[which(coef(t.lasso)[-1] !=0)] #non-zero coefficients

pred.lasso=stabel.pred(X=X[, vs.t.lasso, drop = FALSE],Y=Y,newX = as.matrix(test_data[, vs.t.lasso, drop = FALSE]),newY = test_data$Response,
                       method = "method.NNLS", SL.library=c("SL.randomForest", "SL.svm","SL.lda"), family = "binomial", 
                       thr.prob=NULL, use.youden = TRUE, nfolds = 4,target.specificity=0.985, seed = 10)

pred.avg$predition.accuracy; pred.RF$predition.accuracy; pred.lasso$predition.accuracy
pred.avg$auc; pred.RF$auc; pred.lasso$auc
pred.avg$sensitivity; pred.RF$sensitivity; pred.lasso$sensitivity
pred.avg$specificity; pred.RF$specificity; pred.lasso$specificity
pred.avg$sensitivity.at.specificity; pred.RF$sensitivity.at.specificity; pred.lasso$sensitivity.at.specificity






stabel.vs(X=X[, sig_var], Y=Y, cutoff = cutoff, B = 50, bestlam.Lasso = 0.0003, 
          maxit.Lasso = 1e+05, bestlam.sparseSVM = 0.68, maxit.sparseSVM = 1000, 
          gamma = 0.1, family = "binomial", dfmax = 10, ntree=300, mcAdj = TRUE, 
          maxRuns = 300, pValue = 0.2, fit.fun = c("Lasso", "sparseSVM", "RF"),
          comb.method = "average", seed = 10)
stabel.pred(X = X[, vs.avg, drop = FALSE], Y=Y, 
            newX = as.matrix(test_data[, vs.avg, drop = FALSE]),
            newY = test_data$Response, method = "method.NNLS", 
            SL.library = c("SL.randomForest", "SL.svm", "SL.lda"), 
            family = "binomial", thr.prob = NULL, use.youden = TRUE, 
            nfolds = 4, target.specificity = 0.985, seed = 10)
