
## Installation
```r
install.packages("devtools")
library(devtools)
devtools::install_github("apustat/stabel")
```

## Loading all dependent libraries
```{r, warning=FALSE, echo=TRUE, warning=FALSE, message=FALSE, results='hide'}
library(stabel)
library(MASS)
library(glmnet)
library(Boruta)
library(sparseSVM)
library(SuperLearner)
library(pROC)
```

## Background
In this document, we introduce STABEL, a variable selection and prediction approach that combines stability selection and ensemble learning to identify predictive biomarkers in both low- and high-dimensional datasets. The stability selection framework ensures robustness by controlling the family-wise error rate, while the ensemble learning approach integrates multiple algorithms to improve predictive accuracy. We refer to this method as Stability Selection with Ensemble Learning (STABEL). To evaluate and compare model performance, we provide functionality for cross-validation and detailed metrics, enabling users to assess predictive performance comprehensively across various datasets. We use a simulated dataset generated by $stabel.sim$ function from our package $stabel$. 

## Data Generation

The function $stabel.sim$ generates a synthetic dataset from a multivariate normal distribution under a specified correlation structure. In this example, the dataset consists of $n=200$ observations and $p=50$ predictors. The true coefficients $(β)$ are defined such that the first 5 predictors have an effect size of 1, while the remaining 45 predictors have no effect (coefficients set to 0). The response variable is generated under a binomial family $(family="binomial")$, with a mean vector $(\mu)$ of zeros for all 50 predictors. The correlation structure among the predictors is specified as independent $(corstr="Independent")$, meaning there is no correlation among the predictors. Currently, the function supports only independent, exchangeable, and AR(1) correlation structures. Users can also use their customized correlation matrix by specifying the correlation matrix through the custom.corr argument. However, both $custom.corr$ and either $corstr$ or $rho$ can not be specified. Additionally, the function allows for reproducibility by setting a random seed. Please note that this $cv.stabel$ is not required for further analysis. You can use your own data if you have any.

```{r, echo=TRUE, warning=FALSE, message=FALSE}
n=200 #number of samples
p=50 #number of predictors
data <- stabel.sim(n=n,p=p, beta=c(rep(1,5), rep(0, p-5)), family = "binomial", mu = rep(0, p),
                  rho = NULL, corstr = "Independent", custom.corr = NULL, seed = 1)
X <- data$X
Y <- data$Y
dim(X)
head(X[,1:5])
head(Y)
```

## Cross-validation
We begin by dividing the dataset into training and testing subsets, where the training set comprises 70% of the data, leaving the remaining 30% as the testing set. This partition ensures that the training data is used to develop the models, while the testing data is reserved for evaluating their performance. For model tuning, we utilize the $cv.stabel$ function, which conducts cross-validation specifically for the Lasso and sparseSVM methods. The cross-validation process iteratively evaluates model performance across a range of regularization parameter $(\lambda)$ values, ultimately identifying the optimal $\lambda$ for each method. It is important to note that the Random Forest (RF) approach does not require cross-validation for tuning since it does not rely on a regularization parameter like Lasso and sparseSVM.

```{r}
train <- sample(1:dim(X)[1], dim(X)[1]*0.7, replace=FALSE) #70% data in the training set
train_X <- X[train,]
train_Y <- Y[train]
test_X <- X[-train,]
test_Y <- Y[-train]
cv <- cv.stabel(X=train_X, Y=train_Y, family="binomial", fit.fun=c("Lasso", "sparseSVM"), 
          eval.metric.Lasso="mse", eval.metric.sparseSVM="me", nfolds=5, seed=1)

```

## Variable Selection
After determining the optimal values of the tuning parameters, we can apply our STABEL method for variable selection using the $stabel.vs$ function. It is important to note that variable selection will be performed exclusively on the training data. In the following code, we combine Lasso and sparseSVM to obtain the final set of selected variables. This function currently supports combination methods: $"average"$, which calculates the average of the empirical selection probabilities, and $"union"$, which combines the selected subsets by taking both common and uncommon. As demonstrated, stability selection with Lasso identifies 4 variables, missing one true signal, while stability selection with sparseSVM selects 6 variables, missing one true signal and including 2 noise variables. In contrast, STABEL successfully selects all 5 true signals without including any noise variables. This highlights its advantage, particularly when the objective is to identify only the true biomarkers associated with a specific disease.

```{r, warning=FALSE, message=FALSE}
bestlam.Lasso <- cv$bestlam.Lasso
bestlam.sparseSVM <- cv$bestlam.sparseSVM
vs <- stabel.vs(X=train_X, Y=train_Y, cutoff=0.90, B=50, bestlam.Lasso=bestlam.Lasso, maxit=1e+05, 
          family="binomial", max.iter=1000, dfmax=51, gamma=0.1, 
          bestlam.sparseSVM=bestlam.sparseSVM, maxRuns=100, pValue=0.01, 
          fit.fun=c("Lasso", "sparseSVM"), comb.method="average", seed=1)
vs$selected.variables$Lasso
vs$selected.variables$sparseSVM
vs$final.selected.set
```

The $stabel.vs$ function provides an easy way to obtain the selection probabilities for each method as well as their combined probabilities. These probabilities can then be visualized through straightforward plotting. The resulting plot clearly demonstrates that the true signals have the highest selection probabilities.

```{r, warning=FALSE, message=FALSE}
avg.sel.prob=rowSums(cbind(vs$selection.probabilities$Lasso, vs$selection.probabilities$sparseSVM))/2
plot(avg.sel.prob, xlab="Variable index", ylab="Empirical selection probability")
```

## Prediction

The $stabel.pred$ function facilitates prediction using the superlearner algorithm, which integrates multiple base learners to enhance predictive performance. This function takes the selected variables from the $stabel.vs$ function and applies the superlearner framework to build a robust prediction model. By combining the strengths of different algorithms, the superlearner ensures optimal predictive accuracy, leveraging cross-validation to assign appropriate weights to each learner. Users can customize the set of base learners to suit their specific needs, making the $stabel.pred$ function a flexible and powerful tool for prediction tasks, especially when dealing with high-dimensional data or complex variable relationships.

```{r, warning=FALSE, message=FALSE}
pred <- stabel.pred(X=train_X[, vs$final.selected.set], Y=train_Y, newX=test_X[, vs$final.selected.set], 
                    newY=test_Y, method = "method.NNLS", family = "binomial", nfolds= 5, thr.prob=0.22,
                    SL.library=c("SL.randomForest", "SL.glm", "SL.svm", "SL.lda"), 
                    target.specificity=NULL, target.sensitivity=NULL, seed=1)
#coords(pred$roc.object, x="best", input="threshold", best.method="youden")
plot(pred$roc.object, legacy.axes=TRUE, col = "black", xlab="1-Specificity", lty=1,type="l", lwd=1,
     cex.lab=1, main="ROC curve")
```

The $stabel.pred$ function provides comprehensive prediction performance metrics tailored to the type of outcome. For binary outcomes, it calculates metrics such as confusion matrix, sensitivity, specificity, prediction accuracy, AUC, F1-score, precision, positive predictive value, and negative predictive value. For continuous outcomes, it offers MSE, RMSE, and MAE. It is also possible to obtain either predicted probabilities or predicted responses for further use. Additionally, the function enables the plotting of an ROC curve using the $roc.object$ generated during prediction. While the default threshold probability for categorizing binary outcomes is set to 0.5, users can customize this threshold, incorporating values derived from methods like the Youden index for improved classification.
