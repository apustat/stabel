#' Prediction using stabel
#'
#' The \code{stabel.pred} function implements prediction using ensemble learning by combining multiple predictive models to improve overall accuracy and robustness.
#' @param X The input matrix of predictors (features) used to train the model.
#' @param Y The outcome vector corresponding to the input matrix.
#' @param newX A matrix of new samples. If provided, predictions will be made for these new observations. Should have the same column names as \code{X}.
#' @param newY The true outcomes for \code{newX}, used for evaluating the model's performance on test data. If provided, metrics like AUC or MSE can be calculated by comparing predictions to \code{newY}.
#' @param method The method used for combining models in the ensemble. For example, the default method \code{"method.NNLS"} indicates non-negative least squares. Can be changed to other ensemble methods supported by the \code{SuperLearner} package.
#' @param SL.library A list of algorithms/models to be included in the ensemble learning. Check the \code{SuperLearner} package for more details.
#' @param family Currently allows \code{"gaussian"} or \code{"binomial"} to describe the error distribution. Link function information will be ignored and should be contained in the method argument below. Default is \code{"binomial"}.
#' @param nfolds The number of folds used for cross-validation within the ensemble learning. Default is 10.
#' @param thr.prob The threshold probability for classifying predictions into binary classes. Default is 0.5. It should be set to \code{NULL} if family is not binomial.
#' @param use.youden Youden index as threshold value or cutoff point for classifying predictions into binary classes. Default is \code{TRUE}. Both \code{thr.prob} and \code{use.youden} cannot be used together.
#' @param target.specificity The desired specificity at which sensitivity will be evaluated. It is optional.
#' @param target.sensitivity The desired sensitivity at which specificity will be evaluated. It is optional.
#' @param nfolds Number of folds for cross-validation. Default is 10. 
#' @param seed A random seed for reproducibility of results.
#' @keywords stabel
#' @export
#' @return A list containing the following components:
#' \item{method}{Same as above.}
#' \item{family}{Same as above.}
#' \item{nfolds}{Same as above.}
#' \item{cvRisk}{ Cross-validated risk for each candidate algorithm in the SuperLearner ensemble.}
#' \item{coef}{Calculated weight given to each candidate algorithm based on cvRisk.}
#' \item{roc.object}{ROC object used for plotting the ROC curve or computing confidence intervals of the AUC (available for \code{family = "binomial"}).}
#' \item{auc}{Area under the ROC curve (AUC) for classification models (available for \code{family = "binomial"}).}
#' \item{predicted.probs}{Predicted probabilities for binary classification (available for \code{family = "binomial"}).}
#' \item{predicted.responses}{Predicted responses for continuous outcomes (available for \code{family = "gaussian"}).}
#' \item{threshold.probability}{The threshold probability used for classification (available for \code{family = "binomial"}).}
#' \item{confusion.matrix}{Confusion matrix summarizing the classification results (available for \code{family = "binomial"}).}
#' \item{predition.accuracy}{The overall accuracy of the predictions.}
#' \item{sensitivity}{Sensitivity (true positive rate) of the predictions.}
#' \item{specificity}{Specificity (true negative rate) of the predictions.}
#' \item{precision}{Precision (positive predictive value) of the predictions.}
#' \item{negative.predictive.value}{Negative predictive value of the predictions.}
#' \item{positive.predictive.value}{Positive predictive value of the predictions.}
#' \item{false.positve.rate}{False positive rate of the predictions.}
#' \item{false.negative.rate}{False negative rate of the predictions.}
#' \item{F1.score}{F1 score, the harmonic mean of precision and sensitivity.}
#' \item{target.specificity}{User-specified target specificity.}
#' \item{target.sensitivity}{User-specified target sensitivity.}
#' \item{sensitivity.at.specificity}{Sensitivity calculated at the specified target specificity.}
#' \item{specificity.at.sensitivity}{ Specificity calculated at the specified target sensitivity.}
#' \item{mse}{Mean squared error of the predictions (available for \code{family = "gaussian"}).}
#' \item{rmse}{Root mean squared error of the predictions (available for \code{family = "gaussian"}).}
#' \item{mae}{Mean absolute error of the predictions (available for \code{family = "gaussian"}).}
#' @examples
#' stabel.pred()


stabel.pred <- function(X, Y, newX=NULL, newY=NULL, method = "method.NNLS", SL.library, family = "binomial", nfolds= 10, thr.prob=0.5,
                        use.youden=FALSE, target.specificity=NULL, target.sensitivity=NULL, seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Error check: Ensure both thr.prob and use.youden are not used simultaneously
  if (!is.null(thr.prob) && use.youden) {
    stop("Error: Both 'thr.prob' and 'use.youden' cannot be used simultaneously. Set one to NULL or FALSE.")
  }
  # Check if X and newX are matrices
  if (!is.matrix(X)) {
    stop("Error: 'X' must be a matrix.")
  }
  if (!is.null(newX) && !is.matrix(newX)) {
    stop("Error: 'newX' must be a matrix if provided.")
  }
  sl.model = SuperLearner::SuperLearner(Y = Y, X = data.frame(X), newX = data.frame(if (is.null(newX)) X else newX),
                                        method = method, family = family, SL.library, cvControl = list(V=nfolds))
  #sl.pred=predict(sl.model, newX, type="response")
  predicted.values <- as.vector(sl.model$SL.predict)
  
  # Use Y as newY if newY is NULL
  actual.Y <- if (is.null(newY)) Y else newY
  
  if (family == "binomial") {
    ROC <- pROC::roc(as.vector(actual.Y), predicted.values)
    AUC <- auc(ROC)
    
    # Determine threshold probability
    youden.threshold <- if (use.youden) {
      coords(ROC, x = "best", input = "threshold", best.method = "youden")$threshold
    } else {
      NULL
    }
    
    final.threshold <- if (use.youden) {
      youden.threshold
    } else {
      thr.prob
    }
    
    # Ensure a valid threshold is used for binary classification
    threshold <- if (use.youden) youden.threshold else thr.prob
    
    pred.response=ifelse(predicted.values > threshold, 1,0)
    confusion.matrix=table(actual.Y, pred.response)
    TN <- confusion.matrix[1, 1]
    FP <- confusion.matrix[1, 2]
    FN <- confusion.matrix[2, 1]
    TP <- confusion.matrix[2, 2]
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    sensitivity <- TP / (TP + FN) # True Positive Rate
    specificity <- TN / (TN + FP) # True Negative Rate
    precision <- TP / (TP + FP)  # Positive Predictive Value
    npv <- TN / (TN + FN)        # Negative Predictive Value
    fpr <- FP / (FP + TN)        # False Positive Rate
    fnr <- FN / (FN + TP)        # False Negative Rate
    f1.score <- 2 * (precision * sensitivity) / (precision + sensitivity)
    
    # Compute sensitivity at specificity and specificity at sensitivity if defined
    sensitivity_at_specificity <- if (!is.null(target.specificity)) {
      coords(ROC, x = target.specificity, input = "specificity", ret = "sensitivity")$sensitivity
    } else {
      NULL
    }
    
    specificity_at_sensitivity <- if (!is.null(target.sensitivity)) {
      coords(ROC, x = target.sensitivity, input = "sensitivity", ret = "specificity")$specificity
    } else {
      NULL
    }
  } else if (family == "gaussian") {
    # For Gaussian family: predicted.values are numeric predictions
    residuals <- actual.Y - predicted.values
    mse <- mean(residuals^2)
    rmse <- sqrt(mse)
    mae <- mean(abs(residuals))
    
    ROC <- AUC <- conf.int <- NULL
    youden.threshold <- final.threshold <- NULL
    accuracy <- sensitivity <- specificity <- precision <- npv <- NULL
    fpr <- fnr <- f1.score <- sensitivity_at_specificity <- specificity_at_sensitivity <- NULL
  }
  
  return(list(
    method = method,
    family=family,
    nfolds=nfolds,
    cvRisk=sl.model$cvRisk,
    coef=sl.model$coef,
    roc.object=ROC, #if someone want to draw ROC curve or CI of AUC, they would need it
    auc=AUC,
    predicted.probs = if (family == "binomial") as.vector(predicted.values) else NULL,
    predicted.responses = if (family == "gaussian") as.vector(predicted.values) else NULL,
    threshold.probability = if (family == "binomial") (if (use.youden) youden.threshold else thr.prob) else NULL,
    confusion.matrix = if (family == "binomial") confusion.matrix else NULL,
    predition.accuracy=accuracy,
    sensitivity=sensitivity,
    specificity=specificity,
    precision=precision,
    negative.predictive.value=npv,
    positive.predictive.value=npv,
    false.positve.rate=fpr,
    false.negative.rate=fnr,
    F1.score=f1.score,
    target.specificity=target.specificity,
    target.sensitivity=target.sensitivity,
    sensitivity.at.specificity = sensitivity_at_specificity,
    specificity.at.sensitivity = specificity_at_sensitivity,
    mse = if (family == "gaussian") mse else NULL,
    rmse = if (family == "gaussian") rmse else NULL,
    mae = if (family == "gaussian") mae else NULL
  ))
}

