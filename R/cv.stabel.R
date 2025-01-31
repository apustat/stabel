#' Cross-validation for stabel
#'
#' The \code{cv.stabel} function performs cross-validation for stabel, allowing the use of either Lasso, sparseSVM, or both for variable selection and model evaluation.
#' @param X Input matrix.
#' @param Y Outcome vector.
#' @param family A description of the error distribution and link function to be used in the model. Even though Lasso works for both \code{"gaussian"} and \code{"binomial"} family, \code{sparseSVM} works only for \code{"binomial"} family. Default is \code{"binomial"}.
#' @param fit.fun List the name of the planned variable selection methods. Default is \code{c("Lasso", "sparseSVM")}.  It can be changed to either only \code{"Lasso"} or \code{"sparseSVM"}.
#' @param eval.metric.Lasso Loss to use for cross-validation. Currently five options, not all available for all models. Check \code{glmnet} package for more details. If \code{fit.fun} includes \code{"Lasso"}, then \code{eval.metric.Lasso} must be specified,  otherwise it will throw an error.
#' @param eval.metric.sparseSVM The metric used to choose optimal lambda. Current version only supports "me": misclassification error. If fit.fun includes \code{"sparseSVM"}, then \code{eval.metric.sparseSVM} must be specified,  otherwise it will throw an error.
#' @param nfolds Number of folds for cross-validation. Default is 10. 
#' @param seed A random seed for reproducibility of results.
#' @keywords stabel
#' @export
#' @return A list containing the following components:
#' \item{bestlam.Lasso}{The optimal regularization parameter (\eqn{\lambda}) for the Lasso model selected through cross-validation.}
#' \item{bestlam.sparseSVM}{The optimal regularization parameter (\eqn{\lambda}) for the sparseSVM model selected through cross-validation.}
#' \item{family}{Same as above.}
#' \item{eval.metric.Lasso}{The evaluation metric used to determine the optimal model for Lasso.}
#' \item{eval.metric.sparseSVM}{The evaluation metric used to determine the optimal model for sparseSVM.}
#' \item{nfolds}{Same as above.}
#' @examples
#' cv.stabel()

cv.stabel <- function(X, Y, family="binomial", fit.fun=c("Lasso", "sparseSVM"), eval.metric.Lasso=NULL, eval.metric.sparseSVM=NULL, nfolds=10, seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (!is.matrix(X)) {
    stop("Error: 'X' must be a matrix.")
  }
  
  valid_methods <- c("Lasso", "sparseSVM")
  if (!all(fit.fun %in% valid_methods)) {
    stop("Invalid fit.fun specified. Choose from 'Lasso', 'sparseSVM', or both.")
  }
  
  # Validate family compatibility
  if ("Lasso" %in% fit.fun && !family %in% c("binomial", "gaussian")) {
    stop("Lasso supports only family='binomial' or 'gaussian'.")
  }
  if ("sparseSVM" %in% fit.fun && family != "binomial") {
    stop("sparseSVM supports only family='binomial'.")
  }
  
  # Validate eval.metric for each method
  if ("Lasso" %in% fit.fun && is.null(eval.metric.Lasso)) {
    stop("eval.metric.Lasso must be specified when fit.fun includes 'Lasso'.")
  }
  if ("sparseSVM" %in% fit.fun && is.null(eval.metric.sparseSVM)) {
    stop("eval.metric.sparseSVM must be specified when fit.fun includes 'sparseSVM'.")
  }
  
  # Initialize results
  bestlam.Lasso <- NULL
  bestlam.sparseSVM <- NULL
  # Apply Lasso if chosen
  if ("Lasso" %in% fit.fun) {
    cvfit.Lasso <- glmnet::cv.glmnet(x=as.matrix(X), y= Y, family = family, alpha=1, type.measure=eval.metric.Lasso, nfolds=nfolds) #finding the best value of lambda
    bestlam.Lasso <- cvfit.Lasso$lambda.min
  }
  if ("sparseSVM" %in% fit.fun) {
    cvfit.sparseSVM <- sparseSVM::cv.sparseSVM(X=as.matrix(X), y=Y, ncores = 1, eval.metric = eval.metric.sparseSVM, nfolds = nfolds, seed=seed, trace = FALSE)
    bestlam.sparseSVM <- cvfit.sparseSVM$lambda.min
  }
  return(list(
    bestlam.Lasso=bestlam.Lasso,
    bestlam.sparseSVM=bestlam.sparseSVM,
    family=family,
    eval.metric.Lasso=eval.metric.Lasso,
    eval.metric.sparseSVM=eval.metric.sparseSVM,
    nfolds=nfolds
  ))
}

#bestlam= cv.stabel(X=X, Y=Y, family="binomial", eval.metric.sparseSVM = "me",nfolds=5, seed=100)
