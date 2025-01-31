#' Variable selection using stabel
#'
#' The \code{stabel.vs} function performs variable selection using stability selection in combination with three statistical methods: Lasso, sparseSVM, and Random Forest (RF). It implements the stability selection framework proposed by Shah and Samworth (2013), considering no assumption, to identify robust subsets of variables. The function then combines the selected subsets from these methods through ensemble learning. One can use any individual models or any pair or all three.
#' @param X Input matrix.
#' @param Y Outcome vector.
#' @param cutoff Threshold of selection probability. Default is 0.6.
#' @param B Number of subsamples. Total 2B number of subsamples. Default is 50.
#' @param bestlam.Lasso Cross-validated lambda for Lasso. If fit.fun includes \code{"Lasso"}, then \code{bestlam.Lasso} must be specified,  otherwise it will throw an error.
#' @param maxit.Lasso Maximum number of passes over the data for all lambda values. Default is \eqn{10^5}
#' @param family A description of the error distribution and link function to be used in the model. Even though Lasso and RF work for both \code{"gaussian"} and \code{"binomial"}, sparseSVM works only for \code{"binomial"} family. Default is \code{"binomial"}.
#' @param dfmax Limit the maximum number of variables in the model for both Lasso and sparseSVM. Default is p+1, p is the number of predictors.
#' @param maxit.sparseSVM Maximum number of iterations for sparseSVM. Default is 1000.
#' @param gamma  The tuning parameter for huberization smoothing of hinge loss. Default is 0.1.
#' @param bestlam.sparseSVM Cross-validated lambda for sparseSVM. If fit.fun includes \code{"sparseSVM"}, then \code{bestlam.sparseSVM} must be specified; otherwise it will throw an error.
#' @param ntree Number of trees in the forest.
#' @param mcAdj If set to TRUE, a multiple comparisons adjustment using the Bonferroni method will be applied. Default is \code{TRUE}.
#' @param maxRuns Number of importance source runs. You may increase it to resolve attributes left Tentative. Default is 100.
#' @param pValue Confidence level for RF. Default is 0.01.
#' @param fit.fun Allows users to perform variable selection using stability selection with \code{"Lasso"}, \code{"sparseSVM"}, \code{"RF"} or any combination of these methods. Default is \code{c("Lasso", "sparseSVM")}.
#' @param comb.method For combining two or variable selection methods. Currently allows either \code{"average"} which takes an average of the empirical selection probability or \code{"union"} which takes a union of the selected subsets. If \code{fit.fun} includes only one method, then \code{comb.method} should be specified as \code{NULL}, otherwise it will throw an error.
#' @param seed A random seed for reproducibility of results.
#' @keywords stabel
#' @export
#' @return A list containing the following components:
#' \item{selection.probabilities}{A list of multiple vectors of selection probabilities for each variable, indicating their importance across subsamples.}
#' \item{selected.variables}{A list of multiple vectors of selected variables based on the defined cutoff criteria.}
#' \item{final.selected.set}{The final set of variables selected after combining results from multiple methods.}
#' \item{cutoff}{Same as above.}
#' \item{B}{Same as above.}
#' \item{bestlam.Lasso}{The optimal regularization parameter (\eqn{\lambda}) for the Lasso model, determined via cross-validation.}
#' \item{maxit.Lasso}{Same as above.}
#' \item{family}{Same as above.}
#' \item{maxit.sparseSVM}{Same as above.}
#' \item{dfmax}{Same as above.}
#' \item{gamma}{Same as above.}
#' \item{bestlam.sparseSVM}{The optimal regularization parameter (\eqn{\lambda}) for the sparseSVM model, determined via cross-validation.}
#' \item{ntree}{Same as above.}
#' \item{mcAdj}{Same as above.}
#' \item{maxRuns}{Same as above.}
#' \item{pValue}{Same as above.}
#' \item{fit.fun}{Same as above.}
#' \item{comb.method}{Same as above.}
#' \item{seed}{Same as above.}
#' @examples
#' stabel.vs()



stabel.vs=function(X, Y, cutoff=0.6, B=50, bestlam.Lasso=NULL, maxit.Lasso=1e+05, family="binomial", maxit.sparseSVM=1000, 
                   dfmax=p+1, gamma=NULL, bestlam.sparseSVM=NULL, ntree = 200, mcAdj = TRUE, maxRuns=100, pValue=0.01, 
                   fit.fun=c("Lasso", "sparseSVM"), comb.method="average", seed=NULL){
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (!is.matrix(X)) {
    stop("Error: 'X' must be a matrix.")
  }
  
  p <- ncol(X)
  n <- nrow(X)
  variable_names <- colnames(X)
  
  valid_methods <- c("Lasso", "sparseSVM", "RF")
  
  valid_families <- c("binomial", "gaussian")
  
  # Check if family is NULL
  if (!is.null(family)) {
    # Validate the family argument
    if (!family %in% valid_families) {
      stop("Error: 'family' must be one of ", paste(valid_families, collapse = ", "), " or NULL")
    }
  } else {
    # If family is NULL, keep it as NULL
    family <- NULL
  }
  
  
  # Check if all methods in fit.fun are valid
  if (!all(fit.fun %in% valid_methods)) {
    stop("Invalid input in fit.fun: All methods must be one or more of the following valid options: 'Lasso', 'sparseSVM', 'RF'.")
  }
  
  # Proceed if all methods are valid
  message("Selected methods: ", paste(fit.fun, collapse = ", "))
  
  # Check for incompatible combinations and adjust fit.fun
  if (!is.null(family) && family == "gaussian" && "sparseSVM" %in% fit.fun) {
    message("Warning: 'sparseSVM' is not supported for family='gaussian'. It will be dropped from the algorithm.")
    fit.fun <- setdiff(fit.fun, "sparseSVM")
  }
  
  if ("Lasso" %in% fit.fun) {
    if (is.null(bestlam.Lasso)) stop("Error: 'bestlam.Lasso' must be specified for fit.fun = 'Lasso'")
    if (is.null(family)) stop("Error: 'family' must be specified for fit.fun = 'Lasso'")
  }
  
  if ("sparseSVM" %in% fit.fun) {
    if (is.null(bestlam.sparseSVM)) stop("Error: 'bestlam.sparseSVM' must be specified for fit.fun = 'sparseSVM'")
    if (is.null(gamma)) stop("Error: 'gamma' must be specified for fit.fun = 'sparseSVM'")
  }
  
  if ("RF" %in% fit.fun) {
    if (is.null(maxRuns)) stop("Error: 'maxRuns' must be specified for fit.fun = 'RF'")
    if (is.null(pValue)) stop("Error: 'pValue' must be specified for fit.fun = 'RF'")
  }
  
  # Valid values for comb.method
  valid_comb_methods <- c("average", "union")
  
  # Check the number of methods in fit.fun
  if (length(fit.fun) == 1) {
    if (!is.null(comb.method)) {
      stop("Error: 'comb.method' must be NULL when 'fit.fun' includes only one method.")
    }
    comb.method <- NULL
  } else if (length(fit.fun) > 1) {
    if (is.null(comb.method)) {
      stop("Error: 'comb.method' must be specified when 'fit.fun' includes more than one method.")
    }
    comb.method <- match.arg(comb.method, choices = valid_comb_methods)  # Validate comb.method
  } else {
    stop("Error: 'fit.fun' must include at least one method.")
  }
  
  selected.variables <- list()
  selection.probabilities <- list()
  
  # Apply Lasso if chosen
  if ("Lasso" %in% fit.fun) {
    sel.var.lasso.vec <- c() 
    for (i in 1:B){
      set.seed(seed+i)
      CPSS.lasso <- sample(1:n, n*0.5, replace=FALSE) #subsampling the training dataset
      yy.lasso1 <- Y[CPSS.lasso]
      xx.lasso1 <- as.matrix(X[CPSS.lasso, ])
      #lasso.mod1= cv.glmnet(x=xx.lasso1, y=yy.lasso1, family=family, alpha=1, maxit=maxit.Lasso, dfmax=dfmax) 
      lasso.mod1 <- glmnet::glmnet(x=xx.lasso1, y=yy.lasso1, family=family, alpha=1, lambda=bestlam.Lasso, maxit=maxit.Lasso, intercept=TRUE) 
      sel.var.lasso1 <- which(coef(lasso.mod1)[-1] !=0) #non-zero coefficients
      
      yy.lasso2 <- Y[-CPSS.lasso]
      xx.lasso2 <- as.matrix(X[-CPSS.lasso, ])
      #lasso.mod2= cv.glmnet(x=xx.lasso2, y=yy.lasso2, family=family, alpha=1, maxit=maxit.Lasso, dfmax=dfmax) 
      lasso.mod2 <- glmnet::glmnet(x=xx.lasso2, y=yy.lasso2, family=family, alpha=1, lambda=bestlam.Lasso, maxit=maxit.Lasso, intercept=TRUE)
      sel.var.lasso2 <- which(coef(lasso.mod2)[-1] !=0) #non-zero coefficients
      sel.var.lasso <- c(sel.var.lasso1, sel.var.lasso2)
      sel.var.lasso.vec <- c(sel.var.lasso.vec, sel.var.lasso) # keeping the record of selected variables (column)
    }
    
    sel.var.lasso.tab <- table(factor(sel.var.lasso.vec, levels = 1:p))
    sel.var.lasso.df <- as.data.frame(sel.var.lasso.tab)
    sel.prob.lasso <- sel.var.lasso.df$Freq/(2*B)
    names(sel.prob.lasso) <- variable_names
    sel.lasso <- unname(which(sel.prob.lasso > cutoff))
    selected.variables$Lasso <- sel.lasso
    selection.probabilities$Lasso <- sel.prob.lasso
  }
  
  ##############stability selection with Sparse SVM
  #y=Output vector. Currently the function only supports binary output and converts the output into +1/-1 coding internally
    
    if ("sparseSVM" %in% fit.fun) {
    sel.var.svm.vec <- c() 
    for (j in 1:B){
      set.seed(seed+j)
      CPSS.svm <- sample(1:n, n*0.5, replace=FALSE) # subsampling the training dataset
      yy.svm1 <- Y[CPSS.svm]
      xx.svm1 <- as.matrix(X[CPSS.svm, ])
      svm.mod1 <- sparseSVM::sparseSVM(X=xx.svm1, y=yy.svm1, alpha=1, gamma=gamma, lambda=bestlam.sparseSVM, max.iter=maxit.sparseSVM, dfmax=dfmax) #dfmax=p+1, Upper bound for the number of nonzero coefficients.
      sel.var.svm1 <- which(coef(svm.mod1)[-1] !=0) #all estimated coefficients, excluding intercept
      
      yy.svm2 <- Y[-CPSS.svm]
      xx.svm2 <- as.matrix(X[-CPSS.svm, ])
      svm.mod2 <- sparseSVM::sparseSVM(X=xx.svm2, y=yy.svm2, alpha=1, gamma=gamma, lambda=bestlam.sparseSVM, max.iter=maxit.sparseSVM, dfmax=dfmax)
      sel.var.svm2 <- which(coef(svm.mod2)[-1] !=0) #all estimated coefficients, excluding intercept
      sel.var.svm <- c(sel.var.svm1, sel.var.svm2)
      sel.var.svm.vec <- c(sel.var.svm.vec, sel.var.svm) # keeping the record of selected variables (column)
    }
    
    sel.var.svm.tab <- table(factor(sel.var.svm.vec, levels = 1:p))
    sel.var.svm.df <- as.data.frame(sel.var.svm.tab)
    sel.prob.svm <- sel.var.svm.df$Freq/(2*B)
    names(sel.prob.svm) <- variable_names
    sel.svm <- unname(which(sel.prob.svm > cutoff))
    selected.variables$sparseSVM <- sel.svm
    selection.probabilities$sparseSVM <- sel.prob.svm
  }
  
  ##############stability selection with RF
  if ("RF" %in% fit.fun) {
    sel.var.rf.vec <- c()
    for (k in 1:B){
      set.seed(seed+k)
      CPSS.rf <- sample(1:n, n*0.5, replace=FALSE) # subsampling the training dataset
      yy.rf1 <- Y[CPSS.rf]
      xx.rf1 <- as.matrix(X[CPSS.rf, ])
      rf.mod1 <- Boruta::Boruta(x=xx.rf1, y=yy.rf1,  doTrace = 0, ntree = ntree, mcAdj = mcAdj, maxRuns = maxRuns, pValue = pValue)
      sel.var.rf1 <- which(rf.mod1$finalDecision == "Confirmed") # keeping the record of selected variables (column)
      
      yy.rf2 <- Y[-CPSS.rf]
      xx.rf2 <- as.matrix(X[-CPSS.rf, ])
      rf.mod2 <- Boruta::Boruta(x=xx.rf2, y=yy.rf2, doTrace = 0, ntree = ntree, mcAdj = mcAdj, maxRuns = maxRuns, pValue = pValue)
      sel.var.rf2 <- which(rf.mod2$finalDecision == "Confirmed") # keeping the record of selected variables (column)
      
      sel.var.rf <- c(sel.var.rf1, sel.var.rf2)
      sel.var.rf.vec <- c(sel.var.rf.vec, sel.var.rf) # combining the selected variables (columns) from he above two models
    }
    
    sel.var.rf.tab <- table(factor(sel.var.rf.vec, levels = 1:p))
    sel.var.rf.df <- as.data.frame(sel.var.rf.tab)
    sel.prob.rf <- sel.var.rf.df$Freq/(2*B)
    names(sel.prob.rf) <- variable_names
    sel.rf <- unname(which(sel.prob.rf > cutoff))
    selected.variables$RF <- sel.rf
    selection.probabilities$RF <- sel.prob.rf
  }
  
  # Combine results based on comb.method
  if (length(fit.fun) == 1) {
    final.selected.set <- selected.variables[[fit.fun]]
  } else {
    if(comb.method =="union"){
      final.selected.set <- sort(union(union(selected.variables$Lasso, selected.variables$sparseSVM), selected.variables$RF)) 
    } else if (comb.method == "average"){
      final.selected.set <- unname(which((rowSums(cbind(selection.probabilities$Lasso, selection.probabilities$sparseSVM, selection.probabilities$RF))/
                                            length(selection.probabilities)) > cutoff))
    } else {
      stop("Invalid comb.method. Please use 'union' or 'average'.")
    }
  }
  
  # Return a list of outputs
  return(list(
    selection.probabilities = selection.probabilities,
    selected.variables = selected.variables,
    final.selected.set = final.selected.set,
    cutoff = cutoff,
    B = B,
    bestlam.Lasso = bestlam.Lasso,
    maxit.Lasso=maxit.Lasso,
    family = family,
    maxit.sparseSVM = maxit.sparseSVM,
    dfmax = dfmax,
    gamma = gamma,
    bestlam.sparseSVM = bestlam.sparseSVM,
    ntree=ntree,
    mcAdj=mcAdj,
    maxRuns = maxRuns,
    pValue = pValue,
    fit.fun = fit.fun,
    comb.method = comb.method,
    seed = seed
  ))
}    

              