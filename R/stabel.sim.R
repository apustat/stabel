#' Simulating data for stabel
#'
#' The \code{stabel.sim} function generates a simulated dataset with multivariate normal covariates and corresponding outcomes based on specified model parameters. The simulation supports both continuous and binary outcomes and allows for flexible correlation structures.
#' @param n The total number of observations to be generated in the dataset.
#' @param p The number of covariates (or predictors) in the dataset.
#' @param beta A vector of coefficients for the covariates, representing the true relationship between the predictors and the outcome variable. The length should be same as p.
#' @param family Specifies the type of outcome variable to be generated. Options include \code{"gaussian"} for continuous outcomes and \code{"binary"} for binary outcomes. Default is \code{"binomial"}.
#' @param mu A vector of means for the covariates. It specifies the center of the multivariate normal distribution used for generating covariates. Default is a vector of zero.
#' @param rho The correlation parameter that controls strength of the relationship between covariates. Default is \code{NULL} as the default \code{corstr} is \code{"Independent"}.
#' @param corstr The correlation structure to use for generating the covariates. Options include: \code{"Independent"}, \code{"Exchangeable"}, and \code{"AR1"}. Default is \code{"Independent"}.
#' @param custom.corr A user-specified correlation matrix with dimension \code{p x p}. Default is \code{NULL}.
#' @param seed A random seed for reproducibility of results.
#' @keywords stabel
#' @export
#' @return A list containing the following components:
#' \item{X}{A matrix of size \code{n x p} containing the simulated covariates.}
#' \item{Y}{A vector of length \code{n} containing the simulated outcomes. Continuous outcomes for \code{family = "gaussian"} and binary outcomes (0 or 1) for \code{family = "binary"}.}
#' @examples
#' stabel.sim()
#' 
#n, number of samples
#p, number of predictors
#beta, true values of perameters
#outcome, type of outcome, default is continuous
#mu, mean vector of the design matrix, default is a vector of 0 with length of p
#rho, strength  of correlation, default is 0.5
#corstr, correlation structure, default is Independent. Currently captures Independent, Exchangeable, AR1, or Toeplitz
#seed, for reproducibility, default is NULL

stabel.sim <- function(n, p, beta, family = "binomial", mu = rep(0, p), rho = NULL, 
                       corstr = "Independent", custom.corr = NULL, seed = NULL) {
  # Set the seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Validate conflicting input scenarios
  if (!is.null(custom.corr) && (!is.null(rho) || !is.null(corstr))) {
    stop("Both 'custom.corr' and either 'corstr' or 'rho' can not be specified. Please specify only one correlation input.")
  }
  
  if (is.null(custom.corr)) {
    if (is.null(corstr) || !corstr %in% c("Independent", "Exchangeable", "AR1")) {
      stop("Invalid 'corstr'. Choose from 'Independent', 'Exchangeable', or 'AR1'.")
    }
    if (corstr != "Independent" && is.null(rho)) {
      stop("The 'rho' parameter must be provided for correlation structures other than 'Independent'.")
    }
  }
  
  if (!is.null(custom.corr)) {
    # Validate custom correlation matrix
    if (!is.matrix(custom.corr) || nrow(custom.corr) != p || ncol(custom.corr) != p) {
      stop("'custom.corr' must be a square matrix with dimensions equal to the number of covariates (p).")
    }
    Sigma <- custom.corr
  } else {
    # Generate the covariance matrix based on the correlation structure
    Sigma <- switch(corstr,
                    "Independent" = diag(p),
                    "Exchangeable" = matrix(rho, nrow = p, ncol = p) + diag(1 - rho, p),
                    "AR1" = outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    )
  }
  
  # Generate multivariate normal data for X
  X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  
  # Generate error term epsilon for continuous outcomes
  epsilon <- rnorm(n, mean = 0, sd = 1)
  
  # Generate Y based on the type of outcome
  if (family == "gaussian") {
    Y <- X %*% beta + epsilon  # Linear model: Y = Xb + e
  } else if (family == "binomial") {
    probs <- (exp(X %*% beta)) / (1 + exp(X %*% beta))  # Logistic transformation
    Y <- rbinom(n, size = 1, prob = probs)  # Generate binary outcomes
  } else {
    stop("Invalid family type. Choose 'gaussian' or 'binomial'.")
  }
  
  # Return a list with X and Y
  return(list(X = X, Y = Y))
}
