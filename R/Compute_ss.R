#' Perform Univariate Linear Regression Separately for Columns of X
#' 
#' This is a function provided in the package of "susieR", Wang et al (2020) <doi:10.1101/501114>,
#' for performing the univariate linear regression y ~ x separately for each column x of X to generate 
#' summary statistics. Each regression is implemented using .lm.fit(). 
#' The estimated effect size and stardard error for each variable are outputted.
#' @usage compute_summary_statistics(X, y, Z = NULL, center = TRUE,
#'                                   scale = TRUE, return_residuals = FALSE)
#'@param X n by p matrix of regressors.
#'@param y n-vector of response variables.
#'@param Z Optional n by k matrix of covariates to be included in all regresions. If Z is not NULL, the linear effects of covariates are removed from y first, and the resulting residuals are used in place of y.
#'@param center If center = TRUE, center X, y and Z.
#'@param scale If scale = TRUE, scale X, y and Z.
#'@param return_residuals Whether or not to output the residuals if Z is not NULL.
#'@details A list with two vectors containing the least-squares estimates of the coefficients (betahat) and their standard errors (sebetahat). Optionally, and only when a matrix of covariates Z is provided, a third vector residuals containing the residuals is returned. 
#'@examples
#'# Example
#'set.seed(1)
#'n = 400
#'p = 500
#'beta = rep(0,p)
#'beta[1] = 1
#'X = matrix(rnorm(n*p),nrow = n,ncol = p)
#'X = scale(X,center = TRUE,scale = TRUE)
#'y = drop(X %*% beta + rnorm(n))
#'SS=compute_summary_statistics(X,y)

compute_summary_statistics<-function (X, y, Z = NULL, center = TRUE, scale = TRUE, return_residuals = FALSE) {
  y_na = which(is.na(y))
  if (length(y_na)) {
    X = X[-y_na, ]
    y = y[-y_na]
  }
  if (center) {
    y = y - mean(y)
    X = scale(X, center = TRUE, scale = scale)
  }
  else X = scale(X, center = FALSE, scale = scale)
  X[is.nan(X)] = 0
  if (!is.null(Z)) {
    if (center) 
      Z = scale(Z, center = TRUE, scale = scale)
    y = .lm.fit(Z, y)$residuals
  }
    output = matrix(0, ncol(X), 2)
    for (i in 1:ncol(X)) {
      fit = summary(lm(y ~ X[, i]))$coef
      if (nrow(fit) == 2) 
        output[i, ] = as.vector(summary(lm(y ~ X[, i]))$coef[2,1:2])
      else output[i, ] = c(0, 0)
    }
  if (return_residuals && !is.null(Z)) 
    return(list(betahat = output[, 1], sebetahat = output[, 
                                                          2], residuals = y))
  else return(list(betahat = output[, 1], sebetahat = output[, 
                                                             2]))
}
