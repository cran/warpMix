#' Predict the warping parameters.
#'
#' This function predict the warping parameters, using the estimations of those parameters,
#' and fitting a linear mixed effect model on them.
#'
#' @param thetaObs A matrix (size: n * T) corresponding of the estimations of the warping parameters.
#' @param sigmaEpsilon A number, defining the variance of the noise in the linear mixed-
#'        effect model fitted on the warping parameters.
#'
#' @return A list, with theta, a matrix of predicted warping parameters,
#'         sigmaE the covariance of the random effects, and theta0 the mean.
#'
#'

predictionTheta = function(thetaObs,sigmaEpsilon){

  ## Initialization
  thetaObs = t(thetaObs)
  A = dim(thetaObs)
  n = A[1]
  T = A[2]

  ## Compute the prediction
  theta0hat = apply(thetaObs,2, mean)
  sigmaEhat = cov(thetaObs) - sigmaEpsilon * diag(1,T)
  effetAlea = sigmaEhat %*%  solve(cov(thetaObs)) %*% t(thetaObs)

  result = list(theta = effetAlea, sigmaE = sigmaEhat, theta0 = theta0hat)
  return(result)
}
