#' Estimate the warping parameters.
#'
#' This function estimate the warping parameters, knowing the observations and
#' the individual aligned curves.
#'
#' @param t A vector of numbers, corresponding to time points.
#' @param y A matrix of numbers, corresponding to observations (size: T * n).
#' @param splineBasisW A matrix, corresponding to the spline basis for
#'        the warping functions, evaluted in time points.
#' @param indSignal A matrix, corresponding to the individual aligned curves.
#' @param thetaObs A matrix, corresponding to initial values for the warping parameters.
#'
#' @return A list, with theta, a matrix of estimated warping parameters,
#'         and wT, the corresponding warping functions.
#'
#'
estimationTheta = function(t,y, splineBasisW,indSignal, thetaObs) {

  ## Initialization
  n = dim(y)[2]
  mW = dim(splineBasisW)[2]
  thetaInit = thetaObs
  thetaOpt = matrix(rep(0,n*mW),ncol=n)
  warpTime = matrix(rep(0, n* length(t)),ncol = n)

  ## Computing the warping function -> computing the warping parameters theta
  for (i in c(1:n)) {
    MedCurv = indSignal[,i]
    crit = function (theta) criterion(t,y[,i],MedCurv,theta,splineBasisW)
    a = optim(thetaInit[,i],crit)
    thetaOpt[,i]  = a$par
    warpTime[,i] = warpTimeFunction(splineBasisW,thetaOpt[,i] ,t)$warpTime
  }

  result = list(theta = thetaOpt, wT = warpTime)
  return(result)
}
