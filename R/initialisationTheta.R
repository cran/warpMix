#' Initialize the warping parameters.
#'
#' This function initializes the warping parameters
#'
#' @param t A vector of numbers, corresponding to time points.
#' @param y A matrix of numbers, corresponding to observations (size: T * n).
#' @param splineBasisW A matrix, corresponding to the spline basis for
#'        the warping functions, evaluted in time points.
#'
#' @return A list, with theta, a matrix of estimated warping parameters,
#'         and wT, the corresponding warping functions.
#'
#'
#' @import fda
initialisationTheta = function(t,y,splineBasisW) {

  ## Initialization
  n = dim(y)[2]
  mW = dim(splineBasisW)[2]
  thetaInit = rep(0,mW)
  thetaOpt = matrix(rep(0,n*mW),ncol=n)

  ## Initialization of the central curve
   functBoxPlot = fbplot(y,t,plot=FALSE)
   MedCurv = y[,functBoxPlot$medcurv]
   warpTime = matrix(rep(0, n* length(t)),ncol = n)

  ## Computing the warping function -> computing the warping parameters theta
  for (i in c(1:n)) {
    crit = function (theta) criterion(t,y[,i],MedCurv,theta,splineBasisW)
    a = optim(thetaInit,crit)
    thetaOpt[,i]  = a$par
    warpTime[,i] = warpTimeFunction(splineBasisW,thetaOpt[,i] ,t)$warpTime
  }

  result = list(theta0 = thetaOpt, wT = warpTime)
  return(result)
}
