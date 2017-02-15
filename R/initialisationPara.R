#' Initialize the functional parameters (associated to the aligned curves).
#'
#' This function initializes the mean curve, the individual effect U_i, related to
#' aligned curves.
#'
#' @param t A vector of numbers, corresponding to time points.
#' @param y A matrix of numbers, corresponding to observations (size: T * n).
#' @param splineBasisMu A matrix, corresponding to the spline basis for
#'        the global mean function, evaluted in time points.
#' @param splineBasisU A matrix, corresponding to the spline basis for
#'        the individual curves, evaluted in time points.
#' @param warpTime A matrix, corresponding to warping time points.
#'
#' @return A list, with x, aligned curves, alphaMu the coefficients of the mean curve,
#'         sigmaEpsilon the variance of the noise, sigmaU the variance of the random effects,
#'         and indSignal each individual curves.
#'
#' @import lme4
initialisationPara = function(t,y,splineBasisMu,splineBasisU,warpTime){
  ## Initialization
  n = dim(y)[2]
  x = y

  ## warped process X
  for (i in c(1:n)){
    approxX = approxfun(warpTime[,i],y[,i])
    x[,i] = approxX(t)
  }

  ## melting the data to compute the initialisation of parameters in the mixed model
  splineBasisMuMelt=matrix(rep(0,dim(splineBasisMu)[1]*n * dim(splineBasisMu)[2]),
                           ncol = dim(splineBasisMu)[2])
  for (j in c(1:(dim(splineBasisMu)[2]))){
    splineBasisMuMelt[,j] = rep(splineBasisMu[,j],n)}

  splineBasisUMelt=matrix(rep(0,dim(splineBasisU)[1]*n * dim(splineBasisU)[2]),
                          ncol = dim(splineBasisU)[2])
  for (j in c(1: (dim(splineBasisU)[2]))){
    splineBasisUMelt[,j] = rep(splineBasisU[,j],n)}

  xMelt <- melt(x)
  XM = xMelt[,3]
  g = xMelt[,2]

  ## Initialization for the mixed model
  p = dim(splineBasisU)[2]
  f = as.formula(XM ~ -1 + splineBasisMuMelt)
  for (l in 1:p){
    f = update.formula(f,as.formula(paste( "~ . + (splineBasisUMelt[,",eval(l),"] | g)")))
  }

  ## Inference for the mixed effect model
  Res = lmer(f, REML = FALSE)
  ResS = summary(Res)
  alphaMu = ResS$coefficients[1:(dim(ResS$coefficients)[1])]
  sigmaU = as.data.frame(VarCorr(Res))[1:dim(splineBasisUMelt)[2],4]
  sigmaEpsilon = as.data.frame(VarCorr(Res))[dim(splineBasisUMelt)[2]+1,4]

  ## Fitted individual signal mu+U
  Ufitted = fitted(Res)
  uFittedMatrix = matrix(0, nrow = length(t),ncol = n)
  for (i in c(1:n)){
    uFittedMatrix[,i] = Ufitted[((i-1)*length(t)+1):(i*length(t))]
  }

  result = list(x = x, alphaMu = alphaMu, sigmaEpsilon = sigmaEpsilon, sigmaU = sigmaU, indSignal = uFittedMatrix)
  return(result)
}
