#' Estimate the non linear mixed effect functional model.
#'
#' This function returns estimates of parameters in the
#' non linear functional mixed-effect model defined to warp
#' data and estimate the underlying model with mixed effect.
#'
#' Notice that the warping parameters are considered as random effects.
#'
#' @param t A vector of numbers, corresponding to time points.
#' @param y A matrix of numbers, corresponding to observations (size: T * n).
#' @param baseMu A B-spline used to decompose the global mean.
#' @param baseU A B-spline used to decompose the individual effects.
#' @param baseW A B-spline used to decompose the warping functions.
#' @param sigmaEpsilonTilde A number, defining the variance of the noise in the linear mixed-
#'        effect model fitted on the warping parameters.
#' @param threshold A number, defining the threshold of convergence.
#' @param nIte Maximum number of iterations
#'
#' @return A list, with fonct, functional quantities (indCurvAlign
#'         the individual aligned curves, warping the warping functions and theta,
#'         the warping parameters), para the estimates of parameters (alphaMu,
#'         sigmaU,theta0, sigmaTheta, sigmaEpsHat), dist the criterion computed to
#'         reach the convergence, and others other values (successAlphaMu,
#'         initTheta, initPara, CPUtime).
#' @examples
#' T = seq(0.5,0.841,length.out = 9)
#' n = 10
#' t = c(qnorm(T),1)
#' mu = cos(2*pi*t+pi/2)
#' library(fda)
#' baseMu = create.bspline.basis(c(0,max(t)), norder = 2, breaks = seq(0,1,0.5))
#' splineBasisMu = eval.basis(t,baseMu)
#' alphaMu = Data2fd(mu,argvals = t, baseMu)$coef
#' muApprox = (splineBasisMu) %*% alphaMu
#' baseU = create.bspline.basis(c(0,max(t)),norder = 2, breaks = seq(0,1,0.5))
#' mU = baseU$nbasis
#' sigmaU = diag(0.1,mU)
#' library(MASS)
#' alphaU = t(mvrnorm(n,rep(0,mU),sigmaU))
#' splineBasisU = eval.basis(t,baseU)
#' U = splineBasisU %*% alphaU
#' epsilon = t(mvrnorm(n,rep(0,length(t)),0.01 * diag(1,length(t))))
#' X = as.vector(muApprox) + U + epsilon
#' baseW = create.bspline.basis(c(0,max(t)), norder = 2, breaks = c(0,0.6,1))
#' mW = baseW$nbasis
#' splineBasisW = eval.basis(t,baseW)
#' theta = t(mvrnorm(n,rep(0,mW),diag(0.1,mW + 1e-3 * diag(1,mW))))
#' wtheta = matrix(rep(0,n*length(t)),ncol = n)
#' for (i in c(1:n)){
#'  wtheta[,i] = warpTimeFunction(splineBasisW,theta[,i],t)$warpTime
#' }
#' Y = matrix(0, nrow = length(t), ncol = n)
#' for (i in c(1:n)){
#'  y = approxfun(wtheta[,i],X[,i])
#'  Y[,i] = y(t)
#' }
#' warpMix(t,Y,baseMu, baseU, baseW, nIte = 2)
#'
#'
#'
#' @import fda
#' @import fields
#' @import MASS
#' @import reshape2
#' @import lme4
#' @importFrom stats approxfun as.formula cov fitted integrate optim update.formula
#' @export
#'
warpMix = function (t,y,baseMu,baseU,baseW,sigmaEpsilonTilde = 10^-3, threshold = 10^-3, nIte = 100){

  t0 = Sys.time()

  ##basis
  splineBasisMu = eval.basis(t,baseMu)
  splineBasisU = eval.basis(t,baseU)
  splineBasisW = eval.basis(t,baseW)

  ## dimension
  n = dim(y)[2]
  mW = dim(splineBasisW)[2]
  dimT = dim(y)[1]
  theta0 = rep(0,dim(splineBasisW)[2])
  ite=1;
  diffT = rep(0,length(t)-1)
  for (j in c(1:length(t)-1)){
    diffT[j] = t[j+1]-t[j]
  }

  ## initialization of theta (warping parameter)
  initTheta = initialisationTheta(t,y,splineBasisW)
  thetaIt = initTheta$theta0
  warptimeInit = initTheta$wT

  ## computation of theta0 and sigmaTheta
  theta0 = apply(thetaIt,1,mean)
  sigmaTheta = cov(t(thetaIt))

  ## initialization of alphaMu,sigmaEpsilon,alphaU,sigmaU
  initPara = initialisationPara(t,y,splineBasisMu,splineBasisU,warptimeInit)
  alphaMu = initPara$alphaMu
  sigmaU = initPara$sigmaU
  sigmaEpsilon = initPara$sigmaEpsilon
  indSignalChap = initPara$indSignal
  muChap = splineBasisMu %*% alphaMu

  d1=1
  d2 =1
  it = 1
  distGlob = 1
  diffCrit = 1
  Crit = 0
  crit = 1
  it = 0
  vectAlphaMu = matrix(0,ncol = 150,nrow = dim(splineBasisMu)[2])
  a = theta0
  warpTime = warptimeInit
  diff = rep(0,length(t)-1)
  somme = rep(0,n)
  cpt = 0
  while (((cpt <5)&&(it<nIte))){
    it = it+1

    ## observations for theta (warping parameter)
    est = estimationTheta(t,y, splineBasisW,indSignalChap, thetaIt)
    thetaObs = est$theta
    wThetaObs = est$wT

    ## global criterion
    for (i in c(1:n)){
      approxY = approxfun(wThetaObs[,i],y[,i])
      x = approxY(t)
      for (j in 1:(length(t)-1)){
        diff[j] = (x[j] - indSignalChap[j,i])^2 * (t[j+1]-t[j])
      }
      somme[i] = sqrt(sum(diff))
    }

    ## computation of theta by mixed effect modelling
    pred = predictionTheta(thetaObs,sigmaEpsilonTilde)
    thetaPred = pred$theta
    sigmaTheta = pred$sigmaE + sigmaEpsilonTilde * diag(1,dim(thetaPred)[1])
    for (i in c(1:n)){
      warpTime[,i] = warpTimeFunction(splineBasisW,thetaPred[,i],t)$warpTime
    }
    thetaIt = thetaPred
    theta0 = apply(thetaIt,1,mean)

    ## global criterion
    for (i in c(1:n)){
      approxY = approxfun(warpTime[,i],y[,i])
      x = approxY(t)
      for (j in 1:(length(t)-1)){
        diff[j] = (x[j] - indSignalChap[j,i])^2 * (t[j+1]-t[j])
      }
      somme[i] = sqrt(sum(diff))
    }

    ##for alphaMu,sigmaEpsilon,alphaU,sigmaU
    warpTime = wThetaObs
    para = majPara(t,y,splineBasisMu,splineBasisU,warpTime)
    alphaMu = para$alphaMu
    indSignalChap = para$indSignal
    vectAlphaMu[,it] = alphaMu

    ## global criterion
    for (i in c(1:n)){
      approxY = approxfun(warpTime[,i],y[,i])
      x = approxY(t)
      for (j in 1:(length(t)-1)){
        diff[j] = (x[j] - indSignalChap[j,i])^2 * (t[j+1]-t[j])
      }
      somme[i] = sqrt(sum(diff))
    }
    diffCrit = abs(crit - sum(somme))
    crit = sum(somme)
    if (diffCrit < threshold){
      cpt = cpt+1
    }
    else cpt = 0
    Crit = c(Crit,crit)

    ## Update parameters
    muChap = splineBasisMu %*% alphaMu
    sigmaU = para$sigmaU
    sigmaEpsHat = para$sigmaEpsilon
  }
  t1 = Sys.time()
  dT = t1-t0
  fonct = list(indCurvAlign = indSignalChap, warping = warpTime,theta = thetaIt)
  para = list(alphaMu = alphaMu, sigmaU = sigmaU, theta0 = theta0,
              sigmaTheta = sigmaTheta, sigmaEpsHat=sigmaEpsHat)
  dist = list(GlobalCrit = Crit)
  others = list(successAlphaMu = vectAlphaMu[, 1:it], initTheta = initTheta, initPara = initPara, CPUtime = dT)
  result = list(fonction = fonct,
                parameters = para,
                dist = dist,
                others = others)
  return(result)
}
