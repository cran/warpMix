#' Compute the warped time points
#'
#' This function returns warped time points for a known warping parameter theta.
#'
#' @param splineBasisW A matrix, corresponding to the spline basis for
#'        the warping functions, evaluted in time points.
#' @param theta A matrix, corresponding to initial values for the warping parameters.
#' @param t A vector of numbers, corresponding to time points.
#'
#' @export
#' @return A vector, corresponding to the warped time points.
#'

warpTimeFunction = function(splineBasisW,theta,t){
  ## Tests
  if (is.vector(t) == FALSE || is.matrix(splineBasisW) == FALSE || is.vector(theta) == FALSE)
    stop(paste(sQuote("t"), "must be a vector,"," the other parameters must be matrices"))

  if ((length(t) == dim(splineBasisW)[1]) == FALSE)
    stop(paste(sQuote("splineBasisW"), "must have as much rows as the length of ", sQuote("t")))

  ## Initialization
  T = length(t)
  mW = dim(splineBasisW)[2]

  ## Functions h tilde and h
  hTilde = splineBasisW %*% theta
  hTildeFun = approxfun(t,hTilde)
  h = hTilde - integrate(hTildeFun,0,1,stop.on.error = FALSE)$value

  ## Function w
  hFun = approxfun(t,exp(h))
  num1 = rep(0,length(t))
  num = num1

  for (j in c(1:(length(t)-1))){
    num1[j+1] = integrate (hFun, t[j], t[j+1])$value
  }
  num = cumsum(num1)
  w = num/num[length(t)]

  result = list(warpTime = w)
  return(result)
}
