#' Compute the empirical L_2 distance related to the warping function.
#'
#' This function returns the empirical L_2 between two functions, the
#' first one being warped.
#'
#' @param t A vector of numbers, corresponding to time points.
#' @param f A vector of numbers, corresponding to the evaluated function.
#' @param g A vector of numbers, corresponding to the evaluated function.
#' @param theta A vector of warping parameters.
#' @param splineBasisW A matrix, corresponding to the spline basis for
#'        the warping functions, evaluted in time points.
#'
#' @return A list, with crit the distance.
#'
#'
#'
criterion = function (t,f,g,theta,splineBasisW){
  ## initialization
  T = length(t)
  diff = rep(0,T-1)

  ## warpTime
  A = warpTimeFunction(splineBasisW,theta,t)
  approxF = approxfun(A$warpTime,f)

  ## Computation of the criterion
  for (j in c(1:(T-1))) {
    diff[j] = (approxF(t[j]) - g[j])^2 * (t[j+1]-t[j])
  }
  crit = sqrt(sum(diff))

  result = list(crit = crit)
  return(result)
}
