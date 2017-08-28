

#' Calculate a posterior density based on a binomial likelihood
#'
#' @param prior Prior
#' @param X Number of successes 
#' @param N Number of patients 
#'
#' @return The posterior density function
#' @export
#'
#'
#'
calc.posterior <- function(prior, X, N){
  pr <- function(p) prior(p,X)
  k <- adaptIntegrate(function(p) pr(p)*dbinom(X,N,p),
                      lowerLimit = 0,
                      upperLimit = 1)$integral
  function(p) pr(p)*dbinom(X,N,p)/k
}
