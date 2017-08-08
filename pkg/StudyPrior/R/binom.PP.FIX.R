#' Fixed Weight Power Prior for Binomial Data
#'
#' @param x historical events
#' @param n historical trials
#' @param d fixed weight
#' @param p.prior.a shape1 parameter for beta prior on p
#' @param p.prior.b shape2 parameter for beta prior on p
#' @param verbose Print messages
#'
#' @return Probability density function for paramater p
#' @export
#'
#'
binom.PP.FIX <- function(x, n, d,  p.prior.a=1, p.prior.b=1, verbose){
  ds <- d
  
    function(p,...) dbeta(p,p.prior.a + sum(x*d), p.prior.b + sum(d*(n-x)))
}
