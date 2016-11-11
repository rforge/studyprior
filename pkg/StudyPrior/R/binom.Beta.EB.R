#' Empirical Bayes Beta Prior for Binomial Data
#'
#' @param x historical events
#' @param n historical trials
#'
#' @return A function of the probability parmater p
#' @export
#'
#' @examples
#'
binom.Beta.EB <- function(x, n, verbose=FALSE){
  n.hist <- length(x)

  opMmu <- function(par) prod(
    mapply(VGAM::dbetabinom.ab,
           x=x,
           size=n,
           shape1=par[1]*par[2],
           shape2=par[1]*(1-par[2])))

  op <- optimr::optimr(c(sum(x), mean(x/n)),
                       opMmu,
                       lower=c(0,0),
                       upper=c(sum(n), 1),
                       method="L-BFGS-B",
                       control=list(fnscale=-1*0.1^n.hist))

  VP(op)

  f <- function(p,...) dbeta(p, op[[1]]*op[[2]],op[[1]]*(1-op[[2]]))

  return(f)

}
