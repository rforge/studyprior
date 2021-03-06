#' Full Bayes Power Prior for Binomial Data
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param length  Number of points to evaluate density at
#' @param mc.cores Number of cores for parallel
#' @param verbose Print messages
#' @param p.prior.a shape1 parameter for beta prior on probability
#' @param p.prior.b shape2 parameter for beta prior on probability
#' @param dprior Density function for prior on weights
#'
#' @return A function of the probability parmater p
#' @export
#'
#'
binom.PP.FB <- function(x, n, verbose=FALSE, length=30, dprior, mc.cores=1,  p.prior.a=1, p.prior.b=1){
  n.hist <- length(x)

  ddbinom <- function(x, size, prob, delta) dbinom(x,size,prob)^delta

  cex <- function(d, x, n) integrate(function(p) sapply(p, function(PROB)
    prod(mapply(ddbinom, x=x, size=n, delta=d, prob=PROB))*dbeta(PROB, p.prior.a,p.prior.b)),
    lower=0,upper=1 )

  p <- seq(0, 1, len=length)

  if(missing(dprior)) dprior <- function(d) (d>=0 & d<=1)*1

  dens <-mclapply(p, function(PROB, verbose){
    VP(PROB)
    adaptIntegrate(function(d) prod(dprior(d))*
                     prod(mapply(ddbinom, x=x, size=n, delta=d, prob=PROB))*dbeta(PROB, p.prior.a,p.prior.b)/
                     cex(d,x,n)$value,
                   lowerLimit = rep(0, n.hist),
                   upperLimit = rep(1, n.hist),
                   maxEval = 5000)$integral
    },
    verbose=verbose,
    mc.cores=mc.cores)


  # dens <-mclapply(p, function(PROB, verbose){
  #   VP(PROB)
  #   R2Cuba::vegas(ndim=n.hist,
  #                 ncomp=1,
  #                 integrand =function(d) prod(dprior(d))*
  #                   prod(mapply(ddbinom, x=x, size=n, delta=d, prob=PROB))/
  #                   cex(d,x,n)$value,
  #                 lower=rep(0, n.hist),
  #                 upper=rep(1, n.hist))$value
  # },
  # verbose=TRUE,
  # mc.cores=mc.cores)

  f <- splinefun(x=p, y=dens)

  k <- integrate(f, 0,1)$value

  f <- splinefun(x=p, y=unlist(dens)/k)

  g <- function(p) ifelse(0<=p&p<=1, f(p),0)

  return(g)


}
