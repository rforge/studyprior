#' Full Bayes Power Prior for Binomial Data
#'
#' @param x historical events
#' @param n historical trials
#' @param tau.prior optional prior on heterogeneity parameter
#'
#' @return A function of the probability parmater p
#' @export
#'
#' @examples
#'
binom.PP.FB <- function(x, n, verbose=FALSE, length=30, dprior, mc.cores=1,  p.prior.a=1, p.prior.b=1){
  n.hist <- length(x)

  ddbinom <- function(x, size, prob, delta) dbinom(x,size,prob)^delta

  cex <- function(d, x, n) integrate(function(p) sapply(p, function(PROB)
    prod(mapply(ddbinom, x=x, size=n, delta=d, prob=PROB))*dbeta(PROB, p.prior.a,p.prior.b)),
    lower=0,upper=1 )

  p <- seq(0, 1, len=length)

  if(missing(dprior)) dprior <- function(d) (d>=0 & d<=1)*1

  dens <-parallel::mclapply(p, function(PROB, verbose){
    VP(PROB)
    cubature::adaptIntegrate(function(d) prod(dprior(d))*
                     prod(mapply(ddbinom, x=x, size=n, delta=d, prob=PROB))*dbeta(PROB, p.prior.a,p.prior.b)/
                     cex(d,x,n)$value,
                   lower=rep(0, n.hist),
                   upper=rep(1, n.hist),
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
