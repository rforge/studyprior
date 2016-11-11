#' Normal power prior with full Bayes by Monte Carlo
#'
#' @param x numeric vector of results
#' @param sd standard deviations
#' @param verbose
#' @param length
#' @param dprior
#' @param mc.cores
#' @param samples
#' @param focus list(from, to, length) to evaluate the marginal density when fitting the curve
#'
#' @return
#' @export
#'
#' @examples
normal.PP.FB.MC <- function(x, sd, verbose=FALSE, length=30, d.prior.a=1, d.prior.b=1, mc.cores=1, samples=10000, focus){
  n.hist <- length(x)

  ddnorm <- function(x, mean, sd, delta) dnorm(x,mean,sd)^delta

  cex <- function(d, x, sd) integrate(function(p) sapply(p, function(MEAN)
    prod(mapply(ddnorm, x=x, sd=sd, delta=d, mean=MEAN))),
    lower=-Inf,upper=Inf )

  p <- seq(min(x)-max(sd)*6, max(x)+max(sd)*6, len=length)
  if(!missing(focus)){
    p2 <- unlist(sapply(focus, function(f) seq(from=f[1], to=f[2], length.out = f[3])))
    p <- unique(sort(c(p,p2)))
  }



  # if(missing(dprior)) dprior <- function(d) (d>=0 & d<=1)*1

  dens <- mclapply(p, function(MEAN){
    if(verbose) print(MEAN)

    eval.f <- function(d){ #prod(dprior(d))*
      prod(mapply(ddnorm, x=x, sd=sd, delta=d, mean=MEAN))/
      cex(d,x,sd)$value}

    mean(apply(
      matrix(rbeta(samples*n.hist,d.prior.a,d.prior.b),ncol=n.hist),
      1,eval.f))
  },mc.cores=mc.cores)


  f <-  splinefun(smooth.spline(x=p, y=unlist(dens)))

  g <- function(p) {
    z <- f(p)
    ifelse(0 > z, 0, z)
  }

  k <- integrate(g, min(p),max(p))$value

  z <- smooth.spline(x=p, y=unlist(dens)/k)

  f <- splinefun(z$x, ifelse(z$y<0, 0, z$y))

  g <- function(p) {
    z <- f(p)
    ifelse(0 > z, 0, z)
  }

  return(g)

}
