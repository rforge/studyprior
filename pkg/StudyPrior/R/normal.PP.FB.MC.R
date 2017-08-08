#' Normal power prior with full Bayes by Monte Carlo
#'
#' @param x numeric vector of results
#' @param sd standard deviations
#' @param verbose Print messages
#' @param length Number of points to evaluate density at
#' @param mc.cores Number of cores for parallel
#' @param samples Number of Monte Carlo samples
#' @param focus List of triples specifying regions to focus on, eg peaks. Specified as list(c(lower, upper, length))
#' @param d.prior.a shape1 parameter for beta prior on weights
#' @param d.prior.b shape2 parameter for beta prior on weights
#'
#' @return Density function 
#' @export
#'
normal.PP.FB.MC <- function(x, sd, verbose=FALSE, length=30, d.prior.a=1, d.prior.b=1, mc.cores=1, samples=10000, focus){
  n.hist <- length(x)

  ddnorm <- function(x, mean, sd, delta) dnorm(x,mean,sd)^delta

  # cex <- function(d, x, sd) integrate(function(p) sapply(p, function(MEAN)
  #   prod(mapply(ddnorm, x=x, sd=sd, delta=d, mean=MEAN))),
  #   lower=-Inf,upper=Inf )



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
      cex(d,x,sd)}

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


cex <- function(d,x,sd){
  n.hist <- length(x)
  h <- seq(n.hist)
  r <- 1/sd^2
  exp.part <- 1/(-2*sum(d*r)) * sum(sapply(h, function(i) sum(d[i]*d[h[-i]]*r[i]*r[h[-i]]*(x[i]-x[h[-i]])^2)))
  prod.part <- prod(sqrt(r/2/pi)^d)

a <-   prod.part * sqrt(2*pi / sum(d*r)) * exp(exp.part)
print(a)
return(a)
}

#
# cex2 <- function(d, x, sd) integrate(function(p) sapply(p, function(MEAN)
#   prod(mapply(ddnorm, x=x, sd=sd, delta=d, mean=MEAN))),
#   lower=-Inf,upper=Inf )
#
# cex(c(.4,1),c(1,1.4),c(.3,.1))
# cex2(c(.4,1),c(1,1.4),c(.3,.1))
