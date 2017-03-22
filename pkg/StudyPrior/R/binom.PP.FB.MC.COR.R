#' Title
#'
#' @param x
#' @param n
#' @param verbose
#' @param length
#' @param d.prior.a shape1 parameter for beta prior on weight
#' @param d.prior.b shape2 parameter for beta prior on weight
#' @param mc.cores
#' @param samples
#' @param p.prior.a shape1 parameter for beta prior on probability
#' @param p.prior.b shape2 parameter for beta prior on probability
#' @param focus
#'
#' @return
#' @export
#'
#' @importFrom mvtnorm rmvnorm
#'
#' @examples
binom.PP.FB.MC.COR <- function(x, n, verbose=FALSE, length=30, d.prior.cor=0, p.prior.a=1, p.prior.b=1, mc.cores=1, samples=10000, focus ){
  n.hist <- length(x)

  #Where to calculate density at
  p <- seq(0.0001, .9999, len=length)
  if(!missing(focus)){
    p2 <- unlist(sapply(focus, function(f) seq(from=f[1], to=f[2], length.out = f[3])))
    p <- unique(sort(c(p,p2)))
  }



  dens <- mclapply(p, function(PROB){
    # print(PROB)

    eval.f <- function(d) dbeta(PROB, p.prior.a+sum(d*x), p.prior.b+sum(d*(n-x)))

    #correlation matrix
    sigma <- matrix(d.prior.cor, ncol=n.hist, nrow = n.hist)
    diag(sigma) <- rep(1, n.hist)

    mean(apply(
      pnorm(rmvnorm(samples, mean=rep(0,n.hist), sigma=sigma)),
      1,eval.f))
  },mc.cores=mc.cores)

  dens2 <- density(rep(p,unlist(dens)*10000),from=0,to=1, bw=.02)
  # f <- splinefun(smooth.spline(x=p, y=pmax(dens,0)))


  # k <- integrate(f, 0,1)$value

  # f <- splinefun(x=p, y=unlist(dens)/k)

  # g <- function(p,X) ifelse(0<=p&p<=1, f(p),0)

  g <- approxfun(dens2, yleft=0, yright=0)
  formals(g) <- alist(v=,...=)

  return(g)

}
