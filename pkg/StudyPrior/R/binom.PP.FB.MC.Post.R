#' #'  Posterior of Fully Bayes Power Prior
#' #'
#' #' @param x number of historical successes
#' #' @param n number historical patients
#' #' @param length  Number of points to evaluate density at
#' #' @param mc.cores Number of cores for parallel
#' #' @param samples Number of Monte Carlo samples
#' #' @param focus List of triples specifying regions to focus on, eg peaks. Specified as list(c(lower, upper, length))
#' #' @param verbose Print messages
#' #' @param d.prior.a shape1 parameter for beta prior on weight
#' #' @param d.prior.b shape2 parameter for beta prior on weight
#' #' @param p.prior.a shape1 parameter for beta prior on probability
#' #' @param p.prior.b shape2 parameter for beta prior on probability
#' #' @param X Number of new successes
#' #' @param N Number of new patients
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' binom.PP.FB.MC.Post <- function(x, n, X, N, verbose=FALSE, length=30, d.prior.a=1, d.prior.b=1,p.prior.a=1, p.prior.b=1, mc.cores=1, samples=10000, focus){
#'   n.hist <- length(x)
#' 
#'   #Where to calculate density at
#'   p <- seq(0.0001, .9999, len=length)
#'   if(!missing(focus)){
#'     p2 <- unlist(sapply(focus, function(f) seq(from=f[1], to=f[2], length.out = f[3])))
#'     p <- unique(sort(c(p,p2)))
#'   }
#' 
#' 
#'   # if(missing(dprior)) dprior <- function(d) (d>=0 & d<=1)*1
#' 
#'   dens <- mclapply(p, function(PROB){
#'     # print(PROB)
#' 
#'     eval.f <- function(d) {
#'       # d <- c(d,1)
#'       # x <- c(x,X)
#'       # n <- c(n, N)
#'       # prod(mapply(ddbinom, x=x, size=n, delta=d, prob=PROB))/
#'       # cex(d,x,n)$value
#'       dbeta(PROB, p.prior.a+sum(d*x), p.prior.b+sum(d*(n-x)))* PROB^X*(1-PROB)^(N-X)
#'       }
#' 
#' 
#'     mean(apply(
#'       matrix(rbeta(samples*n.hist, d.prior.a, d.prior.b),ncol=n.hist),
#'       1,eval.f))
#'   },mc.cores=mc.cores)
#' 
#'   f <- splinefun(smooth.spline(x=p, y=dens))
#' 
#'   k <- integrate(f, 0,1)$value
#' 
#'   f <- splinefun(x=p, y=unlist(dens)/k)
#' 
#'   g <- function(p,X) ifelse(0<=p&p<=1, f(p),0)
#' 
#'   return(g)
#' 
#' }
