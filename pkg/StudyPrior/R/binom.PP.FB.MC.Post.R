# #' Title
# #'
# #' @param x
# #' @param n
# #' @param verbose
# #' @param length
# #' @param d.prior.a shape1 parameter for beta prior on weight
# #' @param d.prior.b shape2 parameter for beta prior on weight
# #' @param mc.cores
# #' @param samples
# #' @param focus
# #'
# #' @return
# #' @export
# #'
# #' @examples
# binom.PP.FB.MC.Post <- function(x, n, X, N, verbose=FALSE, length=30, d.prior.a=1, d.prior.b=1, mc.cores=1, samples=10000, focus){
#   n.hist <- length(x)
#
#   ddbinom <- function(x, size, prob, delta) dbinom(x,size,prob)^delta
#
#   cex <- function(d, x, n) integrate(function(p) sapply(p, function(PROB)
#     prod(mapply(ddbinom, x=x, size=n, delta=d, prob=PROB))),
#     lower=0,upper=1 )
#
#   p <- seq(0, 1, len=length)
#   if(!missing(focus)){
#     p2 <- unlist(sapply(focus, function(f) seq(from=f[1], to=f[2], length.out = f[3])))
#     p <- unique(sort(c(p,p2)))
#   }
#
#
#
#   # if(missing(dprior)) dprior <- function(d) (d>=0 & d<=1)*1
#
#   dens <- mclapply(p, function(PROB){
#     # print(PROB)
#
#     eval.f <- function(d) {
#       d <- c(d,1)
#       x <- c(x,X)
#       n <- c(n, N)
#       prod(mapply(ddbinom, x=x, size=n, delta=d, prob=PROB))/
#       cex(d,x,n)$value}
#
#     mean(apply(
#       matrix(rbeta(samples*n.hist, d.prior.a, d.prior.b),ncol=n.hist),
#       1,eval.f))
#   },mc.cores=mc.cores)
#
#   f <- splinefun(smooth.spline(x=p, y=dens))
#
#   k <- integrate(f, 0,1)$value
#
#   f <- splinefun(x=p, y=unlist(dens)/k)
#
#   g <- function(p,X) ifelse(0<=p&p<=1, f(p),0)
#
#   return(g)
#
# }
