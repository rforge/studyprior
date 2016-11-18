#' Empirical Bayes Power Prior for Binomial Data
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
binom.PP.EB <- function(x, n, X, N, verbose=FALSE, mc.cores=1){
  #if X isn't specified we calculate it for all, set flag too
  if(missing(X)) {
    X <- 0:N
    X.only <- FALSE
  } else x.only=TRUE

 dists <-
   mclapply(mc.cores=mc.cores,
            X,
            function(X){

    n.hist <- length(x)

    ddbinom <- function(x, size, prob, delta) dbinom(x,size,prob)^delta


    lik.d <- function(d) VGAM::dbetabinom.ab(X, N, 1+sum(d*x), 1+sum(d*(n-x)))

    # opd <- optimr::optimr(par = rep(.005, n.hist),
    #                       fn = lik.d,
    #                       lower=rep(0, n.hist),
    #                       upper=rep(1, n.hist),
    #                       method = "L-BFGS-B",
    #                       control=list(maximize=TRUE,
    #                                    fnscale=1.0e-20))
    opd <- BB::spg(par = rep(.005, n.hist),
               fn = lik.d,
               lower=rep(0, n.hist),
               upper=rep(1, n.hist),
               control=list(maximize=TRUE,
                            trace=FALSE))


    if(opd$convergence!=0) print(opd)

    d <- opd$par
    print(d)

    f <- Vectorize(function(p) prod(mapply(ddbinom, x=x, size=n, delta=d, prob=p)))
    k <- integrate(f, 0,1)
    VP(k)
    g <- function(p) f(p)/k$value

    return(g)
  })
  f <- function(p,X) do.call(dists[[X+1]], list(p=p))
 return(f)
}
