#' Empirical Bayes Power Prior for Binomial Data
#'
#' @param x historical events
#' @param n historical trials
#' @return A function of the probability parmater p
#' @export
#'
#' @examples
#'
binom.PP.EB.Beta <- function(x, n, X, N, verbose=FALSE){

  n.hist <- length(x)

  ddbinom <- function(x, size, prob, delta) dbinom(x,size,prob)^delta


  lik.d <- function(d) VGAM::dbetabinom.ab(X, N, 1+sum(d*x), 1+sum(d*(n-x)))

  elik <- function(d, choose.d, mpar, Mpar) {
    # this.d <- rep(1,n.hist)
    # this.d[choose.d] <- d[choose.d]
    # this.d *
      VGAM::dbetabinom.ab(X, N, 1+sum(d*x), 1+sum(d*(n-x))) *
      prod(dbeta(d,mpar*Mpar,Mpar-mpar*Mpar))
    }


  A <- 0.5
  B <- 2

  for(i in 1:10){

    # exp.d <- sapply(1:n.hist, function(D) {
    #   plot(D)
    #   cubature::adaptIntegrate(elik,
    #                            lowerLimit = rep(0,n.hist),
    #                            upperLimit = rep(1, n.hist),
    #                            choose.d = D,
    #                            m=A, M=B)})
    # new.d <- unlist(exp.d[1,])
    exp.d <-
      optimr::optimr(par=rep(0.5, n.hist),
                     fn = elik,
                     lower = rep(0,n.hist),
                     upper = rep(1, n.hist),
                     method = "L-BFGS-B",
                     control=list(maximize=TRUE, trace=1),
                     choose.d = 1:n.hist,
                     mpar=A, Mpar=B)

    new.d <- unlist(exp.d[1:n.hist])
    print(new.d)
    max.ab<- optimr::optimr(c(0.5,2),
                            function(par) prod(dbeta(new.d, par[1]*par[2],(1-par[1])*par[2])),
                            lower=c(1/sum(x),2),
                            upper=c(1-1/sum(x), sum(x)),
                            method = "L-BFGS-B",
                            control=list(maximize=TRUE))

    B <-  max.ab[[2]]

    if(max.ab[[1]] < 1/B) {
      A <- 1/B
    } else if(max.ab[[1]] > 1-1/B) {
      A <-1- 1/B
    } else A <- max.ab[[1]]


   print(c(A,B))
  }
    opd <- optimr::optimr(par = rep(0.6, n.hist),
                          fn = lik.d,
                          lower=rep(0, n.hist),
                          upper=rep(1, n.hist),
                          method = "L-BFGS-B",
                          control=list(maximize=TRUE))

  VP(opd)
  d <- unlist(opd[1:n.hist])

  f <- Vectorize(function(p) prod(mapply(ddbinom, x=x, size=n, delta=d, prob=p)))
  k <- integrate(f, 0,1)
  VP(k)
  g <- function(p) f(p)/k$value




return(g)

}
