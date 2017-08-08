#' Empirical Bayes Power Prior for Binomial Data
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param X vector of new successes to computer prior for
#' @param N number of new patients
#' @param verbose TRUE
#' @param mc.cores number of cores for parallel
#' @param p.prior.a shape1 parameter of initial beta prior for successes
#' @param p.prior.b shape2 parameter of initial beta prior for successes 
#'
#' @return A function of the probability parmater p
#' @export
#'
#'
binom.PP.EB <- function(x, n, X, N, verbose=FALSE, mc.cores=1, p.prior.a=1, p.prior.b=1){
  #if X isn't specified we calculate it for all, set flag too
  if(missing(X)) {
    X <- 0:N
    X.only <- FALSE
  } else X.only=TRUE

 ds <-
   mclapply(mc.cores=mc.cores,
            X,
            function(X){


    lik.d <- function(d) VGAM::dbetabinom.ab(X, N, p.prior.a+sum(d*x), p.prior.b+sum(d*(n-x)))

    # opd <- optimr::optimr(par = rep(.005, n.hist),
    #                       fn = lik.d,
    #                       lower=rep(0, n.hist),
    #                       upper=rep(1, n.hist),
    #                       method = "L-BFGS-B",
    #                       control=list(maximize=TRUE,
    #                                    fnscale=1.0e-20))
    opd <- BB::spg(par = rep(.005, length(x)),
               fn = lik.d,
               lower=rep(0, length(x)),
               upper=rep(1, length(x)),
               control=list(maximize=TRUE,
                            trace=FALSE))


    if(opd$convergence!=0) print(opd)

    d <- opd$par



    return(d)
  })

  f <- function(p,X) {
    d <- ds[[X+1]]
    dbeta(p,p.prior.a + sum(x*d), p.prior.b + sum(d*(n-x)))
  }
 return(f)
}
