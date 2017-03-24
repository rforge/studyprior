#' Full Bayesian Meta-Analytic Prior for Binomial Data
#'
#' @param x historical events
#' @param n historical trials
#' @param tau.prior optional prior on heterogeneity parameter
#' @param verbose
#'
#' @return A function of the probability parmater p
#' @export
#'
#'
binom.MAP.FB <- function(x, n, tau.prior, verbose=FALSE){
  n.hist <- length(x)
  dat <- data.frame(x=c(x,NA), n=c(n,NA), z=1:(n.hist+1))

  if(missing(tau.prior)){
    prior <-  list(prior= "logtnormal", param=c(0,1))
  } else {
    prior <- tau.prior
  }

  formula <- x ~ 1 + f(z, model="iid",
                       hyper = list(theta = prior))

  result <- INLA::inla(formula,
                       data = dat,
                       family = "binomial",
                       control.fixed = list(mean.intercept = 0, prec.intercept = 1/1000),
                       Ntrials=n,
                       control.predictor = list(compute=TRUE, link=1))


  # marg.string <- paste0(capture.output(dput(result$marginals.fitted.values[n.hist+1])),
                        # collapse = "")


#
#   fun.string <- paste0(
#     "function(p) {
#     marg <-",marg.string,"
#     inla.dmarginal(p, marg)
#   }", collapse="")
#
# f <- eval(parse(text=fun.string))



 f <-  splinefun(result$marginals.fitted.values[[n.hist+1]][,1],
                 result$marginals.fitted.values[[n.hist+1]][,2])

  rm(dat, formula, n, n.hist, prior, result, tau.prior, verbose, x)


   function(p,x) {
    dens <- rep(0,length(p))
    i <- which(0<p&p<1)
    dens[i] <- f(p[i])
    dens
  }

}
