#' Full Bayesian Meta-Analytic Prior for Binomial Data
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

  f <- inla.dmarginal

  A <- result$marginals.fitted.values[[n.hist+1]][,1]
  B <- result$marginals.fitted.values[[n.hist+1]][,2]

  formals(f) <- list(x=NULL,
                     marginal=list(x=A,y=B),
                     log=FALSE)

  g <- function(p,X) f(p)
  return(g)
}
