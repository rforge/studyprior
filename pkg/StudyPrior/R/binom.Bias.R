#' Full Bayes Bias Prior for Binomial Data
#'
#' @param x historical events
#' @param n historical trials
#' @param tau.prior optional prior on bias parameter
#'
#' @return A function of the probability parmater p
#' @export
#'
#' @examples
#'
binom.Bias.FB <- function(x, n, verbose=FALSE, length=30, dprior, mc.cores=1){
  n.hist <- length(x)

  dat <- data.frame(x, n, z=1:n.hist)

  # if(missing(tau.prior)){
    prior <-  list(prior= "logtnormal", param=c(0,1))
  # } else {
  #   prior <- tau.prior
  # }

  formula <- x ~ 1 + f(z, model="iid",
                       hyper = list(prec = prior))

  result <- INLA::inla(formula,
                       data = dat,
                       family = "binomial",
                       control.fixed = list(mean.intercept = 0, prec.intercept = 1/1000),
                       Ntrials=n,
                       control.predictor = list(compute=TRUE, link=1))

k <- integrate(Vectorize(function(p) prod(sapply(result$marginals.fitted.values, inla.dmarginal, x=p))),
               lower=0,
               upper=1)$value


  f <- function(p) sapply(p, function(p) prod(sapply(result$marginals.fitted.values, inla.dmarginal, x=p))/k)



}
