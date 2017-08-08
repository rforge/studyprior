
#' Calculate power
#'
#' @param prior Prior to calculate posterior
#' @param prob.range Range of values to calculate over
#' @param length Number of values to calculate for
#' @param n.binom.control Number of patients in new trial's contral arm
#' @param n.binom.treatment Number of patients in new trial's treatment arm
#' @param treatment.difference Predefined treatment difference
#' @param level Significance level for calculation of sig.mat
#' @param sig.mat Precalculated significance matrix (See sig.matrix)
#' @param mc.cores Number of cores for parallel
#'
#' @return Vector or power values
#' @export
#'

calc.power <- function(prior, prob.range=c(.5,.9), length=20, n.binom.control=30,
                       n.binom.treatment=n.binom.control, treatment.difference=0.1, level=0.95, sig.mat, mc.cores=1){
# the probability to detect a given true difference
# detection by
  if(treatment.difference+prob.range[2] > 1) stop("Unable to calculate power for difference at upper end of probability range (prob.range+treatment.difference > 1).")

  if(missing(sig.mat)) sig.mat <- sig.matrix(n.binom.control,
                                             n.binom.treatment, level, prior,
                                             treat.beta.prior.par=c(1,1), mc.cores=mc.cores)


  probs <- sapply(seq(prob.range[1],prob.range[2], len=length),
                  function(P){
    matrix(dbinom(0:n.binom.control, n.binom.control, prob=P),nrow=1) %*%
      sig.mat %*%
      matrix(dbinom(0:n.binom.treatment, n.binom.treatment, prob=P+treatment.difference),ncol=1)
  })

  names(probs) <- seq(prob.range[1],prob.range[2], len=length)

  return(probs)
}
