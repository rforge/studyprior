
#' Calculate power
#'
#' @param prior
#' @param prob.range
#' @param length
#' @param n.binom
#'
#' @return
#' @export
#'
#' @examples
calc.power <- function(prior, prob.range=c(.5,.9), length=20, n.binom.control=30, n.binom.treatment=n.binom.control, treatment.difference=0.1, level=0.95, sig.mat, mc.cores=1){
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
