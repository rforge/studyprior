#' Empirical Bayes Power Prior for Binomial Data Esimated Separately
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param X number of new successes
#' @param N number of new patients
#' @param verbose Print messages
#'
#' @return A function of the probability parmater p
#' @export
#'
#'
binom.PP.EB.Sep <- function(x, n, X, N, verbose=FALSE){

  n.hist <- length(x)

  ddbinom <- function(x, size, prob, delta) dbinom(x,size,prob)^delta
  dists <-
    lapply(X, function(X){
      lik.d <- function(d, xh, nh) VGAM::dbetabinom.ab(X, N, 1+(d*xh), 1+(d*(nh-xh)))
      opd <-sapply(1:n.hist, function(i){
        optimize(f = lik.d,
                 lower=0,
                 upper=1,
                 maximum = TRUE,
                 xh=x[i],
                 nh=n[i])$maximum
      })

      VP(opd)
      d <- opd

      f <- Vectorize(function(p) prod(mapply(ddbinom, x=x, size=n, delta=d, prob=p)))
      k <- integrate(f, 0,1)
      VP(k)
      g <- function(p) f(p)/k$value

      return(g)
})
  ds <- lapply(dists, function(D) get(x="d", pos=environment(D)))
  f <- function(p,X) do.call(dists[[X+1]], list(p=p))
  return(f)
}
