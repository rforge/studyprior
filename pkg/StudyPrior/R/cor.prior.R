#'
#' #' Correlation prior
#' #'
#' #' @param delta
#' #' @param rho
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' cor.prior <- function(delta, rho, rooti){
#'
#'
#'   k <- length(delta)
#'
#'   y <- qnorm(delta)
#'
#'   det <- sqrt((2*pi)^k) * exp(sum(y^2) / 2)
#'
#'   # f <- function(i) -y[i]*(-y[i] -(k-1)*y[i]*rho + rho*sum(x[-i])) / (rho-1)*(1+rho*(k-1))
#'   #
#'   # det * exp(0.5*sum(sapply(seq(k), f )))/
#'   #    sqrt( (-1)^(k-1) * (1 + (k-1) * rho) * (rho -  1)^(k-1)  )
#'   sigma <- matrix(rho, ncol=k, nrow=k)
#'   diag(sigma) <- 1
#'
#'   if(missing(rooti)) rooti <- backsolve(chol(sigma),diag(k))
#'
#'   quads <- colSums((crossprod(rooti,y))^2)
#'   exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads) * det
#'
#' }
#'
#' # Code from bayesm
#' # seen at gallery.rcpp.org/articles/dmvnorm_arma
#' #
#' # dMvn <- function(X,mu,Sigma)
#' # k <- ncol(X)
#' # rooti <- backsolve(chol(Sigma),diag(k))
#' # quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
#' # return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))
