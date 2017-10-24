#' Full Bayes Power Prior with Correlated Weight Parameters using Rao-Blackwell
#'
#' @param x number of historical successes
#' @param n number historical patients
#' @param verbose Print messages
#' @param mixture.size Number of components in the mixture. If d.prior.cor=0, will be rounded to nearest value a^length(x) to make the grid evenly spaced.
#' @param mc.cores Number of cores for parallel
#' @param p.prior.a shape1 parameter for initial beta prior on probability
#' @param p.prior.b shape2 parameter for initial beta prior on probability
#' @param d.prior.cor Correlation parameter for transformed multivariate normal prior on weights
#'
#' @return A density function
#' @export
#'
#' @importFrom mvtnorm rmvnorm
#'
binom.PP.FB.COR <- function(x, n, verbose=FALSE, mixture.size=1000, d.prior.cor=0, p.prior.a=1, p.prior.b=1, mc.cores=1){
  n.hist <- length(x)
  
  
  if(d.prior.cor==1){#pooled case
    sumx <- sum(x)
    sumnx <- sum(n)-sumx
    D <- seq(0,1,len=mixture.size)
    
    pars <- sapply(D, function(d) c(d*sumx+1, d*sumnx+1) )
    
    mix <- create.mixture.prior("beta", pars, weights=rep(1/mixture.size,mixture.size))
    
    
  } else if(d.prior.cor==0){#independent case
    leng <- round(mixture.size^(1/n.hist))
    mixture.size <- leng^n.hist
    
    D <- expand.grid(rep(list(seq(0,1,len=leng)),n.hist))
    pars <- apply(D,1, function(d) c(d%*%x+1, d%*%(n-x)+1) )
    
    mix <- create.mixture.prior("beta", pars, weights=rep(1/mixture.size,mixture.size))
    
    
  } else{ #anything else
    sigma <- matrix(d.prior.cor, ncol=n.hist, nrow = n.hist)
    diag(sigma) <- rep(1, n.hist)
    
    D <- pnorm(rmvnorm(mixture.size, mean=rep(0,n.hist), sigma=sigma) )
    
    pars <- apply(D,1, function(d) c(d%*%x+1, d%*%(n-x)+1) )
    
    mix <- create.mixture.prior("beta", pars, weights=rep(1/mixture.size,mixture.size))
  }
  
  
  f <- function(p,X) eval.mixture.prior(p, mix)
  
  return(f)
  
}
