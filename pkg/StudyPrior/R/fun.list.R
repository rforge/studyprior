#' Create mixture model
#'
#' @param type Type of mixture, "normal" or "beta"
#' @param pars Parameter values
#' @param weights Mixture weights
#'
#' @return A \code{mixture.prior} object
#' @export
#'
create.mixture.prior <- function(type, pars, weights){
  degree <- length(weights)

  fl <- rep(list(switch(type, "normal" = dnorm,"beta" = dbeta)),
            degree)

  attr(fl,"type") <- type
  attr(fl,"par.type") <- switch(type,
                                "normal" = c("mean","sd"),
                                "beta" = c("shape1","shape2") )

  fl <- update.mixture.prior(fl, pars=pars, weights=weights)
  class(fl) <- "mixture.prior"
  return(fl)
}

#' Test if object is a mixture.prior
#'
#' @param x Object 
#'
#' @return TRUE if \code{x} is a \code{mixture.prior}, else FALSE
#' @export
#'
is.mixture.prior <- function(x) {
  inherits(x, "mixture.prior")
}

#' Update parameters or weights in mixture model
#'
#' @param ... Unused
#' @param object mixture.prior to update
#' @param pars New paramers 
#' @param weights New weights
#'
#' @return A \code{mixture.prior} object with updated weights and parameters
#' @export
#' @method update mixture.prior
update.mixture.prior <- function(object, ..., pars, weights ){
  
  # args <- eval(substitute(alist(...)))
  # pars <- args$pars
  # weights <- args$weights
  
  if(!missing(pars)){
    par.names <- attr(object, "par.type")
    split.pars <- split(unlist(pars), rep(1:length(object), each=length(par.names)))
    attr(object,"pars") <- unlist(pars)
    for(i in 1:length(object)){
      tmp <- formals(object[[i]])
      tmp[par.names] <- split.pars[[i]]
      formals(object[[i]]) <- tmp
    }
  }
  if(!missing(weights)){
    attr(object,"weights") <- weights/sum(weights)
  }

  return(object)
}



#' Get density of mixture
#'
#' @param x Vector of values to evaluate density at
#' @param mixture.prior Mixture prior
#' @param weights Optional weights different to those in mixture.prior
#' @param subset Vector of values specifying which mixture components should be included, eg c(1,2)
#'
#' @return Densities at of prior at x 
#' @export 
#'
eval.mixture.prior <- function(x, mixture.prior, weights, subset){
  if(missing(subset)) subset <- 1:length(mixture.prior)
  X <- rep(0,length(x))

  if(missing(weights)) weights <- attr(mixture.prior, "weights")
  w <- weights/sum(weights)

  for(i in subset){
    X <- X + w[i]*do.call(mixture.prior[[i]], list(x=x))
  }
  X
}

#' Plot mixture model
#'
#' @param x Vector of values to evaluate density at
#' @param mixture.prior Mixture prior
#' @param stack Plot the mixture components as stacked regions
#' @param lines.only Draw lines only (to be used over an existing plot)
#' @param ... Additional arguments to \code{\link{lines}}
#'
#' @export
#'

plot.mixture.prior <- function(mixture.prior, x, stack=FALSE, lines.only=FALSE, ...){

  if(missing(x)) x <- seq(0,1, length=200)
      
  Y <- eval.mixture.prior(x, mixture.prior)
  if(!lines.only) plot(x,Y, type='n', ...)

  if(stack==TRUE){
    degree <- length(mixture.prior)
    w <- attr(mixture.prior, "weights")
    z <- order(w, decreasing = TRUE)

    for(i in 1:(degree-1)){

        lines(x,eval.mixture.prior(x, mixture.prior,subset=z[1:i]), type='l', ...)
    }
  }
  lines(x,Y, type='l', ...)
}

#' Calculate the equivalent sample size of mixture model
#' 
#' @param mixture.prior Mixture prior object
#'
#' @return Sample size
#' @export
#'
ess.mixture.prior <- function(mixture.prior){
  p <- attr(mixture.prior,"pars")
  w <- attr(mixture.prior,"weights")

  sum(matrix(p, nrow=2) %*% w)
}


#' Calculate the mean of mixture model
#'
#' @param x Mixture prior object
#' @param ... Ignored
#'
#' @return Mean of prior
#' @export
#'

mean.mixture.prior <- function(x, ...){
  
  m <- apply(matrix(attr(x,"pars"),nrow=2),2, function(y) y[1]/sum(y))

  w <- attr(x,"weights")

  sum(m*w)
}

#' Calculate the variance of mixture model
#'
#' @param mixture.prior Mixture prior object
#'
#' @return Variance of prior
#' @export
#'
var.mixture.prior <- function(mixture.prior){
  m <- apply(matrix(attr(mixture.prior,"pars"),nrow=2),2, function(x) x[1]/sum(x))
  v <- apply(matrix(attr(mixture.prior,"pars"),nrow=2),2, function(x) x[1]*x[2]/( sum(x)^2 * (sum(x)+1) ))
  w <- attr(mixture.prior,"weights")

  if(length(v)==1) return(v)

  mw <- m*w
  #by law of total variance for partitioned space
  sum(v*w) +
    sum(m^2*(1-w)*w) -
    2* sum( unlist(sapply(2:length(m), function(i) sapply(1:(i-1), function(j) mw[i]*mw[j]))))
}



#' Calculate posterior mixture model
#'
#' @param xs Number of successes in new trial
#' @param ns Number of patients in new trial
#' @param mixture.prior Mixture prior object
#'
#' @return A \code{mixture.prior} object
#' @export
#'
posterior.mixture.prior <- function(xs, ns, mixture.prior){
  pars <- matrix(attr(mixture.prior,'pars'),nrow=2)
   
  ps <- as.vector(pars) + c(xs,ns-xs)

  # s1 <- as.vector(pars[1,])
  # s2 <- as.vector(pars[2,])

  # lik <- mapply(dbetabinom.ab, x=xs, size=ns, shape1=s1, shape2=s2)
  lik <- dbetabinom.ab(x=xs, size=ns, shape1=as.vector(pars[1,]), shape2=as.vector(pars[2,]))
  ws <- attr(mixture.prior,'weights')*lik
  ws <- ws/sum(ws)

  flp <- update.mixture.prior(mixture.prior,
                         pars=ps,
                         weights= ws)
  return(flp)
}


#' Print mixture model
#'
#' @param x Mixture prior object
#' @param ... Ignored
#'
#' @export
#'
print.mixture.prior <- function(x, ...){
  
  w <- attr(x,"weights")
  p <- attr(x,"pars")
  par.type <- attr(x,"par.type")
  df <- as.data.frame(matrix(p, byrow = TRUE, ncol = length(par.type)))
  po <- rowSums(df)*w
  df <- cbind(df,w, po)
  colnames(df) <- c(par.type, "weights", "ESS")
  rownames(df) <- NULL
  cat("\n",attr(x,"type"),"mixture\n")
  print(df)
  cat("\nTotal sample size: ",ess.mixture.prior(x),"\n")
}


#' Sample from the mixture distribution
#'
#' @param n Number of samples
#' @param mixture.prior Mixture prior object
#'
#' @return A vector of n samples from \code{mixture.prior}
#' @export
#'
sample.mixture.prior <- function(n, mixture.prior){
  w <- attr(mixture.prior, "weights")
  p <- attr(mixture.prior, "pars")
  degree <- length(weights)
  each.dist <- sample(1:degree, size=n, replace=TRUE, prob=w)
  unlist(
    sapply(1:degree, function(d){
    rbeta(sum(each.dist==d), shape1=p[2*d-1], shape2=p[2*d])
  }))
}



#' Quantile function
#'
#' @param p Probability
#' @param mixture.prior Mixture prior object
#'
#' @return Vector of quantiles
#' @export
#'
q.mixture.prior <- function(p, mixture.prior){
  sapply(p, function(p){
    # optimise( function(q){(p.mixture.prior(q,mixture.prior)-p)^2},
              # interval=c(0,1))$minimum
    uniroot(function(q) (p.mixture.prior(q,mixture.prior)-p),
            interval=c(0,1))$root
  })
}



#' Probabilty function
#'
#' @param q Quantile
#' @param mixture.prior Mixture prior object
#'
#' @return Vectore of probabilities
#' @export
#'
p.mixture.prior <- function(q, mixture.prior){
  w <- attr(mixture.prior,"weights")
  L <- length(w)
  p <- split(attr(mixture.prior,"pars"),rep(1:L, each=2))

  sapply(q, function(q){
    sum(sapply(1:L, function(i){
      w[i]*pbeta(q, p[[i]][1],p[[i]][2])
    }))
  })
}



