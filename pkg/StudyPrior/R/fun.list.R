#' Create mixture model
#'
#' @param type
#' @param pars
#' @param weights
#'
#' @return
#' @export
#'
#' @examples
create.fun.list <- function(type, pars, weights){
  degree <- length(weights)

  fl <- rep(list(switch(type, "normal" = dnorm,"beta" = dbeta)),
            degree)

  attr(fl,"type") <- type
  attr(fl,"par.type") <- switch(type,
                                "normal" = c("mean","sd"),
                                "beta" = c("shape1","shape2") )

  fl <- update.fun.list(fl, pars, weights)
  class(fl) <- "mixture.list"
  return(fl)
}


#' Update parameters or weights in mixture model
#'
#' @param fun.list
#' @param pars
#' @param weights
#'
#' @return
#' @export
#'
#' @examples
update.fun.list <- function(fun.list, pars, weights ){
  if(!missing(pars)){
    par.names <- attr(fun.list, "par.type")
    split.pars <- split(unlist(pars), rep(1:length(fun.list), each=length(par.names)))
    attr(fun.list,"pars") <- pars
    for(i in 1:length(fun.list)){
      tmp <- formals(fun.list[[i]])
      tmp[par.names] <- split.pars[[i]]
      formals(fun.list[[i]]) <- tmp
    }
  }
  if(!missing(weights)){
    attr(fun.list,"weights") <- weights/sum(weights)
  }

  return(fun.list)
}



#' Get density of mixture
#'
#' @param x
#' @param fun.list
#' @param weights
#' @param subset
#'
#' @return
#' @export
#'
#' @examples
eval.fun.list <- function(x, fun.list, weights, subset){
  if(missing(subset)) subset <- 1:length(fun.list)
  X <- rep(0,length(x))

  if(missing(weights)) weights <- attr(fun.list, "weights")
  w <- weights/sum(weights)

  for(i in subset){
    X <- X + w[i]*do.call(fun.list[[i]], list(x=x))
  }
  X
}

#' Plot mixture model
#'
#' @param x
#' @param fun.list
#' @param stack
#' @param lines.only
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.fun.list <- function(x, fun.list, stack=FALSE, lines.only=FALSE, ...){

  Y <- eval.fun.list(x, fun.list)
  if(!lines.only) plot(x,Y, type='n', ...)

  if(stack==TRUE){
    degree <- length(fun.list)
    w <- attr(fun.list, "weights")
    z <- order(w, decreasing = TRUE)

    for(i in 1:(degree-1)){

        lines(x,eval.fun.list(x, fun.list,subset=z[1:i]), type='l', ...)
    }
  }
  lines(x,Y, type='l')
}

#' Calculate the effective sample size of mixture model
#'
#' @param fun.list
#'
#' @return
#' @export
#'
#' @examples
ess.fun.list <- function(fun.list){
  p <- attr(fun.list,"pars")
  w <- attr(fun.list,"weights")

  sum(matrix(p, nrow=2) %*% w)
}


#' Calculate the mean of mixture model
#'
#' @param fun.list
#'
#' @return
#' @export
#'
#' @examples
mean.fun.list <- function(fun.list){
  m <- apply(matrix(attr(fun.list,"pars"),nrow=2),2, function(x) x[1]/sum(x))

  w <- attr(fun.list,"weights")

  sum(m*w)
}

#' Calculate the variance of mixture model
#'
#' @param fun.list
#'
#' @return
#' @export
#'
#' @examples
var.fun.list <- function(fun.list){
  m <- apply(matrix(attr(fun.list,"pars"),nrow=2),2, function(x) x[1]/sum(x))
  v <- apply(matrix(attr(fun.list,"pars"),nrow=2),2, function(x) x[1]*x[2]/( sum(x)^2 * (sum(x)+1) ))
  w <- attr(fun.list,"weights")

  if(length(v)==1) return(v)

  mw <- m*w
  #by law of total variance for partitioned space
  sum(v*w) +
    sum(m^2*(1-w)*w) -
    2* sum( unlist(sapply(2:length(m), function(i) sapply(1:(i-1), function(j) mw[i]*mw[j]))))
}



#' Calculate posterior mixture model
#'
#' @param xs
#' @param ns
#' @param fun.list
#'
#' @return
#' @export
#'
#' @examples
posterior.fun.list <- function(xs, ns, fun.list){
  pars <- attr(fun.list,'pars')
  weights <- attr(fun.list,'weights')
  degree <- length(weights)

  ps <- as.vector(matrix(pars,nrow=2) + c(xs,ns-xs))

  s1 <- as.vector(matrix(pars,nrow=2)[1,])
  s2 <- as.vector(matrix(pars,nrow=2)[2,])

  lik <- mapply(VGAM::dbetabinom.ab, x=xs, size=ns, shape1=s1, shape2=s2)
  ws <- weights*lik / sum(weights*lik)

  flp <- update.fun.list(fun.list,
                         pars=ps,
                         weights= ws)
  return(flp)
}


#' Print mixture model
#'
#' @param fun.list The mixture model to print
#'
#' @return
#' @export
#'
print.fun.list <- function(fun.list){
  w <- attr(fun.list,"weights")
  p <- attr(fun.list,"pars")
  par.type <- attr(fun.list,"par.type")
  df <- as.data.frame(matrix(p, byrow = TRUE, ncol = length(par.type)))
  po <- rowSums(df)*w
  df <- cbind(df,w, po)
  colnames(df) <- c(par.type, "weights", "ESS")
  rownames(df) <- NULL
  cat("\n",attr(fun.list,"type"),"mixture\n")
  print(df)
  cat("\nEffective sample size: ",ess.fun.list(fun.list))
}


#' Sample from the mixture distribution
#'
#' @param n number of samples
#' @param fun.list
#'
#' @return A vector of n observations
#' @export
#'
sample.fun.list <- function(n, fun.list){
  w <- attr(fun.list, "weights")
  p <- attr(fun.list, "pars")
  degree <- length(weights)
  each.dist <- sample(1:degree, size=n, replace=TRUE, prob=w)
  unlist(
    sapply(1:degree, function(d){
    rbeta(sum(each.dist==d), shape1=p[2*d-1], shape2=p[2*d])
  }))
}



#' Quantile function
#'
#' @param p
#' @param fun.list
#'
#' @return
#' @export
#'
#' @examples
q.fun.list <- function(p, fun.list){
  sapply(p, function(p){
    optimise( function(q){(p.fun.list(q,fun.list)-p)^2},
              interval=c(0,1))$minimum
  })
}



#' Probabilty density function
#'
#' @param q
#' @param fun.list
#'
#' @return
#' @export
#'
#' @examples
p.fun.list <- function(q, fun.list){
  w <- attr(fun.list,"weights")
  L <- length(w)
  p <- split(attr(fun.list,"pars"),rep(1:L, each=2))

  sapply(q, function(q){
    sum(sapply(1:L, function(i){
      w[i]*pbeta(q, p[[i]][1],p[[i]][2])
    }))
  })
}
