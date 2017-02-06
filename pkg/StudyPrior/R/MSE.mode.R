

#' Calculate mean squared error based on the mode of the posterior
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
calc.MSE.mode <- function(prior, prob.range=c(.5,1), length=20, n.binom=30, mc.cores=1){

  P <- seq(prob.range[1],prob.range[2],len=length)




  MSE.for.x <- parallel::mclapply(0:n.binom, function(Xs){


    if(inherits(prior, "function")){
      post <- function(p,g=1) prior(p,Xs)*dbinom(x=Xs, size=n.binom, prob=p)/g
      p <- optimise(post, interval=c(0,1), maximum=TRUE)

      #   square error
      sq.err <- (p$maximum-P)^2

      return(sq.err)

    } else if(inherits(prior, "mixture.list")){
      post.list <- posterior.fun.list(Xs, n.binom, prior)
      return((mean.fun.list(post.list)-P)^2 + var.fun.list(post.list))
    } else if(inherits(prior, "list")){
      post.list <- posterior.fun.list(Xs, n.binom, prior[[Xs+1]])
      return((mean.fun.list(post.list)-P)^2 + var.fun.list(post.list))
    }


  }, mc.cores = mc.cores)

  # return(MSE.for.x)

  MSE.for.x <- matrix(unlist(MSE.for.x), nrow=length)


  sapply(seq_along(P), function(i){
    sum(MSE.for.x[i,] * dbinom(0:n.binom, size=n.binom, prob=P[i]))})
}

