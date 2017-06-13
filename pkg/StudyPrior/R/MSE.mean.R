

#' Calculate mean squared error based on the mean of the posterior
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
calc.MSE.mean <- function(prior, prob.range=c(.5,1), length=20, n.binom=30, mc.cores=1, posterior){

  P <- seq(prob.range[1],prob.range[2],len=length)

  if(missing(posterior)){
  MSE.for.x <- parallel::mclapply(0:n.binom, function(Xs){
    if(inherits(prior, "function")){
      post <- function(p,g=1) prior(p,Xs)*dbinom(x=Xs, size=n.binom, prob=p)/g
      # G <-  adaptIntegrate(post, lower=0, upper=1, maxEval = 2e5)$integral

      f <- splinefun(smooth.spline(seq(0.001,0.999,len=1000), pmax(0,post(seq(.001,.999,len=1000)))))
      # print(Xs)
      K <- adaptIntegrate(f, lower=0, upper=1, maxEval = 2e5)$integral

      post.mean <- adaptIntegrate(function(p) p*f(p)/K,0,1,  maxEval=2e5)$integral
        return((post.mean - P)^2)

    } else if(inherits(prior, "mixture.list")){
      post.list <- posterior.fun.list(Xs, n.binom, prior)
      return((mean.fun.list(post.list)-P)^2 )
    } else if(inherits(prior, "list")){
      post.list <- posterior.fun.list(Xs, n.binom, prior[[Xs+1]])
      return((mean.fun.list(post.list)-P)^2 )
    }
  }, mc.cores = mc.cores)
  } else if(!missing(posterior)){
    MSE.for.x <- parallel::mclapply(posterior, function(post){
      post.mean <- adaptIntegrate(function(p) p*post(p),0,1,  maxEval=2e5)$integral
      return((post.mean - P)^2)
      }, mc.cores = mc.cores)
  }

  # return(MSE.for.x)

  MSE.for.x <- matrix(unlist(MSE.for.x), nrow=length)


  sapply(seq_along(P), function(i){
    sum(MSE.for.x[i,] * dbinom(0:n.binom, size=n.binom, prob=P[i]))})
}

