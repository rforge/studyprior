

#' Calculate bias
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
calc.bias <- function(prior, prob.range=c(.5,1), length=20, n.binom=30, posterior, mc.cores=1){

  P <- seq(prob.range[1],prob.range[2],len=length)
  if(missing(posterior)){
    Bias.for.x <- sapply(0:n.binom, function(Xs){

      if(inherits(prior, "function")){
        post <- function(p,g=1) prior(p,Xs)*dbinom(x=Xs, size=n.binom, prob=p)/g
        f <- splinefun(smooth.spline(seq(0.001,.999,len=1000), pmax(0,post(seq(0.001,.999,len=1000)))))
        # print(Xs)
        K <- adaptIntegrate(f, lower=0, upper=1, maxEval = 2e5)$integral

        # probability * square error
        Ep <- function(p, true.p) f(p)/K *p

        return(sapply(P, function(true.p){
          adaptIntegrate(Ep,0,1, true.p=true.p, maxEval=1e4)$integral - true.p})
        )

      } else if(inherits(prior, "mixture.list")){
        return(mean.fun.list(fun.list=posterior.fun.list(Xs, n.binom, prior))-P)

      } else if(inherits(prior, "list")){
        return(mean.fun.list(fun.list=posterior.fun.list(Xs, n.binom, prior[[Xs+1]]))-P)
      }
    })
  } else if(!missing(posterior)){
    Bias.for.x <- parallel::mclapply(posterior, function(post){
      print("Biasing!")
      print(post)
      Ep <-  adaptIntegrate(function(p) post(p)*p, 0,1, maxEval=1e4)$integral

      return(sapply(P, function(true.p){ Ep - true.p})
      )
      }, mc.cores = mc.cores)
  }

  Bias.for.x <- matrix(unlist(Bias.for.x), nrow=length)

  sapply(seq_along(P), function(i){
   sum(Bias.for.x[i,] * dbinom(0:n.binom, size=n.binom, prob=P[i]))})
}

