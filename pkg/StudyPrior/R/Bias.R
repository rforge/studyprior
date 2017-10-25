#TODO ADD SUPPORT FOR MODE TO PRIOR ONLY CALCULATIONS

#' Calculate bias
#'
#' @param prior Prior to calculate posterior. Specify posterior instead if available.
#' @param prob.range Range of values to calculate MSE over
#' @param length Number of values to calculate MSE for
#' @param n.binom Number of patients in new trial
#' @param mc.cores Number of cores for parallel
#' @param posterior Posterior density 
#' @param type Either "mean" or "mode" of the posterior to use as the estimate
#'
#' @return A vector of bias values
#' @export
#'
calc.bias <- function(prior, prob.range=c(.5,1), length=20, n.binom=30, posterior, mc.cores=1, type="mean"){


  P <- seq(prob.range[1],prob.range[2],len=length)
  if(missing(posterior)){
    Bias.for.x <- parallel::mclapply(0:n.binom, function(Xs){

      if(inherits(prior, "function")){
        post <- function(p,g=1) prior(p,Xs)*dbinom(x=Xs, size=n.binom, prob=p)/g
        f <- splinefun(smooth.spline(seq(0.001,.999,len=1000), pmax(0,post(seq(0.001,.999,len=1000)))))
        # print(Xs)
        K <- adaptIntegrate(f, lowerLimit = 0, upperLimit = 1, maxEval = 2e5)$integral

        # probability * square error
        Ep <- function(p, true.p) f(p)/K *p

        return(sapply(P, function(true.p){
          adaptIntegrate(Ep,0,1, true.p=true.p, maxEval=1e4)$integral - true.p})
        )

      } else if(inherits(prior, "mixture.prior")){
        return(mean.mixture.prior(x=posterior.mixture.prior(Xs, n.binom, prior))-P)

      } else if(inherits(prior, "list")){
        return(mean.mixture.prior(x=posterior.mixture.prior(Xs, n.binom, prior[[Xs+1]]))-P)
      }
    }, mc.cores=mc.cores)
    
  } else if(!missing(posterior)){
    
    if(inherits(posterior[[1]], "mixture.prior")){
      #for a list of mixtures
      Bias.for.x <- parallel::mclapply(posterior, function(post){
        return((mean(post)-P) ) #mean.mixture.prior()
      }, mc.cores = mc.cores)
    } else {
      #Methods for a list of posterior functions  
      Bias.for.x <- parallel::mclapply(posterior, function(post){
        print("Biasing!")
        print(post)
        if( type == "mean") Ep <-  adaptIntegrate(function(p) post(p)*p, 0,1, maxEval=1e4)$integral
        else if(type == "mode") Ep <-  optimize(function(p) post(p), c(0,1), maximum = TRUE)$maximum
        
        return(sapply(P, function(true.p){ Ep - true.p})
        )
      }, mc.cores = mc.cores)
    }
  }

  Bias.for.x <- matrix(unlist(Bias.for.x), nrow=length)

  sapply(seq_along(P), function(i){
   sum(Bias.for.x[i,] * dbinom(0:n.binom, size=n.binom, prob=P[i]))})
}

