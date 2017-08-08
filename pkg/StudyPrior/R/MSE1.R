#' Calculate mean squared error
#'
#' @param prior Prior to calculate posterior. Specify posterior instead if available.
#' @param prob.range Range of values to calculate MSE over
#' @param length Number of values to calculate MSE for
#' @param n.binom Number of patients in new trial
#' @param mc.cores Number of cores for parallel
#' @param posterior Posterior density 
#'
#' @return A vector of error values
#' @export
#'
calc.MSE <- function(prior, prob.range=c(.5,1), length=20, n.binom=30, mc.cores=1, posterior){

  P <- seq(prob.range[1],prob.range[2],len=length)

  if(missing(prior) &  missing(posterior)) stop("prior or posterior must be specified")
  
  if(missing(posterior)){
    MSE.for.x <- parallel::mclapply(0:n.binom, function(Xs){
      if(inherits(prior, "function")){
        post <- function(p,g=1) prior(p,Xs)*dbinom(x=Xs, size=n.binom, prob=p)/g
        f <- splinefun(smooth.spline(seq(0.001,0.999,len=1000), pmax(0,post(seq(.001,.999,len=1000)))))
        # print(Xs)
        K <- adaptIntegrate(f, lowerLimit = 0, upperLimit = 1, maxEval = 2e5)$integral
        #  probability * square error
        sq.err <- function(p, true.p) f(p)/K * (p-true.p)^2

        return(sapply(P, function(true.p){
          adaptIntegrate(sq.err,0,1, true.p=true.p, maxEval=2e5)$integral})
        )

      } else if(inherits(prior, "mixture.list")){
        post.list <- posterior.mixture.prior(Xs, n.binom, prior)
        return((mean.mixture.prior(post.list)-P)^2 + var.mixture.prior(post.list))
      } else if(inherits(prior, "list")){
        post.list <- posterior.mixture.prior(Xs, n.binom, prior[[Xs+1]])
        return((mean.mixture.prior(post.list)-P)^2 + var.mixture.prior(post.list))
      }


    }, mc.cores = mc.cores)
  } else if(!missing(posterior)){
    MSE.for.x <- parallel::mclapply(posterior, function(post){
        sq.err <- function(p, true.p) post(p) * (p-true.p)^2
        return(sapply(P, function(true.p){
          adaptIntegrate(sq.err,0,1, true.p=true.p, maxEval=2e5)$integral})
        )}, mc.cores = mc.cores)
  }


  MSE.for.x <- matrix(unlist(MSE.for.x), nrow=length)


  sapply(seq_along(P), function(i){
   sum(MSE.for.x[i,] * dbinom(0:n.binom, size=n.binom, prob=P[i]))})
}

