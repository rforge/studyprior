
#' Title
#'
#' @param n.control Number of subjects in the control group
#' @param n.treatment
#' @param level Significance level to test (default 0.975)
#' @param prior Prior on binomial probability parameter of the control group
#' @param treat.beta.prior.par Paramters for beta prior on treatment group probability parameter
#' @param mc.cores Number of cores to use for mclapply
#'
#' @return A matrix of TRUE/FALSE values of size n.control+1 x n.treatment+1
#' @export
#'

sig.matrix <- function(n.control, n.treatment, level=0.975, prior, treat.beta.prior.par=c(1,1), mc.cores=1, check.xt, check.xs,debug=FALSE) {
lapply
  ZZ.list <-   mclapply(check.xs,
    function(Xs){
      # print(Xs)
      post <-
        if(inherits(prior, "function")){
          post <- function(p,g=1) prior(p,Xs)*dbinom(x=Xs, size=n.control, prob=p)/g
          f <- splinefun(smooth.spline(seq(0,1,len=1000), pmax(0,post(seq(0,1,len=1000)))))

          K <- adaptIntegrate(f, lower=0, upper=1, maxEval = 2e5)$integral
          function(p,g=K) f(p)/g
          # formals(post) <- alist(p = , g = K)
          # K <- K *integrate(post, lower=0, upper=1)$value
        } else if(inherits(prior, "mixture.list")){
          post.list <- posterior.fun.list(Xs, n.control, prior)
          function(p) eval.fun.list(p, post.list)
        } else if(inherits(prior, "list")){
          post.list <- posterior.fun.list(Xs, n.control, prior[[Xs+1]])
          function(p) eval.fun.list(p, post.list)
        }
    ZZ <-unlist(lapply(check.xt, function(xT){
      # print(paste0('.',xT))
      res <-  try({
        if(debug ) browser()

        unsure <- TRUE
        #start quick and dirty
        tol <- 0.05

        this.int <- NA

        while(unsure){
          # this.int.res <-  adaptIntegrate(
          #   function(p0) {
          #     pmax(0,
          #         dbeta(p0,
          #               treat.beta.prior.par[1]+xT,
          #               treat.beta.prior.par[2] + n.treatment-xT) -
          #           post(p0))
          #     },
          #   lower=0,upper=1, rel.tol = tol, subdivisions = 200)

          this.int.res <-  adaptIntegrate(
            function(p0) {
              pbeta(p0,
                    treat.beta.prior.par[1]+xT,
                    treat.beta.prior.par[2] + n.treatment-xT,
                    lower.tail = FALSE) *
                post(p0)
            },
            lower=0,upper=1, tol = tol, maxEval = 2e5)
          this.int.res$value <- this.int.res$integral
          this.int.res$abs.error <- this.int.res$integral*this.int.res$error

          this.int <- this.int.res$value

          if( (this.int - this.int.res$abs.error) < level & (this.int + this.int.res$abs.error)  > level){
            # print("We were very close! Trying again -------------------")
            # print(this.int)
            tol <- tol/5
            # print(paste0("Tol: ",tol))
          } else {
            unsure <- FALSE
          }

        }

        this.int > level
      })
      if(inherits(res,"try-error")) browser()#save(xT, post, file=paste0('ERROR-2-',Xs,'-',xT,'.Rda'))
      return(res)
    }) )

    # print(ZZ)
# browser()
    ZZ
  }
  , mc.cores=mc.cores, mc.preschedule = FALSE
  )

  matrix(unlist(ZZ.list),
         byrow=TRUE,
         nrow=n.control+1,
         ncol=n.treatment+1,
         dimnames = list(Control = 0:n.control,
                         Treatment = 0:n.treatment)
  )

}
