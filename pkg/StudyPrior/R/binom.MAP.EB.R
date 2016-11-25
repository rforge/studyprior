#' Empirical Bayes Meta-Analytic Prior for Binomial Data
#'
#' @param x historical events
#' @param n historical trials
#' @param tau.prior optional prior on heterogeneity parameter
#'
#' @return A function of the probability parmater p
#' @export
#'
#' @examples
#'
binom.MAP.EB <- function(x, n, X, N, verbose=FALSE, upper, VAR, mc.cores){
  n.hist <- length(x)
  n.new <- 0

 f <- mcmapply(mc.cores=mc.cores,
               X=X, N=N,
    FUN=function(X,N, upper=upper, VAR=VAR, verbose=verbose){

      if(!(missing(X)|missing(N))) {
        x <- c(x,X)
        n <- c(n,N)
        n.new <- 1
      }

      dat <- data.frame(x=x, n=n, z=1:(n.hist+n.new))



      if(missing(upper)) upper <- 4

      logiflatAB <- paste0("expression:
  upper = ",upper,";
  sd = exp(-x/2);
  A = (sd > upper ) ? 0.0000001 : 1/upper;
  B = abs(exp(-x/2)/2);
  logdens = log(A) + log(B);
  return(logdens)", collapse="")

      prior <-  list(prior= logiflatAB, initial=-3)
      # prior <-  list(prior= "logtnormal", param=c(0,1))
      # prior <-  list(prior= "flat", param=numeric(0))


      formula <- x ~ 1 + f(z, model="iid",
                           hyper = list(theta = prior))
      ##########################################


      if(missing(VAR)){
        result <- INLA::inla(formula,
                             data = dat,
                             family = "binomial",
                             control.fixed = list(mean.intercept = 0, prec.intercept = 1/1000),
                             Ntrials=n
                             # ,
                             # control.predictor = list(compute=TRUE, link=1)
                             )
        # ,
                             # verbose = verbose,
                             # control.inla = list(int.strategy = "eb"))

        result <- INLA::inla.hyperpar(result)

        # plot(INLA::inla.tmarginal(function(x) 1/x^0.5,
        #                           result$marginals.hyperpar[[1]],
        #                           n=200), xlim=c(0,0.0001))

        mode_tau <-  inla.mmarginal(INLA::inla.tmarginal(function(x) 1/x^.5,
                                                         result$marginals.hyperpar[[1]],
                                                         n=2000)[50:1950,])
        #
        #   mode_tau <- optimize(function(x){
        #     INLA::inla.dmarginal(x,
        #                    marg=INLA::inla.tmarginal(function(x) x^-.5,
        #                                        result$marginals.hyperpar[[1]],
        #                                        n=2000))},
        #   interval = c(0,upper),
        #   maximum = TRUE)

        # print(mode_tau)
      }
      # VXN <- ifelse(missing(VAR), var(log((x/n)/(1-x/n))), VAR)
      VXN <- mode_tau^2

      formulaEB = x ~ 1 + f(z, model="iid",
                            hyper = list(theta = list(fixed=TRUE,
                                                      initial=log(1/VXN)#1/mode_tau$maximum
                            )))

      dat <- data.frame(x=c(x[1:n.hist],NA), n=c(n[1:n.hist],NA), z=1:(n.hist+1))

      resultEB = INLA::inla(formulaEB,
                            data = dat,
                            family = "binomial",
                            control.fixed = list(mean.intercept = 0, prec.intercept = 1/1000),
                            Ntrials=n,
                            control.predictor = list(compute=TRUE, link=1))

      #   plot( resultEB$marginals.fixed$`(Intercept)`)
      #   points( resultEB$marginals.fitted.values$fitted.predictor.4)
      #   plot( resultEB$marginals.fitted.values$fitted.predictor.4)
      #
      #
      #   # for some reason the link is wrong
      #   plot(inla.tmarginal(function(x) 1/(1+exp(-x)), resultEB$marginals.fixed$`(Intercept)`),xlim=c(0,1))
      # points(inla.tmarginal(function(x) 1/(1+exp(-x)), resultEB$marginals.fitted.values$fitted.predictor.4))
      ind <- round(resultEB$marginals.fitted.values[[n.hist+1]][,1],7)%in% c(0,1)

      A <- resultEB$marginals.fitted.values[[n.hist+1]][!ind,1]
      B <- resultEB$marginals.fitted.values[[n.hist+1]][!ind,2]

      f <- INLA::inla.dmarginal

      formals(f) <- list(x=NULL,
                         marginal=list(x=A,y=B),
                         log=FALSE)
      f
})

  g <- function(p,X) f[[X+1]](p)
  return(g)
}
