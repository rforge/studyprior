
#' Approximate distribution with mixture model
#'
#' @param distr The density function of the function to be approximated
#' @param type Use beta or normal densities in the mixture
#' @param range The range of the density function to approximate
#' @param starts Initial parameters for the search
#' @param length.fit Number of points to fit to
#' @param min.degree Minimum size of mixture
#' @param max.degree Maximum size of mixture
#' @param robust Include a robust component? (FALSE or final weight between 0 and 1)
#' @param do.plot Show the components of the mixture in a plot
#' @param return.value Convergence threshold
#'
#' @return Returns a mixture.prior object
#' @export
#'

conj.approx <- function(distr,
                        type=c("beta","normal"),
                        max.degree = 3,
                        min.degree=max.degree,
                        robust=FALSE,
                        range=c(0,1),
                        do.plot=FALSE,
                        return.value=0.2,
                        starts,
                        length.fit=100
                        ){

  x <- rep(seq(range[1],range[2], length.out = length.fit), each=1)
  y <- distr(x)
  dat <- data.frame(x,y)[-c(1:4,(length.fit-4):length.fit),]

  dat.quick <- data.frame(  x = rep(seq(range[1],range[2], length.out = length.fit), each=1),
                            y = distr(x))[-c(1:3,(23):25),]

  plot.x <- rep(seq(range[1],range[2], length.out = length.fit*4), each=1)
  plot.y <- distr(plot.x) 
  if(do.plot) plot(plot.x, plot.y, type='l', col=2, xlab=expression(theta), ylab="Density")

  type <- match.arg(type)

  ddistr <- switch(type,
                   "beta" = dbeta,
                   "normal" = dnorm )

  values <- rep(9999,length=max.degree)
  results <- list()

  degrees <- seq(min.degree,max.degree)
  for( degree in degrees){

    params <- paste0(
      switch(type,"normal" = c("mean","sd"),"beta" = c("shape1","shape2")),
      rep(1:degree, each = 2))
    ws <- paste0("w", rep(1:degree))

    #############################################################################

    start.list <- numeric(length(c(ws,params)))
    names(start.list) <- c(ws,params)
    start.list[1:degree] <- degree:1/sum(1:degree) #start with equal weights

    mode.d <- x[which.max(y)]

    fudge <- if(mode.d!=0) round(1/mode.d) else 1

    # if(missing(starts)){
      start.list[-c(1:degree)] <-
        switch(type,
               "beta" = c(rbind(mode.d*(fudge+1:degree)/(1-mode.d), fudge + 1:degree)),
               "normal" = c(rbind(mode.d+rnorm(degree,0,0.1),1:degree/degree/10 )))
    # } else start.list <- starts
    #############################################################################
    lower.list <-  switch(type,
                          "beta" = c(rep(0.0075, degree), rep(0.001, length(params))),
                          "normal" =c(rep(0, degree),rep(c(0,0),degree)))

    #############################################################################

    upper.list <-  switch(type,
                          "beta" = c(rep(1, degree), rep(1000, length(params))),
                          "normal" =c(1/1:degree,rep(c(1,5),degree)))

    #############################################################################
    fl <- create.mixture.prior(type,
                          pars = matrix(unlist(start.list)[-c(1:degree)],ncol=2),
                          weights = unlist(start.list)[1:degree])


    opt.env <- new.env()
    assign("x", x, envir=opt.env)
    assign("y", y, envir=opt.env)
    assign("mixture.prior", fl, envir=opt.env)

    if(missing(starts)) starts <- start.list

    opt <-
      optimr::optimr(unlist(starts),
                     function(PAR){
                       sum((dat$y-eval.mixture.prior(dat$x, update.mixture.prior(object = fl,
                                                                                   pars=matrix(PAR[-(1:degree)],ncol=2),
                                                                                   weights=PAR[1:degree])))^2)/length.fit
                     },
                     lower=unlist(lower.list),
                     upper=unlist(upper.list),
                     method = "L-BFGS-B",
                     # method = "Rvmmin",
                     control=list(trace=FALSE,
                                  maxit=3000)
      )

    if(opt$convergence!=0){

      opt <-
        optimr::optimr(unlist(start.list),
                       function(PAR){
                         sum((dat$y-eval.mixture.prior(dat$x, update.mixture.prior(object = fl,
                                                                         pars=matrix(PAR[-(1:degree)],ncol=2),
                                                                         weights=PAR[1:degree])))^2)/length.fit
                       },
                       lower=unlist(lower.list),
                       upper=unlist(upper.list),
                       method = "L-BFGS-B",
                       # method = "Rvmmin",
                       control=list(trace=3,
                                    maxit=3000)
        )
    }

    if(opt$convergence!=0){
      start.rand <- mapply(runif, min=unlist(lower.list),max=unlist(upper.list), n=1)
      opt <-
        optimr::optimr(start.rand,
                       function(PAR){
                         sum((dat$y-eval.mixture.prior(dat$x, update.mixture.prior(object = fl,
                                                                         pars=matrix(PAR[-(1:degree)],ncol=2),
                                                                         weights=PAR[1:degree])))^2)/length.fit
                       },
                       lower=unlist(lower.list),
                       upper=unlist(upper.list),
                       method = "L-BFGS-B",
                       # method = "Rvmmin",
                       control=list(trace=3,
                                    maxit=3000)
        )
    }

    # print(paste(degree, ": ", opt$value))

    if(opt$convergence==0){
      results[[degree]] <- opt
      values[[degree]] <- opt$value
      #if we're smaller than the threshold return immediately
      if(opt$value  < return.value) break
    } else{
      # browser()
      results[[degree]] <- opt
      values[[degree]] <- .Machine$double.xmax
    }

  }#end the for loop over degrees

  # browser()
  id <- which.min(values)
  opt <- results[[id]]
  degree <- (1:max.degree)[id]

  #add robust component
  if(robust){
    fl <- create.mixture.prior(type, 
                               pars=rbind(matrix(opt$par[-(1:degree)],ncol=2), c(1,1)),
                               weights=c(opt$par[1:degree]/sum(c(opt$par[1:degree]))*(1-robust),
                                         robust))
  } else fl <- create.mixture.prior(type, 
                                    pars=matrix(opt$par[-(1:degree)],ncol=2),
                                    weights=opt$par[1:degree])
  
  #do plot if necessary
  if(do.plot) plot.mixture.prior(plot.x, fl, stack=TRUE, lines.only=TRUE)

  return(fl)
}

conj.approx2 <- conj.approx