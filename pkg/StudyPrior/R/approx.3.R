
#' Approximate distribution with mixture model
#'
#' @param distr
#' @param type
#' @param degree
#' @param range
#' @param starts
#' @param length.fit
#'
#' @return
#' @export
#'
#' @examples
conj.approx2 <- function(distr, type=c("beta","normal"), min.degree=max.degree, max.degree = 3, return.value=0.2, range=c(0,1), starts, length.fit=100, do.plot=FALSE, robust=FALSE){

  x <- rep(seq(range[1],range[2], length.out = length.fit), each=1)
  y <- distr(x) #+ rnorm(length(x), 0, 0.05)
  dat <- data.frame(x,y)[-c(1:4,(length.fit-4):length.fit),]

  dat.quick <- data.frame(  x = rep(seq(range[1],range[2], length.out = length.fit), each=1),
                            y = distr(x))[-c(1:3,(23):25),]

  if(do.plot) plot(x,y, type='l', col=2)

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
    fl <- create.fun.list(type,
                          unlist(start.list)[-c(1:degree)],
                          unlist(start.list)[1:degree])


    opt.env <- new.env()
    assign("x", x, envir=opt.env)
    assign("y", y, envir=opt.env)
    assign("fun.list", fl, envir=opt.env)

    if(missing(starts)) starts <- start.list

    opt <-
      optimr::optimr(unlist(starts),
                     function(PAR){
                       sum((dat$y-eval.fun.list(dat$x, update.fun.list(fun.list = fl,
                                                                                   pars=PAR[-(1:degree)],
                                                                                   weights=PAR[1:degree])))^2)
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
                         sum((dat$y-eval.fun.list(dat$x, update.fun.list(fun.list = fl,
                                                                         pars=PAR[-(1:degree)],
                                                                         weights=PAR[1:degree])))^2)
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
                         sum((dat$y-eval.fun.list(dat$x, update.fun.list(fun.list = fl,
                                                                         pars=PAR[-(1:degree)],
                                                                         weights=PAR[1:degree])))^2)
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

if(robust){
  fl <- create.fun.list(type, pars=c(opt$par[-(1:degree)],1,1), weights=c(opt$par[1:degree],0.0001))
} else fl <- create.fun.list(type, pars=opt$par[-(1:degree)], weights=opt$par[1:degree])

  if(do.plot) plot.fun.list(x, fl, stack=TRUE, lines.only=TRUE)

  return(fl)
}
