#' Mixture distribution
#'
#' @param distr
#' @param conj.distr
#' @param degree
#' @param range
#'
#' @return
#' @export
#'
#' @examples
conj.approx <- function(distr, type=c("beta","normal"), degree = 3, range=c(0,1), starts){

    x <- rep(seq(range[1],range[2], length.out = 100), each=1)
  y <- distr(x) #+ rnorm(length(x), 0, 0.05)
  dat <- data.frame(x,y)[-c(1,100),]
  plot(x,y, type='l', col=2)

  type <- match.arg(type)

  ddistr <- switch(type,
                   "beta" = "dbeta",
                   "normal" = "dnorm" )


  weights <- switch(as.character(degree),
                    "1" = "w1",
                    "2" = c("w1/(w1+w2)","w2/(w1+w2)"),
                    "3" = c("w1/(w1+w2+w3)","w2/(w1+w2+w3)","w3/(w1+w2+w3)"))

  ws <- paste0("w",1:degree)
  params <- paste0(c("a","b"),rep(1:degree,each=2))

  optim.formula <-
  paste(optim.formula <- "y ~ ",
        paste0(
          sapply(1:degree, function(i){
            paste0(weights[i],"*",ddistr,"(x,",
                   paste(params[(2*i-1):(2*i)], collapse=","),
                   ")"
            )}),
          collapse=" + "))

  start.list <- numeric(length(c(ws,params)))
  names(start.list) <- c(ws,params)
  start.list[1:degree] <- degree:1/sum(1:degree) #start with equal weights

  mode.d <- x[which.max(y)]

  fudge <- round(1/mode.d)

  if(missing(starts)){
    start.list[-c(1:degree)] <-
      switch(type,
             "beta" = c(rbind(mode.d*(fudge+1:degree)/(1-mode.d), fudge + 1:degree)),
             "normal" = c(rbind(mode.d+rnorm(degree,0,0.1),1:degree/degree/10 )))
  }

  lower.list <-  switch(type,
                        "beta" = c(rep(0, degree), rep(1, length(params))),
                        "normal" =c(rep(0, degree),rep(c(0,0),degree)))



  upper.list <-  switch(type,
                        "beta" = c(rep(1, degree), rep(100, length(params))),
                        "normal" =c(1/1:degree,rep(c(1,5),degree)))

  fit <- minpack.lm::nlsLM(formula(optim.formula),
              data=dat,
              start=start.list,
              lower=lower.list,
              upper=upper.list,
              trace = TRUE,
              # algorithm = 'port',
              control=nls.control(max=100)
              )

  # opt2 <-
  #   nlmrt::nlxb(formula(optim.formula),
  #               data=data.frame(x,y),
  #               start=unlist(start.list),
  #               lower=unlist(lower.list),
  #               upper=unlist(upper.list),
  #               trace = TRUE)

  #
  # optim(start.list,
  #   function(X) sum((y -  X[1]/(X[1]+X[2])*dbeta(x,X[3],X[4]) + X[2]/(X[1]+X[2])*dbeta(x,X[5],X[6]))^2),
  #   #(y -  w1/(w1+w2)*dbeta(x,a1,b1) + w2/(w1+w2)*dbeta(x,a2,b2))^2
  #       lower = lower.list,
  #   upper = upper.list,
  #   method="L-BFGS-B"
  # )

  lines(dat$x,fit$m$predict(dat$x))

  p <- fit$m$getPars()
  p[1:degree] <- p[1:degree]/sum(p[1:degree])
  fit$m$setPars(newPars=p)

  fun.list <- create.fun.list(type, pars=p[-c(1:degree)], weights=p[c(1:degree)])

  return(list(fun.list, fit))
}
