#' Plot the mixture model
#'
#' @param pars Parameters of each distribution component
#' @param weights Vector of weights for each component
#' @param dist Distribution
#'
#' @return
#' @export
#'
#' @examples
mix.plot <- function(pars, weights, dist, split.f){
  degree <- length(weights);

  if(missing(split.f)){
    if(length(pars)%% degree !=0) stop("Don't know how to apply parameters")
    par.split <- split(pars, rep(1:degree, each=length(pars)/degree))
  } else par.split <- split(pars, split.f)

  component <- function(d, par.list) as.call(c(list(d), as.list(par.list)))

  calls <- lapply(1:degree, function(i) component(d=dist, par.list = par.split[[i]] ))
x <- seq(0,1,len=11)
  lapply(1:degree, function(i) eval(calls[[i]],envir=environment(x)))

}
