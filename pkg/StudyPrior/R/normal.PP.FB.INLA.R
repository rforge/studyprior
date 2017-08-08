#' Normal power prior with full Bayes by INLA
#'
#' @param x numeric vector of results
#' @param sd standard deviations
#' @param verbose Print messages
#'
#' @return A density function
#' @export
#'

normal.PP.FB.INLA <- function(x, sd, verbose=FALSE ){
  n.hist <- length(x)

  A<-1
  B<-1
  prior.function <- function(theta) ifelse(theta < log(1),
                                           log(dbeta(exp(theta),A,B)*exp(theta)),
                                           -theta-10000)
  # ifelse(theta < log(1), exp(theta), -theta-10000)

  lprec = seq(-100, 10, len=1000)
  prior.table = paste(c("table:", cbind(lprec, prior.function(lprec))),
                      sep = "", collapse = " ")


  Y <- matrix(NA,n.hist,n.hist)
  diag(Y) <- x


  fam <- rep("gaussian", n.hist)
  con.fam <- c(rep(list(list(hyper=list(theta=list(prior=prior.table, initial=log(.1) )))), n.hist))

  sca <- 1/sd^2

  ifb <- inla(Y~1,
              data=list(Y=Y),
              family=fam,
              scale= sca,
              control.family = con.fam,
              control.inla = list(strategy='laplace', int.strategy='eb'),
              verbose=TRUE)

  summary(ifb)






  g <- function(p,X) ifelse(0<=p&p<=1, f(p),0)

  return(ifb)

}
