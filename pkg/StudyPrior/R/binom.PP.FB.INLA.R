#' #' Title
#' #'
#' #' @param x
#' #' @param n
#' #' @param verbose
#' #' @param length
#' #' @param d.prior.a shape1 parameter for beta prior on weight
#' #' @param d.prior.b shape2 parameter for beta prior on weight
#' #' @param mc.cores
#' #' @param samples
#' #' @param focus
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' binom.PP.FB.INLA <- function(x, n, X,N,verbose=FALSE, length=30, d.prior.a=1, d.prior.b=1, mc.cores=1){
#'   n.hist <- length(x)
#' 
#' A<-1
#' B<-1
#'   prior.function <- function(theta) ifelse(theta < log(1),
#'                                            log(dbeta(exp(theta),A,B)*exp(theta)),
#'                                            -theta-10000)
#'   # ifelse(theta < log(1), exp(theta), -theta-10000)
#' 
#'   lprec = seq(-100, 10, len=1000)
#'   prior.table = paste(c("table:", cbind(lprec, prior.function(lprec))),
#'                       sep = "", collapse = " ")
#' 
#'   # Y <- matrix(NA,n.hist+1+1, n.hist+1+1)
#'   Y <- matrix(NA,n.hist,n.hist)#+1, n.hist+1)
#'   diag(Y) <-c(x/n)#, X/N) #-log(1/c(x/n,0.5)-1)
#' 
#' 
#'   fam <- rep("beta", n.hist)#+1)# +1)
#'   con.fam <- c(rep(list(list(hyper=list(theta=list(prior=prior.table, initial=log(.1) )))), n.hist)
#'                # ,
#'                # list(list(hyper=(list(theta=list(fixed=TRUE, initial=0)))))
#'                )
#'   # ,               list(list(hyper=(list(theta=list(fixed=TRUE, initial=0))))))
#' 
#'   sca <- c(n)#, 2*n.hist+2)#,N)
#' 
#'   ifb <- inla(Y~1,
#'               data=list(Y=Y),
#'               family=fam,
#'               scale= sca,
#'               control.family = con.fam,
#'               control.inla = list(strategy='laplace', int.strategy='eb'),
#'               verbose=TRUE)
#' 
#'   summary(ifb)
#' 
#' 
#' 
#' 
#' 
#' 
#'   g <- function(p,X) ifelse(0<=p&p<=1, f(p),0)
#' 
#'   return(ifb)
#' 
#' }
