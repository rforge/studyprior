#' #' Extension of the binomial distribution to allow real values for size and x
#' #'
#' #' @param x
#' #' @param size
#' #' @param prob
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' dbinomc <- function(x, size, prob){
#'   if(x>size) stop("x greater than size")
#'   if(x<0) stop("x must be positive")
#' 
#'   A <- x*log(prob) + (size-x)*log(1-prob)
#'   B <- lgamma(size+1)-lgamma(x+1)-lgamma(size-x+1)
#' 
#'   logdens <- A+B
#'   return(exp(logdens))
#' 
#' }
