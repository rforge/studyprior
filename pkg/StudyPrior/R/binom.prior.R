

#' Prior distribution for binomial data
#'
#' @param x Vector of historical events
#' @param n Vector of historical sample sizes
#' @param type Which prior to provide
#' @param ... Parameters to provide to prior
#'
#' @return A function
#' @export
#'
#' @examples
binom.prior <- function(type = c("MAP.FB", "MAP.EB","PP.FB","PP.EB","PP.EB.Beta","Beta.EB","PP.EB.Sep"),
                        x, n,
                        verbose=FALSE,
                        ...){

  type <- match.arg(type)

  f.to.call <- switch(type,
                      "MAP.FB" = binom.MAP.FB,
                      "MAP.EB" = binom.MAP.EB,
                      "PP.FB" = binom.PP.FB.MC,
                      "PP.EB" = binom.PP.EB,
                      "PP.EB.Beta" = binom.PP.EB.Beta,
                      "PP.EB.Sep" = binom.PP.EB.Sep,
                      "Beta.EB" = binom.Beta.EB
                      )

  do.call(f.to.call, args=list(x=x,n=n, verbose=verbose, ...))

}
