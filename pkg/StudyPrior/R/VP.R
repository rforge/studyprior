#' Print messages if verbose=TRUE in calling function.
#'
#' @param message
#'
#' @return
#' @export
#'
#' @examples
VP <- function(message){
  # if(parent.frame()$verbose)
    print(message)
}
