#' Print messages if verbose=TRUE in calling function.
#'
#' @param message Message to print
#'
#' @export
#'
VP <- function(message){
  # if(parent.frame()$verbose)
    print(message)
}
