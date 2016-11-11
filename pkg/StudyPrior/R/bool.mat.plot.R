#' Plot a boolean matrix in two colours
#'
#' @param mat Matrix with TRUE or FALSE entries
#' @param col.pair Vector of colours. First entry is used for FALSE and second for TRUE
#' @param ... extra parameters to plot
#'
#' @return Nothing
#' @export
#'

bool.mat.plot <- function(mat, col.pair=c('darkred','lightblue'),...){
  x <- ncol(mat)
  y <- nrow(mat)

  plot(x=c(0,x), y=c(0,y), type='n', yaxs='i', xaxs='i',...)


  C <- rep(seq.int(x), times = y)
  R <- rep(seq.int(y), each = x)
  cols <- mapply(function(x,y) mat[x,y], R,C)


  rect((C-1), (R-1),
       C, R,
       col = col.pair[1+cols],
       border=TRUE,
       lty=2,
       ...)

  box()
  #
  # sapply(1:x, function(C){
  #   sapply(1:y, function(R){
  #     rect((C-1), (R-1),
  #          C, R,
  #          col = ifelse(mat[R,C],col.pair[2],col.pair[1]),
  #          border=TRUE,
  #          lty=2,
  #          ...)
  #   })
  # })
  #  box()
}

