#' Rescale image values
#' 
#' @param x image date in 
#' @param scale range of scaled data
#' @return matrix of scaled data
#' @export

rescale <- function(x,scale){((x-min(x))/(max(x)-min(x)))*scale}
