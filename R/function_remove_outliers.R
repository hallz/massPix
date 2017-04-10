#' Remove outliers for image analysis
#' 
#' @param x image data in 
#' @param na.rm remove missing values
#' @param replace.1.min initial value to replace minimum values with
#' @param replace.1.max initial value to replace maximum values with
#' @return matrix of image data with outliers removed
#' @export
remove_outliers <- function(x, na.rm = TRUE,replace.1.min, replace.1.max) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- replace.1.min
  y[x > (qnt[2] + H)] <- replace.1.max
  
  y[x < (qnt[1] - H)] <- min(y)
  y[x > (qnt[2] + H)] <- max(y)
  y
}
