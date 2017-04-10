#' Filter to a matrix subset that includes variables above a threshold of missing values
#' 
#' @param imagedata.in dataframe containing spectral informaion
#' @param thres.filter Defines threshold for proportion of missing values (this is the step number, not the actual proportion)
#' @param offset Defines the number of columns that preceed intensitiy values
#' @param steps Sequence of values between 0 and 1 that define the thresholds of missing values to test
#' @return Matrix of image data containing only variables with a missing value proporiton below thres.filter.
#' @export 
 
filter<-function(imagedata.in, steps=seq(0,1,0.05), thres.filter=11,offset=4){
  print("Starting filter")
  zeros <- zeroperrow(steps, as.matrix(imagedata.in[,(offset+1):ncol(imagedata.in)]))
  image.filtered<-imagedata.in[zeros[[thres.filter]],]
  return(image.filtered)
}
