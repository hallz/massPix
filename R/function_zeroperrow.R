#' determine the number of zeros per row in a matrix
#' @param steps Sequence of values between 0 and 1 that define the thresholds of missing values to test
#' @param matrix Matrix of image data
#' @return List of row indices corresponing to variables that were retained at each thres.filter step.
#' @export 

zeroperrow <- function(steps, matrix){
## filter summary for MVA
indices <-list()
results <- vector()
counter = 1
for(threshold in steps){
  filter<-vector()
  for(i in 1:nrow(matrix)){
    if((length(which(matrix[i,]<1))/ncol(matrix))<=threshold) filter<-c(filter,i)
  }
 print(paste("Step ", counter,": ", length(filter)," records at ",threshold," threshold." ,sep=""))
## record the retained indices at each threshold  
  indices[counter] <- list(filter)
  results<-c(results, length(matrix[filter,1]))
  counter<-counter+1
}
plot(steps, results[1:(length(results))], ylab="Number of ions", xlab="% missing values", main="Coverage of ions across selected pixels",type="s")
results
return(indices)
}

