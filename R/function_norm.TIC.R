#' Normalise data to the TIC
#' @param imagedata.in Product of 'constructImage' function. 
#' @param offset number of columns that preceed image data
#' @export 

norm.TIC<-function(imagedata.in, offset){
  sums <- as.vector(colSums(imagedata.in[,(offset+1):ncol(imagedata.in)], na.rm = T))
  factor <- sums/mean(sums)
  image.norm <- cbind(imagedata.in[,1:offset], t(t(imagedata.in[,5:ncol(imagedata.in)])/factor))
  empty.spectra<-which(factor == 0, arr.ind=TRUE)
  if(length(empty.spectra)>0){
    for(i in 1:length(empty.spectra)){
      for(j in 1:nrow(image.norm)) image.norm[j,empty.spectra[i]+offset] <- 0
    }
  }
  image.norm<-cbind(imagedata.in[,1:offset], image.norm[,(offset+1):ncol(image.norm)])
#  hist(factor, breaks=100)
  return(image.norm)
}

