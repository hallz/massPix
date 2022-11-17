#' bin all m/zs by ppm bucket
#' 
#' @param extracted Product of mzextractor:list containing matrix of m/zs and intensities in the first element. 
#' @param bin.ppm Defines width of the ppm bin.
#' @return List containing two elements. The first is a matrix of extracted m/z and binned m/z for each extracted m/z, the second is a vector of all unique binned m/z values.
#' @export 

peakpicker.bin<-function(extracted, bin.ppm = 12){
  print("Starting peakpicker.bin...")
  spectra <- cbind(round(extracted[[1]],digits=4), "bin.med" = 0)
  marker <- 0
  i = 1
  ptm <- proc.time()
  while(i <= nrow(spectra)){
      ## group ions into common bin
      result <- spectra[ spectra[,1] >= spectra[i,1] & spectra[,1] <= (spectra[i,1]+(bin.ppm * spectra[i,1]) / 1000000),]
      if(class(result)[1]=="matrix" & spectra[i,"bin.med"] == 0){
        spectra[i:(i+nrow(result)-1),2]<-round(median(result[,1]),digits=4)
        i<-i+nrow(result)
      } 
      ## if single ion in bin, label with it's m/z
      else {
        spectra[i,2] <- spectra[i,1]
        i<-i+1
      }
      # progress report
      if(i/nrow(spectra)*100 > marker+1) {
        marker<- marker+1
        print(paste(marker, "% complete", sep=""))
      }
    
  }
  print(proc.time() - ptm)
  bin.spectra <- spectra
  finalmz <- unique(spectra[,2])
  summary <- paste(nrow(spectra), "ions across all pixels binned to",length(finalmz),"bins.", sep=" ")
  print(summary); log<-c(log,summary)
  peaks<-list(bin.spectra, finalmz)
  return(peaks)
}