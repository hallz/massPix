#' extract all m/zs from image data
#' 
#' @param files spectra raw file (consisting of .imzML and .ibd file); multiple files processing in devleopment, at the moment one file at a time
#' @param spectra_dir Defines the path to the spectral files
#' @param imzMLparse path to imzMLConverter
#' @param thres.int Defines if intensity threshold, above which ions are retained
#' @param thres.low Defines the minumum m/z threshold, above which ions will be retained
#' @param thres.high Defines the minumum m/z threshold, below which ions will be retained
#' @return List containing two elements. The first is a numeric vector containing of all unique masses, the second a list of height, width and height.width in pixels for each file.
#' @export 

mzextractor<-function(files,spectra_dir, imzMLparse, thres.int= 10000, thres.low=200, thres.high=1000){
  print("Starting mzextractor...")
  ## load parse and java
  library(rJava)
  .jinit()
  .jaddClassPath(path=imzMLparse)
  
  
  log<-vector()
  sizes<-list()
  all.mzs <- vector()
  for(a in 1:length(files)){
    imzML <- J("imzMLConverter.ImzMLHandler")$parseimzML(paste(spectra_dir,substr(files[a],2,100),sep=''))
    width <- J(imzML, 'getWidth')
    height <- J(imzML, 'getHeight')
    size<-height*width; sizes[[a]]<-list(height, width, size)
    
    # plot(mzs, counts, "l")
    ## determine list of image wide m/z values -> store in vectore s.mz (unique and sorted)
    marker <- 0
    for(i in 1:(height*1)){
      for(j in 1:(width*1)){
        spectrum <- J(imzML, 'getSpectrum', as.integer(j), as.integer(i)) #i is width, j is height  
        mzs <- J(spectrum, 'getmzArray')
        counts <- J(spectrum, 'getIntensityArray')
        scan <- cbind("r.mz"=round(mzs, digits=4),counts)
        f.scan <- scan[scan[,2] > thres.int,,drop=FALSE]
        all.mzs<-rbind(all.mzs, f.scan[,1:2])
        #all.mzs<-c(all.mzs, f.scan[,1])
        
      }
      # progress report
      if((i/height)*100 > marker+5){
        marker<- marker+5; print(paste(marker, "% complete", sep=""))}
    }
    print(paste("mzs from", size, "spectra in", files[a], "now read.", sep=" "))
  } ## end of reading in files
  
  final.size = 0; for(a in 1:length(sizes)) final.size <- final.size+ as.numeric(sizes[[a]][3])
  mz <- all.mzs[,1]
  u.mz<- unique(mz)
  s.mz<- sort(u.mz)
  f.s.mz <- s.mz[s.mz> thres.low & s.mz < thres.high ]
  summary <- paste(length(f.s.mz),"unique ions retained within between",thres.low, "m/z and", thres.high, "m/z from",length(u.mz),"unique ions detected across all pixels.", sep=" ")
  print(summary); log<-c(log,summary)
  output <-list(f.s.mz, sizes)
  return(output)
  
}