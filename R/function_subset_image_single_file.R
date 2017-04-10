#' Generate subset of first image file to improve speed of deisotoping.
#' 
#' @param extracted Product of mzextractor:list containing matrix of m/zs and intensities in the first element. 
#' @param peaks Product of peakpicker.bin:list containing a vector of all unique binned m/z values in the 2nd element. 
#' @param percentage.deiso Defines the proportion of total pixels to select, at random from the first file to produce a subset of the image
#' @param thres.int Defines if intensity threshold, above which ions are retained
#' @param thres.low Defines the minumum m/z threshold, above which ions will be retained
#' @param thres.high Defines the minumum m/z threshold, below which ions will be retained
#' @param files a vector of file names
#' @param spectra_dir Defines the path to the spectral files
#' @param imzMLparse path to imzMLConverter
#' @return matrix of subset of first image file. variables are binned peak m/z, observations are pixels (n = percentage.deiso x size) chosen at random from the first image file.
#' @export 
#' 
subsetImage<-function(extracted,peaks,percentage.deiso=3, thres.int=10000, thres.low=200, thres.high=2000, files,spectra_dir,imzMLparse){
  print("Making image subset...")
  #R parser
  library(rJava)
  .jinit()
  .jaddClassPath(path=imzMLparse)  
  
  sizes<-extracted[[2]]
  bin.spectra <- peaks[[1]]
  finalmz <- peaks[[2]]
  file<-sample(1:length(sizes), 1, replace=F, prob=NULL)
#   imzML <- J("imzMLConverter.ImzMLHandler")$parseimzML(files[file]) 
#   imzML <- J("imzMLConverter.ImzMLHandler")$parseimzML(paste(getwd(),files[a], sep=""))
  imzML <- J("imzMLConverter.ImzMLHandler")$parseimzML(paste(getwd(),substr(files[file],2,100),sep=''))
  
  subset <- sample(1:as.numeric(sizes[[file]][3]), as.numeric(sizes[[file]][3])*(percentage.deiso/100), replace = FALSE, prob = NULL)
  temp.image<-cbind(as.numeric(finalmz),matrix(0, nrow=length(finalmz), ncol=length(subset)))
  marker <- 0
  for(n in 1:length(subset)){
   
    remainder <- subset[n]
    rows <- floor(remainder/as.numeric(sizes[[file]][2]))
    cols <- remainder-(rows*as.numeric(sizes[[file]][2]))
    if(cols==0){
      remainder<-remainder-1
      rows <- floor(remainder/as.numeric(sizes[[file]][2]))
      cols <- remainder-(rows*as.numeric(sizes[[file]][2]))
    }
    
    spectrum <- J(imzML, 'getSpectrum', as.integer(cols), as.integer(rows+1))     ## i is height, j is width
    mzs <- J(spectrum, 'getmzArray')
    counts <- J(spectrum, 'getIntensityArray')
    scan <- cbind("r.mz"=round(mzs, digits=4),counts)
    f.scan <- scan[scan[,2] > thres.int & scan[,1] > thres.low & scan[,1] < thres.high,,drop=FALSE]
    if(length(f.scan)>0){
      for(k in 1:nrow(f.scan)){
        bin.group <- bin.spectra[which(f.scan[k,1] == bin.spectra[,1], arr.ind=T ),2]
        if(length(bin.group)>0) temp.image[which(temp.image[,1] == bin.group),n+1] <- as.numeric(f.scan[k,"counts"])
      }
    }
    if((n/length(subset)*100) > marker+5){
      marker<- marker+5; print(paste(marker, "% done", sep=""))}
  }
  colnames(temp.image) <- c("mzbin", subset)
  return(temp.image)
}