#' construct dataframe of image data
#' 
#' @param extracted Product of mzextractor: list containing a list of the image diminsions in the 2nd element. 
#' @param deisotoped Product of deisotope
#' @param peaks Product of peakpicker.bin:list containing a matrix of all peaks and the corresponding m/z bin in the first element. 
#' @param imzMLparse path to imzMLConverter
#' @param spectra_dir Defines the path to the spectral files
#' @param thres.int Defines if intensity threshold, above which ions are retained
#' @param thres.low Defines the minumum m/z threshold, above which ions will be retained
#' @param thres.high Defines the minumum m/z threshold, below which ions will be retained
#' @param files a vector of file names
#' @return Dataframe of image data. Variables are deisotoped binned m/z, observations are pixels. For image of width w and height h, the number of columns is w x h. The first w columns are from the first row (from left to right), the next w columns are the next row, from left to right and so on. List, first element is has column containing m/z values preceeding image data, second element has 4 columns preceeding image data which include m/z, annotation, isotope status, modification status. 
#' @export 


contructImage<-function(extracted,deisotoped,peaks,imzMLparse, spectra_dir,thres.int=10000, thres.low=200, thres.high=1000,files){
  print("Starting 'constructImage'...")
  library(rJava)
  .jinit()
  .jaddClassPath(path=imzMLparse)  
  
  sizes<-extracted[[2]]
  final.size = 0; for(a in 1:length(sizes)) final.size <- final.size+ as.numeric(sizes[[a]][3])
  bin.spectra <- peaks[[1]]
  d.finalmz <-as.vector(deisotoped[[2]]$mz.obs)
  image<-cbind(as.numeric(d.finalmz),matrix(0, nrow=length(d.finalmz), ncol=final.size))
  
  for(a in 1:length(files)){
    height<-as.numeric(sizes[[a]][1])
    width<-as.numeric(sizes[[a]][2])
    imzML <- J("imzMLConverter.ImzMLHandler")$parseimzML(paste(spectra_dir,substr(files[a],2,100),sep=''))
    marker <- 0
    for(i in 1:height){
      for(j in 1:width){
        spectrum <- J(imzML, 'getSpectrum', as.integer(j), as.integer(i))    
        mzs <- J(spectrum, 'getmzArray')
        counts <- J(spectrum, 'getIntensityArray')
        scan <- cbind("r.mz"=round(mzs, digits=4),counts)
        f.scan <- scan[scan[,2] > thres.int & scan[,1] > thres.low & scan[,1] < thres.high,,drop=FALSE]
        if(length(f.scan)>0){
          for(k in 1:nrow(f.scan)){
            bin.group <- bin.spectra[which(f.scan[k,1] == bin.spectra[,1], arr.ind=T ),2]
            if(length(bin.group)>0) image[which(image[,1] == bin.group),((i-1)*width)+(j+1)] <- as.numeric(f.scan[k,"counts"])
            #if(length(result)>0) image[result,((i-1)*width)+(j+1)]<- agg.mzs[k,"counts"]
          }
        }
      }
      # progress report
      
      if((i/height)*100 > marker+1){
        marker<- marker+1; print(paste(marker, "% done", sep=""))}
    }
  }; summary<- "imzMLcube constructed!"; log<-c(summary)
  
  print(log)
  return(image)
}
