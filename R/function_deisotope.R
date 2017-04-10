#' Generate subset of first image file to improve speed of deisotoping.
#' 
#' @param ppm Tolerance (ppm) within which mass of isotope must be within . 
#' @param no_isotopes Number of isotopes to consider (1 or 2)
#' @param prop.1 Proportion of monoisotope intensity the 1st isotope intensity must not exceed 
#' @param prop.2 Proportion of monoisotope intensity the 2nd isotope intensity must not exceed 
#' @param peaks Product of peakpicker.bin:list containing a vector of all unique binned m/z values in the 2nd element. 
#' @param image.sub Product of subset.image: matrix containing intensity values to average for each binned m/z 
#' @param search.mod Search modifications T/F.
#' @param mod modifications to search eg. c(NL=T, label=F, oxidised=T,desat=T) 
#' @param lookup_mod A dataframe defining modifications
#' @return List of elements. Each element contains a dataframe with columns for m/z, mean intensity, deisotope annotation, modification annotation. Element 1, 2 and 3 have dataframes contain rows for all peaks, deisotoped and isotopes only. 
#' @export 

deisotope<-function(ppm=3, no_isotopes=2, prop.1=0.9,prop.2=0.5, peaks, image.sub, search.mod=F, mod=c(NL=T, label=F, oxidised=T,desat=T), lookup_mod){
  print("Starting deisotoping and difference scanning...")
  counts<-round(rowMeans(image.sub[,2:ncol(image.sub)]),digits=0)
  spectra <- cbind(peaks[[2]], counts,"", "")  
  colnames(spectra) <- c("mz.obs", "intensity", "isotope", "modification" )
  
C13_1 <- 1.003355
C13_2 <- C13_1*2

## set pmm window

## make column to store isotope annotation

## isotope counter
k = 0
m = 0

## run loop to find isotopes for each ion.
for(i in (1:nrow(spectra)-1)){
  ## values of search
  mass <- as.numeric(spectra[i,1])
  intensity <- as.numeric(spectra[i,2])
  ## calculated values    
  offset = (ppm * mass) / 1000000

  
  ## find isotope with ppm filter on isotpe
  search <- round((mass+C13_1),digits = 3)  
  top <- search + offset
  bottom <- search - offset
  result <- spectra[as.numeric(spectra[,"intensity"]) <= (intensity*prop.1) & spectra[,1] >= bottom & spectra[,1] <= top & spectra[,"isotope"] == "",]
  result <- rbind(result,blank1 = "", blank2 = "")

  if(no_isotopes ==2){
  ## find isotope with ppm filter on isotpe
    search <- round((mass+C13_2),digits = 3)  
    top <- search + offset
    bottom <- search - offset
    result_2 <- spectra[as.numeric(spectra[,"intensity"]) <= (intensity*prop.2) & spectra[,1] >= bottom & spectra[,1] <= top & spectra[,"isotope"] == "",]
    result_2 <- rbind(result_2,blank1 = "", blank2 = "")
    
  }

  result_3<-vector()
  if(search.mod != F){
    for(j in 1:nrow(lookup_mod)){
      if(mod[which(lookup_mod[j,"type"] == rownames(as.data.frame(mod)))] == T){
        ## find isotope with ppm filter on isotpe
        search <- round(mass+lookup_mod[j,"mass"],digits = 3)  
        top <- search + offset
        bottom <- search - offset
        if(0 != length(spectra[ spectra[,1] >= bottom & spectra[,1] <= top ,])){
          temp.hits <-rbind(spectra[ spectra[,1] >= bottom & spectra[,1] <= top ,])
          for(l in 1:nrow(temp.hits)){
            res_3_temp<-c(temp.hits[l,1], paste(as.vector(lookup_mod[j,"class"]),"(", row.names(lookup_mod[j,]),")", sep=""))
            result_3 <- rbind(result_3, res_3_temp)
          }
                      
         # res_3_temp <- c(spectra[spectra[,1] >= bottom & spectra[,1] <= top,], paste(as.vector(lookup_mod[j,"class"]), "(", row.names(lookup_mod[j,]),")", sep=""))
          
        }
      }
    }
  }
  result_3 <- rbind(result_3,blank1 = "", blank2 = "")
 
    
  # result<- as.matrix(result)
  ## add intensity filter
  # iso_intensity <- (((mass-380)/8.8)/100)*intensity
  
  if(nrow(result)>2){
    k = k +1
    spectra[i, "isotope"] <- paste(spectra[i, "isotope"]," ", "[",k,"]","[M]", sep="")
    for(j in 1:(nrow(result)-2)){
      indices <- which(spectra == result[j,1], arr.ind=TRUE)
      spectra[indices[,"row"], "isotope"] <- paste(spectra[indices[,"row"],"isotope"]," ","[",k,"]","[M+1]", sep="")
    }
    if(no_isotopes ==2 && nrow(result_2)>2){
      for(j in 1:(nrow(result_2)-2)){
        indices <- which(spectra == result_2[j,1], arr.ind=TRUE)
        spectra[indices[,"row"], "isotope"] <- paste(spectra[indices[,"row"],"isotope"]," ","[",k,"]","[M+2]", sep="")
      }
    }
  }   
  if(nrow(result_3)>2){
     m= m +1
    spectra[i, "modification"] <- paste(spectra[i, "modification"]," ","Precursor[",m,"]", sep="")
    for(j in 1:(nrow(result_3)-2)){
      indices <- which(spectra == result_3[j,1], arr.ind=TRUE)
      spectra[indices[,"row"], "modification"] <- paste(spectra[indices[,"row"],"modification"]," ","Fragment[",m,"] ",result_3[j,ncol(result_3)], sep="")
    }
  }
}  
allpeaks <- as.data.frame(spectra)
deisotoped  <- allpeaks[(grep("\\[M\\+", allpeaks$isotope, invert = T)),]
isotopes  <- allpeaks[(grep("\\[M\\+", allpeaks$isotope, invert = F)),]
results<-list(allpeaks, deisotoped, isotopes)
summary <- paste(length(as.vector(deisotoped$mz.obs)),"monoisotopic peaks retained and", length(as.vector(isotopes$mz.obs)), "C13 isotopes discarded from", length(as.vector(allpeaks$mz.obs)), "detected ions", sep=" ")
print(summary); log<-c(log,summary)
return(results)
}