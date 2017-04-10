#' Performs identification of lipids and adducts
#' 
#' @param ionisation_mode "positive" or "negative" determines which adducts to search for
#' @param deisotoped Product of function 'deisotope': list with dataframe in the second element containing image data for deisotoped data
#' @param adducts vector of adducts to be searched in the library of lipid masses
#' @param ppm.annotate Defines if ppm threshold for which |observed m/z - theotical m/z| must be less than for annotation to be retained
#' @param dbase The product of library - dataframe containing lipid massess.
#' @return Character vector containing annotations
#' @export 


annotate<-function(ionisation_mode, deisotoped, adducts=c(H = T, NH4 = F, Na = T, K = T, dH = F, Cl = F, OAc = F),ppm.annotate=10, dbase){
  print("Starting annotation")
d.finalmz <-as.vector(deisotoped[[2]]$mz.obs)
s1<-dbase
spectra<-cbind(round(as.numeric(d.finalmz), digits=3),d.finalmz)
combined <- vector()
sel.adducts<- vector()
index<- 13 # offset to search only rounded masses in library
if (ionisation_mode=="positive"){
  adducts<- c(H = T, NH4 = F, Na = T, K = T, dH = F, Cl = F, OAc = F)
}
if (ionisation_mode=="negative"){
  adducts<- c(H = F, NH4 = F, Na = F, K = F, dH = T, Cl = T, OAc = F)
}
for(a in 1:length(adducts)){
  if(adducts[a]==T) sel.adducts<- c(sel.adducts, index+a)
}
for(i in 1:nrow(spectra)){
  search<-as.numeric(spectra[i,1])
  offset = (ppm.annotate * search) / 1000000
  top <- search + offset
  bottom <- search - offset
  result<-which( s1[sel.adducts,] >= bottom & s1[sel.adducts,] <= top ,arr.ind=TRUE)
  if(nrow(result)>0){
    for(j in 1:nrow(result))
    {
      col<-result[j,"col"]
      row<-result[j,"row"]; row<-sel.adducts[row]
      ## determine the adduct that was matched, summarising match information from library for matched mass (as 'data')
      ##determine which adduct
      if(row == "14"){adduct <- "protonated"; name.adduct <- "H"}
      if(row == "15"){adduct <- "ammoniated"; name.adduct <- "NH4"}
      if(row == "16"){adduct <- "sodiated"; name.adduct <- "Na"}
      if(row == "17"){adduct <- "potassiated"; name.adduct <- "K"}
      if(row == "18"){adduct <- "deprotonated"; name.adduct <- "-H"}
      if(row == "19"){adduct <- "chlorinated"; name.adduct <- "Cl"}
      if(row == "20"){adduct <- "acetate"; name.adduct <- "OAc"}
      
      a.ppm = round(abs(((as.numeric(spectra[i,2]) - as.numeric(s1[adduct,col])) / as.numeric(spectra[i,2]))* 1000000), digits=1)
      
      ##make vector with summary of match and paired match 
      data<- c(s1[row,col],  s1[adduct,col],spectra[i,2], a.ppm, s1["formula",col], name.adduct,  s1["protonated",col],  s1["FA1",col], s1["FA2",col], s1["FA3",col])
      
      ##make matrix of search results
      combined <- rbind(combined,unlist(data, use.names=F))
    }
  }
}
if(length(combined)>0){
  colnames(combined) <- c("mz.matched",  "mz.matched.lib", "mz.observed", "ppm", "formula", "adduct", "mz.lib.protonated", "FA1", "FA2", "FA3")
  
  ids <- unique.matrix(combined[,c(3, 5,6)])
  annotations<-cbind(d.finalmz, "")
  for(i in 1:nrow(annotations)){
    result <- which(ids[,1] == annotations[i,1], arr.ind = T)
    if(length(result)>0){
      for(j in 1:length(result)){
        annotations[i,2]<-paste(annotations[i,2],"[",ids[result[j],"formula"],"+",ids[result[j],"adduct"],"]", sep="") 
      }       
    } 
  }

#image.roi<-cbind(image.roi[,2:ncol(image.roi)])
#vresults<-cbind(annotations[,2], image.roi[,1],image.roi[,2:ncol(image.roi)])
#colnames(image)<-c("annoation","modification","bin.mz", files)

summary <- paste(length(annotations[annotations[,2] !="",2]),"from",length(as.vector(deisotoped[[2]]$mz.obs)),"monoisotopic peaks were annoated (using accuract mass) with a",ppm.annotate,"ppm tollerance", sep=" "); log<-c(log,summary)
print(summary)

return(annotations[,2])
}
if(length(combined)==0) print("No annotations were made")
}