#' make library of lipid masses
#' 
#' @param ionisation_mode Choose "positive" or "negative", will determine which lipid classes are in database
#' @param sel.class A vector defining classes of lipids to be included in the library
#' @param fixed Defines if one of the SN positions is fixed, default is F
#' @param fixed_FA Defines the name of the fixed FA eg 16, 18.1, 20.4.
#' @param lookup_lipid_class A data.frame defining lipid classes frm library
#' @param lookup_FA A dataframe defining FAs from library
#' @param lookup_element A data.frame defining elements from library
#' @return Dataframe of masses for all combinations of FAs with chosen head groups
#' @export 

makelibrary <-function(ionisation_mode, sel.class, fixed=F, fixed_FA, lookup_lipid_class,lookup_FA, lookup_element){
  print("Making library of lipid masses...")
  
  
  if (ionisation_mode == "positive"){
    sel.class<- c( 
      T, #TG
      T, #DG
      T, #PC
      F, #PA
      T, #PE
      T, #PS
      F, #PG
      F, #PI
      F, #PIP
      F, #PIP2
      F, #PIP3
      T, #LysoPC
      T, #DG-H20
      T, #CE
      F, #FFA
      T, #SM
      T  #Cer
    )
  } 
  
  if (ionisation_mode =="negative"){
    sel.class<- c( 
      F, #TG
      F, #DG
      T, #PC
      T, #PA
      T, #PE
      T, #PS
      T, #PG
      T, #PI
      T, #PIP
      T, #PIP2
      T, #PIP3
      F, #LysoPC
      F, #DG-H20
      F, #CE
      T, #FFA
      F, #SM
      F  #Cer
    )
  }
  
  lookup_lipid_class <- cbind(lookup_lipid_class, sel.class)
  
  # FAs to use in library
  FA_expt <-list('10','12','14','15','16','16.1','17','17.1','18','18.1','18.2','18.3','20.3','20.4','20.5','21','22','22.5','22.6','23','24.1')
  
  library<- numeric()
  for(i in 1:nrow(lookup_lipid_class)){
    if(lookup_lipid_class[i,"sel.class"] == T){
      ## key variables
      rounder = 3 # number of decimals the rounded masses are rounded to.
      ## lipidclass = "TG"
      lipidclass <- row.names(lookup_lipid_class[i,])
      
      
      ##determine how many FAs places to be used for combination and generate combination of FAs
      FA_number <- as.numeric(lookup_lipid_class[lipidclass,"FA_number"])
      if(fixed==TRUE)FAnum <- FA_number - 1 else FAnum <- FA_number                 
      s1<-combn(FA_expt,FAnum)
      
      ## if one place is fixed add this FA to the matrix
      if(fixed==TRUE){
        s1 <- rbind(s1,"fixed"=fixed_FA)
        FAnum <- FAnum +1
      }
      
      ## if sn2 or sn3 does not have FA bind 'empty' FA channel.
      if(FAnum == 1){ s1<-rbind(s1,sn2<-vector(mode="numeric",length=ncol(s1)),sn3<-vector(mode="numeric",length=ncol(s1))) 
                      FAnum = FAnum +2}
      if(FAnum == 2){ s1<-rbind(s1,sn3<-vector(mode="numeric",length=ncol(s1))) 
                      FAnum = FAnum +1}
      
      
      ## label the matrix
      if(FAnum == 3) row.names(s1) <-c("FA1", "FA2","FA3")
      
      ## add rows to matrix for massofFAs and formula
      massofFAs<-vector(mode="numeric",length=ncol(s1))
      s1 <- rbind(s1,massofFAs)
      formula<-vector(mode="numeric",length=ncol(s1))
      s1 <- rbind(s1,formula)
      ##row.names(s1) <-c("FA1", "FA2","FA3", "massofFAs")
      for(i in 1:ncol(s1)){
        
        ## for 3 FAs
        if(FAnum == 3){                     
          FA_1 <- as.character((s1[1,i]))  
          FA_2 <- as.character((s1[2,i]))  
          FA_3 <- as.character((s1[3,i])) 
          s1["massofFAs",i] <- as.numeric((lookup_FA[FA_1, "FAmass"]))+as.numeric((lookup_FA[FA_2, "FAmass"]))+as.numeric((lookup_FA[FA_3, "FAmass"]))
          ##determine the formula
          temp_carbon <- as.numeric((lookup_FA[FA_1, "FAcarbon"]))+as.numeric((lookup_FA[FA_2, "FAcarbon"]))+as.numeric((lookup_FA[FA_3, "FAcarbon"]))
          temp_doublebond <- as.numeric((lookup_FA[FA_1, "FAdoublebond"]))+as.numeric((lookup_FA[FA_2, "FAdoublebond"]))+as.numeric((lookup_FA[FA_3, "FAdoublebond"]))
          s1["formula",i] <- paste(lipidclass,"(",temp_carbon,":",temp_doublebond,")", sep = "")
        }
      }
      
      ## calculate total mass
      totalmass<-vector(mode="numeric",length=ncol(s1))
      s1 <- rbind(s1,totalmass)
      
      for(i in 1:ncol(s1)){
        s1["totalmass",i] <- as.numeric(s1["massofFAs",i]) + as.numeric(as.character(lookup_lipid_class[lipidclass,"headgroup_mass"])) - (as.numeric(lookup_lipid_class[lipidclass,"FA_number"])*as.numeric(lookup_element["H","mass"]))   
      }
      
      ##make rows for charged lipids masses 
      protonated<-vector(mode="numeric",length=ncol(s1))
      ammoniated<-vector(mode="numeric",length=ncol(s1))
      sodiated<-vector(mode="numeric",length=ncol(s1))
      potassiated<-vector(mode="numeric",length=ncol(s1))
      deprotonated<-vector(mode="numeric",length=ncol(s1))
      chlorinated<-vector(mode="numeric",length=ncol(s1))
      acetate<-vector(mode="numeric",length=ncol(s1))  
      s1 <- rbind(s1,protonated, ammoniated, sodiated, potassiated, deprotonated, chlorinated, acetate)
      
      ##calculate charged lipids masses 
      for(i in 1:ncol(s1)){
        s1["protonated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["H","mass"])),digits = 4)  
        s1["ammoniated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["NH4","mass"])),digits = 4) 
        s1["sodiated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["Na","mass"])),digits = 4)  
        s1["potassiated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["K","mass"])),digits = 4)  
        s1["deprotonated",i] <- round((as.numeric(s1["totalmass",i]) - as.numeric(lookup_element["H","mass"])),digits = 4)  
        s1["chlorinated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["Cl","mass"])),digits = 4)  
        s1["acetate",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["CH3COO","mass"])),digits = 4)  
      }
      
      ##make rows for rounded charged lipids masses
      round.protonated<-vector(mode="numeric",length=ncol(s1))
      round.ammoniated<-vector(mode="numeric",length=ncol(s1))
      round.sodiated<-vector(mode="numeric",length=ncol(s1))
      round.potassiated<-vector(mode="numeric",length=ncol(s1))
      round.deprotonated<-vector(mode="numeric",length=ncol(s1))
      round.chlorinated<-vector(mode="numeric",length=ncol(s1))
      round.acetate<-vector(mode="numeric",length=ncol(s1))
      s1 <- rbind(s1,round.protonated, round.ammoniated, round.sodiated, round.potassiated, round.deprotonated, round.chlorinated, round.acetate)
      
      ##calculate rounded charged lipids masses 
      for(i in 1:ncol(s1)){
        s1["round.protonated",i] <- round(as.numeric(s1["protonated",i]),digits = rounder)
        s1["round.ammoniated",i] <- round(as.numeric(s1["ammoniated",i]),digits = rounder)
        s1["round.sodiated",i] <- round(as.numeric(s1["sodiated",i]),digits = rounder)
        s1["round.potassiated",i] <- round(as.numeric(s1["potassiated",i]),digits = rounder) 
        s1["round.deprotonated",i] <- round(as.numeric(s1["deprotonated",i]),digits = rounder) 
        s1["round.chlorinated",i] <- round(as.numeric(s1["chlorinated",i]),digits = rounder) 
        s1["round.acetate",i] <- round(as.numeric(s1["acetate",i]),digits = rounder) 
      }
      
      library<-cbind(library,s1)
    }
  }
  return(library)
}