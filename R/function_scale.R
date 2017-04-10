#' Normalise image data
#' 
#' @param imagedata.in Product of 'constructImage' function. 
#' @param scale.type Mode of scaling Valid argument: "c", "cs", "none" for centre and center + pareto scaling, respectively .
#' @param transform log transform data T/F
#' @param offset number of columns that precede image data
#' @return Centred, scaled, transformed dataframe containing image data.
#' @export 

centreScale<-function(imagedata.in=image.ann, scale.type="cs", transform=F, offset=4){
  
  # matrix of non-transformed data
  if(transform ==F){
    matr <-as.matrix(imagedata.in[,(offset+1):ncol(imagedata.in)])  
  }
  #log transform data to remove heteroscedasticity
  if(transform ==T){
    matr <-as.matrix(log(imagedata.in[,(offset+1):ncol(imagedata.in)] +1)) # log of 0 causes -inf values so +1..zeros (missing values) become 0
  }

  matr[matr == 0] <- NA # replace zeros for NA so they are ommited from center and scaling

  
  if(scale.type =="c"){
  matr <-as.matrix(log(imagedata.in[,(offset+1):ncol(imagedata.in)] +1)) # log of 0 causes -inf values so +1..zeros (missing values) become 0     }
  }
  if(scale.type =="cs"){
  matr<-t(scale(t(matr), center = TRUE, scale=apply(t(matr),2,function(x) sqrt(sd(x,na.rm=TRUE))))) #mean center and pareto scale
  }
  if(scale.type =="pareto"){
  matr <-t(scale(t(matr), center = FALSE, scale=apply(t(matr),2,function(x) sqrt(sd(x,na.rm=TRUE))))) #pareto scale only for slicing=T
  }
  if(scale.type =="none"){
    matr<-matr #no scaling
  }
  # replace NA with value (in this case 0)
  matr[which(is.na(matr))]=0;
  
  imagedata.in<-cbind(imagedata.in[,1:offset], matr)
  return(imagedata.in)
}