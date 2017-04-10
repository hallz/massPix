#' Normalise data to standards
#' @param imagedata.in Product of 'constructImage' function. 
#' @param offset number of columns that preceed image data
#' @param standards vector of row indices corresponding to variables that are standards.
#' @export 

norm.standards<-function(imagedata.in,offset,standards){
for(i in 1:length(standards)) boxplot(as.vector(t(imagedata.in[standards[i],(offset+1):length(imagedata.in)])), main = (paste("distribution of standard",imagedata.in[standards[i],1],".",sep=" "))) 

av <- mean(colMeans(imagedata.in[standards,(offset+1):length(imagedata.in)]))  
norm <- (t(imagedata.in[,(offset+1):ncol(imagedata.in)])/colMeans(imagedata.in[standards,(offset+1):length(imagedata.in)]))
norm <- norm*av
norm.rep <- replace(norm, norm == "NaN" , 0)
norm.rep <- replace(norm.rep, norm == "Inf", 0)
image.norm <- cbind(imagedata.in[,1:offset], t(norm.rep))
image.norm 
}


