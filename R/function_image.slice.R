#' Generate heat map based on spectral information.
#' 
#' @param row indices of row which corresponds to spectral data to plot; use image.norm_short.csv to hel identify row number of interest
#' @param imagedata.in Dataframe containing image data
#' @param scale range of scale that intensity values will be scaled to
#' @param x.cood width of image
#' @param y.cood height of image
#' @param nlevels Graduations of colour scale
#' @param name main name of image
#' @param subname sub name of image
#' @param res.spatial spatial resolution of image
#' @param offset number of columns preceding image data
#' @param rem.outliers Remove intensities that are outliers, valid arguments: "only" (without outliers only) or T (with and without outliers0)
#' @param summary T/F
#' @param title show titles, T/F"
#' @return heatmap image
#' @export

imageSlice<-function(row, imagedata.in, scale,x.cood,y.cood,nlevels,name, subname,offset,res.spatial, rem.outliers,summary, title=T){
  
  slice <- as.numeric(as.vector(as.matrix(imagedata.in[row,(1+offset):ncol(imagedata.in)])))
 
  
  if(rem.outliers != "only"){
    if(summary==F){
      boxplot(slice, main = (paste("distribution with 'outliers' for",imagedata.in[row,1],".",sep=" ")),cex.main = 1) 
      hist(slice, breaks<- seq(min(slice),max(slice),(max(slice)-min(slice))/100), prob=T, main=paste("distribution with 'outliers' for",imagedata.in[row,1],".",sep=" "),cex.main = 1)
      lines(density(slice),col="red",lwd=2) 
    }

    rescaled <-rescale(slice,scale)
    section <-t(matrix(rescaled, nrow=y.cood, ncol=x.cood, byrow=T))
    
    filled.contour(
      x=seq(from=1,to=x.cood,length=x.cood)*res.spatial, 
      y=seq(from=1,to=y.cood,length=y.cood)*res.spatial,
      z=section,
      nlevels=50,
      axes=TRUE,
      asp=1,
      plot.title = title(xlab = "x (micrometers)", ylab = "y (micrometers)", cex.axis= 0.5),
      color.palette=topo.colors
    )
    if(title == T){
      title(main = paste("with outliers", name, sep=" "), cex.main =1)
    }
    title(sub=subname, cex.sub = 1/2)
  }
  
  if(rem.outliers == T || rem.outliers == "only"){
    slice <- remove_outliers(x =slice, na.rm = TRUE,replace.1.min=0, replace.1.max=0) 
    if(summary==F){
      boxplot(slice, main = (paste("distribution without 'outliers' for",imagedata.in[row,1],".",sep=" ")),cex.main = 1) 
      hist(slice, breaks<- seq(min(slice),max(slice),(max(slice)-min(slice))/100), prob=T, main=paste("distribution without 'outliers' for",imagedata.in[row,1],".",sep=" "),cex.main = 1)
      lines(density(slice),col="red",lwd=2) 
    }
    
    rescaled <-rescale(slice,scale)
    section <-t(matrix(rescaled, nrow=y.cood, ncol=x.cood, byrow=T))
    
    filled.contour(
      x=seq(from=1,to=x.cood,length=x.cood)*res.spatial, 
      y=seq(from=1,to=y.cood,length=y.cood)*res.spatial,
      z=section,
      nlevels=50,
      axes=TRUE,
      asp=1,
      plot.title = title(xlab = "x (micrometers)", ylab = "y (micrometers)", cex.axis= 0.5),
      color.palette=topo.colors
    )
    if(title == T){
      title(main = paste("without outliers", name, sep=" "), cex.main = 1)
    }
    title(sub=subname, cex.sub = 1/2)
    
    
  }
}