#' Create heat map for PCA image
#' 
#' Generate heat map based on spectral information.
#' 
#' @param imagedata.in Dataframe containing image data.
#' @param offset number of columns preceeding image data
#' @param PCnum number of PCs to consider
#' @param scale range of scale that intensity values will be scaled to.
#' @param x.cood width of image.
#' @param y.cood height of image.
#' @param nlevels Graduations of colour scale.
#' @param res.spatial spatial resolution of image
#' @param rem.outliers Remove intensities that are outliers, valid arguments: "only", "true".
#' @param summary T/F 
#' @param title show titles, T/F".
#' @return heatmap image
#' @export

imagePca<-function(imagedata.in, offset, PCnum, scale,x.cood,y.cood,nlevels,res.spatial,summary, title=T, rem.outliers="only"){
  pca <-princomp(t(imagedata.in[,(offset+1):ncol(imagedata.in)]), cor = FALSE, scores = TRUE, covmat = NULL)
  for(i in 1:PCnum){
    imageSlice(row=i, imagedata.in=t(pca$scores), scale, x.cood, y.cood, nlevels,name=paste("PC",i,sep=""),subname="", offset=0, res.spatial, rem.outliers, summary,title)
 
    
    library(calibrate)
    labs.all<-as.numeric(as.vector(imagedata.in[,1]))
    perc<-5 # percent of data to label
    y<-cut(pca$loadings[,i],breaks=c(-Inf,quantile(pca$loadings[,i],p=c(perc/100)),quantile(pca$loadings[,i],p=c(1-(perc/100))),Inf),labels=c("low","mid", "long"))
    labs<-labs.all; labs[which(y == "mid")]<-""
    plot(x =as.numeric(as.vector(imagedata.in[,1])), y =pca$loadings[,i],type="n", main=paste("PC",i, " loadings",sep=""), xlab = "m/z", ylab="p[4]")
    lines(x =as.numeric(as.vector(imagedata.in[,1])), y =pca$loadings[,i], type="h") 
    textxy(X =as.numeric(as.vector(imagedata.in[,1])), Y =pca$loadings[,i], labs=labs, cx = 0.5, dcol = "black", m = c(0, 0))
    
  }
}