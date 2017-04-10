#' k-means clustering for imaging processing

#' @param cluster.type Currently only "kmeans" suported
#' @param imagedata.in Dataframe containing image data
#' @param offset columns preceding data
#' @param res.spatial spatial resolution of the image
#' @param width width of image; x.cood
#' @param height height of image; y.cood
#' @param clusters number of desired clusters
#' @return clustered images and cluster centers; writes csv files for cluster centers
#' @export 
#' 

################################################################################################################################

# define function

cluster <- function(cluster.type=cluster.type, imagedata.in = imagedata.in, offset = offset, res.spatial=res.spatial, width = x.cood, height = y.cood, clusters = clusters){
  
  
  # to do k-means clustering based on spectral similarity of pixels
  
  if (cluster.type == "kmeans"){
    k <- kmeans(t(imagedata.in[,(offset+1):ncol(imagedata.in)]), clusters)
    
    
    
    # rearrange your matrix to fit in image space, where nrow=y and ncol=x
    # transpose for the heatmap
    
    k.matrix <- data.frame()
    k.matrix <- matrix(k$cluster, nrow=height, ncol=width, byrow=T)
    t.k.matrix <- t(k.matrix)
    
    
    # plot the heatmap, choose colour scheme and colour according to cluster class
    
    
    filled.contour(
      x=seq(1, width, length=width)*res.spatial,
      y=seq(1, height, length=height)*res.spatial,
      z=t.k.matrix,
      col=grey(seq(0,1,length=10)),
      #color.palette=topo.colors,
      axes=TRUE,
      nlevels=12,
      #col=rainbow(10, alpha=0.5),
      asp=1, 
      plot.title = title(xlab= "x (micrometers)", ylab = "y (micrometers)", cex.axis = 0.5),
      key.title= title("Cluster")
    )
    title(main = paste("k-means clustering -", clusters, "clusters", sep=" "), cex.main =1)
    
    
    
    for (i in 1:clusters){
      colnames(k$centers) <- imagedata.in[,1]
      mz <- as.numeric(colnames(k$centers))
      Intensity <- k$centers[i,]
      plot(mz, Intensity, "h")
      abline(0, 0)
      title(main = paste("Cluster", i, sep=" "))
      
      one.cluster <- t.k.matrix
      one.cluster[one.cluster[,]!=i] <- 0
      
      filled.contour(
        x=seq(1, width, length=width)*res.spatial,
        y=seq(1, height, length=height)*res.spatial,
        z=one.cluster,
        # col=grey(seq(0,1,length=2)),
        col=c("black", sample(rainbow(10),1)),
        #col=c("black", "purple"),
        axes=TRUE,
        nlevels=2,
        asp=1, 
        plot.title = title(xlab= "x (micrometers)", ylab = "y (micrometers)", cex.axis = 0.5),
        
      )
      title(main = paste("k-means clustering -", "Cluster", i, sep=" "), cex.main =1)
      
      write.csv(Intensity, paste("Cluster_",i,".csv"))
    }
    
   
  } 
    
  
}


#################################################################################################
