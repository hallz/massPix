#' Normalise image data
#' 
#' @param imagedata.in Product of 'constructImage' function. 
#' @param norm.type Mode of normalisation. Valid argument: "standards", "TIC", "median", "none".
#' @param standards vector of row indices corresponding to variables that are standards.
#' @param offset number of columns that preceed image data
#' @return Normalised dataframe containing image data.
#' @export 

normalise<-function(imagedata.in=image.ann, norm.type="TIC", standards=NULL, offset=4){
  if(norm.type == "standards"){
  ## from standards
  images.f.n<-norm.standards(imagedata.in,offset,standards)
  imagedata.in <- images.f.n
  }
  if(norm.type == "TIC"){
  ## from TIC
  images.f.n <- norm.TIC(imagedata.in,offset)
  image
  imagedata.in <- images.f.n
  }
  if(norm.type == "median"){
  ## from median
  images.f.n <- norm.median(imagedata.in,offset)
  imagedata.in <- images.f.n
  }
  if(norm.type == "none"){
  ## no normalisation
  imagedata.in<-imagedata.in
  }
  return(imagedata.in)
}