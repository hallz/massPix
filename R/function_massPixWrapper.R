#' wrapper to process image data
#' 
#'
#' @param process T to process data from imzML, F to skip processing
#' @param pca T to conduct PCA image analysis
#' @param slice T to generate two dimensional image of one ion defined in row
#' @param cluster.k T to perform clustering
#' @param ionisation_mode Choose "positive" or "negative"; determines which lipid classes to use when building database
#' @param thres.int Defines if intensity threshold, above which ions are retained
#' @param thres.low Defines the minumum m/z threshold, above which ions will be retained
#' @param thres.high Defines the minumum m/z threshold, below which ions will be retained
#' @param bin.ppm Mass accuracy for binning between spectra
#' @param thres.filter Defines threshold for proportion of missing values (this is the step number, not the actual proportion)
#' @param ppm.annotate Defines if ppm threshold for which |observed m/z - theotical m/z| must be less than for annotation to be retained
#' @param res.spatial spatial resolution of image
#' @param PCnum number of PCs to calculate
#' @param row indices of row which corresponds to spectral data to plot in image slice - use image.norm_short.csv to find row number
#' @param cluster.type At the moment only supports "kmeans"
#' @param clusters Number of clusters to use for clustering
#' @param ppm Tolerance (ppm) within which mass of isotope must be within . 
#' @param no_isotopes Number of isotopes to consider (1 or 2)
#' @param prop.1 Proportion of monoisotope intensity the 1st isotope intensity must not exceed 
#' @param prop.2 Proportion of monoisotope intensity the 2nd isotope intensity must not exceed 
#' @param search.mod Search modifications T/F.
#' @param mod modifications to search eg. c(NL=T, label=F, oxidised=T,desat=T) 
#' @param lookup_mod A dataframe defining modifications
#' @param adducts vector of adducts to be searched in the library of lipid masses
#' @param sel.class A vector defining classes of lipids to be included in the library
#' @param fixed Defines if one of the SN positions is fixed, default is F.
#' @param fixed_FA Defines the name of the fixed FA is fixed is T, e.g. 16, 16.1, 18.2.
#' @param lookup_lipid_class A data.frame defining lipid classes
#' @param lookup_FA A dataframe defining FAs
#' @param lookup_element A dataframe defining elements
#' @param files a vector of file names; currently only supports processing one file at a time
#' @param spectra_dir Defines the path to the spectral files, 
#' @param imzMLparse path to imzMLConverter
#' @param percentage.deiso Defines the proportion of total pixels to select, at random from the first file to produce a subset of the image
#' @param steps Sequence of values between 0 and 1 that define the thresholds of missing values to test
#' @param imagedata.in Dataframe containing image
#' @param norm.type "TIC" or "median" or "standards" or "none", recommend "TIC"
#' @param standards default is NULL, vector of row indices corresponding to variables that are standards
#' @param scale.type options are "cs" center and pareto scale, "c" center only, "pareto" pareto only and "none"
#' @param scale range of scale that intensity values will be scaled to
#' @param transform T/F to perform log transform
#' @param x.cood width of image
#' @param y.cood height of image
#' @param nlevels Graduations of colour scale.
#' @param name main name of image.
#' @param subname sub name of image
#' @param rem.outliers Choose "only" to show only images after removing outliers, T for both with and without outliers
#' @param summary T/F to show extra information
#' @param title show titles, T/F"
#' @param offset number of columns of text preceding data
#' 
#' @return Dataframe of annotated image data. Variables are deisotoped binned m/z, observations are pixels. For image of width w and height h, the number of columns is w x h. The first w columns are from the first row (from left to right), the next w columns are the next row, from left to right and so on.  

#' @export 



massPix <- function(
  process=T,
  pca=T,
  slice=T,
  cluster.k=T,
  ionisation_mode = "positive",
  thres.int = 500000,
  thres.low = 200,
  thres.high = 1000,
  bin.ppm = 12,
  thres.filter = 11,
  ppm.annotate = 12,
  res.spatial = 50,
  PCnum = 5,
  row = 50,
  cluster.type = "kmeans",
  clusters = 6,
  ppm = 3,
  no_isotopes = 2,
  prop.1 = 0.9,
  prop.2 = 0.5,
  search.mod = T,
  mod = c(NL = T, label = F, oxidised = T, desat = T),
  lookup_mod,
  adducts,
  sel.class,
  fixed = F,
  fixed_FA,
  lookup_lipid_class,
  lookup_FA,
  lookup_element,
  files,
  spectra_dir,
  imzMLparse, 
  percentage.deiso = 3,
  steps = seq(0, 1, 0.05),
  imagedata.in=imagedata.in,
  norm.type = "TIC",
  standards = NULL,
  scale.type = "cs",
  transform = F,
  scale = 100,
  x.cood = x.cood,
  y.cood = y.cood,
  nlevels = 50,
  name = "",
  subname = "",
  rem.outliers = "only",
  summary = T,
  title = T,
  offset=4){

######################################################################################################################################




  ### setting up, reading library files


  ## to read and process very large files

  options(java.parameters = "Xmx4g")

  library(rJava)

  
  ## read in library files

  setwd(lib_dir)
  read <- read.csv('lib_FA.csv', sep=",", header=T);lookup_FA <- read[,2:4]; row.names(lookup_FA) <- read[,1]
  read <- read.csv('lib_class.csv', sep=",", header=T);lookup_lipid_class <- read[,2:3]; row.names(lookup_lipid_class) <- read[,1]
  read <- read.csv('lib_element.csv', sep=",", header=T);lookup_element <- read[,2:3]; row.names(lookup_element) <- read[,1]
  read <- read.csv('lib_modification.csv', sep=",", header=T);lookup_mod <- read[,2:ncol(read)]; row.names(lookup_mod ) <- read[,1]


  ### parsing the data and getting x and y dimensions

  imzMLparse<-  paste(home_dir,"imzMLConverter/imzMLConverter.jar",sep="") # spectra files

  setwd(spectra_dir)
  files <- list.files(, recursive = TRUE, full.names = TRUE,pattern = ".imzML")

  setwd(spectra_dir)
  .jinit()
  .jaddClassPath(path=imzMLparse)
  imzML <- J("imzMLConverter.ImzMLHandler")$parseimzML(paste(getwd(),substr(files,2,100),sep=''))
  x.cood <- J(imzML, 'getWidth')
  y.cood <- J(imzML, 'getHeight')



###################################################################################################################################

###main script

  # make library
if (process==T){  
  dbase <-makelibrary(ionisation_mode, sel.class, fixed, fixed_FA, lookup_lipid_class,lookup_FA, lookup_element) 
  print(paste("library containing", ncol(dbase),"entries was made",sep=" "))

  # extract m/z and pick peaks
  setwd(spectra_dir)
  extracted<-mzextractor(files, spectra_dir, imzMLparse,thres.int, thres.low, thres.high)
  peaks <- peakpicker.bin(extracted, bin.ppm) #pick peaks
  
  # perform deisotoping on a subset of the image
  temp.image <-  subsetImage(extracted, peaks, percentage.deiso,thres.int, thres.low, thres.high, files, spectra_dir,imzMLparse)
  temp.image.filtered <- filter(imagedata.in=temp.image, steps, thres.filter, offset = 1)
  deisotoped<-deisotope(ppm, no_isotopes, prop.1, prop.2, peaks=list("",temp.image.filtered[,1]), image.sub=temp.image.filtered, search.mod, mod, lookup_mod)
  
  # perform annotation of lipids using library
  annotated<-annotate(ionisation_mode, deisotoped, adducts,ppm.annotate, dbase)

  # make full image and add lipid ids
  final.image<-contructImage(extracted, deisotoped, peaks, imzMLparse, spectra_dir, thres.int, thres.low, thres.high, files)
  ids<-cbind(deisotoped[[2]][,1],annotated,deisotoped[[2]][,3:4])
  image.ann <- cbind(ids, final.image[,2:ncol(final.image)]) # create annotated image

  # normalise image 
  image.norm<-normalise(imagedata.in=image.ann, norm.type, standards,offset=4)
  write.csv(image.norm[,], "image.norm.csv")
  write.csv(image.norm[,1:3], "image.norm_short.csv")
  image.process = image.norm
} else { image.process <- read.csv("image.norm.csv")
         image.process <- image.process[,2:ncol(image.process)]
}

    # perform PCA if requested
 if(pca==T){
    image.scale<-centreScale(imagedata.in = image.process, scale.type,transform, offset = 4)
    imagedata.in=image.scale
    imagePca(imagedata.in=image.scale, offset=4,  PCnum, scale, x.cood, y.cood, nlevels, res.spatial, summary,title,rem.outliers)
    pca <-princomp(t(imagedata.in[,(offset+1):ncol(imagedata.in)]), cor = FALSE, scores = TRUE, covmat = NULL)
    labs.all<-as.numeric(as.vector(imagedata.in[,1]))
    for (i in 1:PCnum){
      loadings <- pca$loadings[,i]
      loadings <- cbind(loadings, labs.all)
      write.csv(loadings, paste("loadings_PC",i,".csv"))
    }
  }
  
  # make ion slice if requested
  if(slice==T){
    imageSlice(row, imagedata.in=image.process, scale, x.cood, y.cood, nlevels,name=image.process[row,1],subname=image.process[row,2], offset=4, res.spatial, rem.outliers, summary,title)
  }
  
  # perform clustering if requested
  if(cluster.k==T){
    cluster(cluster.type, imagedata.in=image.process, offset, width=x.cood, res.spatial, height=y.cood, clusters)
    }

# returns normalised data file, also printed as a csv
if (process==T) {return(image.norm)} 

}