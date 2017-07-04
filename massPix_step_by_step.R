############################  massPix step by step ##################################################


######## setting up directories #####

## Set your working directory and ensure this has three folders:
## 1. "libraries" with the library files inside
## 2. "data" with your data file to be processed inside: this should be two files (.imzML and .ibd)
## 3. "imzMLConverter" with the imzML converter subfolders and files inside


library(massPix)

home_dir <- paste(getwd(), "/", sep="") 
lib_dir <- paste(home_dir,"libraries",sep="") # library files
spectra_dir <-  paste(home_dir,"data",sep="") # spectra files
imzMLparse <- paste(home_dir,"imzMLConverter/imzMLConverter.jar",sep="") # imzML converter


##### setting parameters #####

  # some important general settings
  ionisation_mode = "positive"
  thres.int = 100000
  thres.low = 200
  thres.high = 1000
  bin.ppm = 10
  thres.filter = 11
  ppm.annotate = 10
  res.spatial = 50
  
  # settings for PCA, clustering and ion slicing 
  PCnum = 5
  row = 12
  cluster.type = "kmeans"
  clusters = 5
  
  # settings for deisotoping
  ppm = 3
  no_isotopes = 2
  prop.1 = 0.9
  prop.2 = 0.5
  
  # some general settings which can probably be left as they are
  search.mod = T
  mod = c(NL = T, label = F, oxidised = F, desat = F)
  fixed = F
  percentage.deiso = 3
  steps = seq(0, 1, 0.05)
  norm.type = "TIC"
  standards = NULL
  scale.type = "cs"
  transform = F
  scale = 100
  nlevels = 50
  name = ""
  subname = ""
  rem.outliers = "only"
  summary = T
  title = T
  offset=4

####################################################################################################################


########## massPix function broken down into composite functions ################


## read and process very large files

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

###  main script ###

  # make library

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


  # perform PCA if requested

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
  
  
  # make ion slice if requested

  imageSlice(row, imagedata.in=image.process, scale, x.cood, y.cood, nlevels,name=image.process[row,1],subname=image.process[row,2], offset=4, res.spatial, rem.outliers, summary,title)
  
  
  # perform clustering if requested

  cluster(cluster.type, imagedata.in=image.process, offset, width=x.cood, res.spatial, height=y.cood, clusters)
    

 

