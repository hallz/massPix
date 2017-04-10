massPix
==========
Processes high resolution mass spectrometry imaging data, performs multivariate statistics (PCA, clustering) and lipid identification.

This is the initial README file for the massPix repository. 

Currently this processes one file at a time.

Supports .imzML files - if your files are currently in this format proceed to "processing". Otherwise see "raw_to_imzML_quickstart.pptx" for instructions.

See also powerpoint presentation with step by step instructions to get started - "massPix_quickstart.pptx"

massPix is a package which can be called by a single function (massPix) using processing_script.R (see "processing" section for details).

An R script has also been provided which shows the composite functions of massPix being used one by one (massPix_step_by_step.R). This may be useful for troubleshooting and/or adapting the source code. See "massPix Step-byStep" section below)


Notes
-------------------
Download files from Github. 

The latest version of imzMLConverter can also be downloaded from http://www.cs.bham.ac.uk/~ibs/imzMLConverter/# (Race et al, J Proteomics, 2012)  - ensure this is saved in massPix project folder, where the "data/" and "libraries/" folders are located. Ensure folder name is "imzMLConverter/".

In addition, the full test_POS.ibd file (~160 MB) must be downloaded from the dropbox link https://www.dropbox.com/s/5ff7yi3z1irqda0/test_POS.ibd?dl=0 - ensure download completes. 

Keep both .imzML and .ibd data files together in the "data/" folder. 



processing
----------

Use the processing_script.R to set your working directory and set processing parameters.

See manual files or scroll down for meaning of each parameter.

If you want to process from scratch chose process=T

Output is image.norm.csv (full normalised image data) and image.norm_short.csv (containing lipid annotations and m/z only).

If you want to perform PCA analysis choose pca=T 

Output is R plots to view/export and loadingsPC.csv files.

If you want to perform clustering choose cluster.k=T 

Output is R plots to view/export and Cluster.csv files.

If you want to view the distribution of a single ion choose slice=T 

The m/z row to be sliced is set in "row". Find this by looking at output file image.norm_short.csv. Output is R plot to view/export.



 ``` 
 ######## setting up directories #####

## Set your home directory and ensure this has three folders:
## 1. "libraries" with the library files inside
## 2. "data" with your data file to be processed inside: this should be two files (.imzML and .ibd)
## 3. "imzMLConverter" with the imzML converter subfolders and files inside
## These can all be downloaded from github

library(massPix)

home_dir <- paste(getwd(), "/", sep="") 
lib_dir <- paste(home_dir,"libraries",sep="") # library files
spectra_dir <-  paste(home_dir,"data",sep="") # spectra files
imzMLparse <- paste(home_dir,"imzMLConverter/imzMLConverter.jar",sep="") # imzML converter


# running the main function and setting options

results <- massPix(# what type of data analysis do you want?
                      process=T, 
                      pca=T,
                      slice=T,
                      cluster.k=T,
                      
                      # some important general settings
                      ionisation_mode = "positive",
                      thres.int = 100000,
                      thres.low = 200,
                      thres.high = 1000,
                      bin.ppm = 10,
                      thres.filter = 11,
                      ppm.annotate = 10,
                      res.spatial = 50,
                      
                      # settings for PCA, clustering and  ion slicing on
                      PCnum = 5,
                      row = 20,
                      cluster.type = "kmeans",
                      clusters = 5,
                      
                      # settings for deisotoping
                      ppm = 3,
                      no_isotopes = 2,
                      prop.1 = 0.9,
                      prop.2 = 0.5,
                                                                  
                      # some general settings which can probably be left as they are
                      search.mod = T,
                      mod = c(NL = T, label = F, oxidised = F, desat = F),
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
                      imagedata.in,
                      norm.type = "TIC",
                      standards = NULL,
                      scale.type = "cs",
                      transform = F,
                      scale = 100,
                      x.cood,
                      y.cood,
                      nlevels = 50,
                      name = "",
                      subname = "",
                      rem.outliers = "only",
                      summary = T,
                      title = T,
                      offset=4)

 ```
 
 
{process}{T to process data from imzML, F to skip processing}

{pca}{T to conduct PCA image analysis}

{slice}{T to generate two dimensional image of one ion defined in row}

{cluster.k}{T to perform clustering}

{ionisation_mode}{Choose "positive" or "negative"; determines which lipid classes to use when building database}

{thres.int}{Defines if intensity threshold, above which ions are retained}

{thres.low}{Defines the minumum m/z threshold, above which ions will be retained}

{thres.high}{Defines the minumum m/z threshold, below which ions will be retained}

{bin.ppm}{Mass accuracy for binning between spectra}

{thres.filter}{Defines threshold for proportion of missing values (this is the step number, not the actual proportion)}

{ppm.annotate}{Defines if ppm threshold for which |observed m/z - theoretical m/z| must be less than for annotation to be retained}

{res.spatial}{spatial resolution of image}

{PCnum}{number of PCs to calculate}

{row}{indices of row which corresponds to spectral data to plot in image slice - use image.norm_short.csv to find row number}

{cluster.type}{At the moment only supports "kmeans"}

{clusters}{Number of clusters to use for clustering}

{ppm}{Tolerance (ppm) within which mass of isotope must be within .}

{no_isotopes}{Number of isotopes to consider (1 or 2)}

{prop.1}{Proportion of monoisotope intensity the 1st isotope intensity must not exceed}

{prop.2}{Proportion of monoisotope intensity the 2nd isotope intensity must not exceed}

{search.mod}{Search modifications T/F.}

{mod}{modifications to search eg. c(NL=T, label=F, oxidised=T,desat=T)}

{lookup_mod}{A dataframe defining modifications}

{adducts}{vector of adducts to be searched in the library of lipid masses}

{sel.class}{A vector defining classes of lipids to be included in the library}

{fixed}{Defines if one of the SN positions is fixed, default is F.}

{fixed_FA}{Defines the name of the fixed FA is fixed is T, e.g. 16, 16.1, 18.2.}

{lookup_lipid_class}{A data.frame defining lipid classes}

{lookup_FA}{A dataframe defining FAs}

{lookup_element}{A dataframe defining elements}

{files}{a vector of file names; currently only supports processing one file at a time}

{spectra_dir}{Defines the path to the spectral files,}

{imzMLparse}{path to imzMLConverter}

{percentage.deiso}{Defines the proportion of total pixels to select, at random from the first file to produce a subset of the image}

{steps}{Sequence of values between 0 and 1 that define the thresholds of missing values to test}

{imagedata.in}{Dataframe containing image}

{norm.type}{"TIC" or "median" or "standards" or "none", recommend "TIC"}

{standards}{default is NULL, vector of row indices corresponding to variables that are standards}

{scale.type}{options are "cs" center and pareto scale, "c" center only, "pareto" pareto only and "none"}

{transform}{T/F to perform log transform}

{scale}{range of scale that intensity values will be scaled to}

{x.cood}{width of image}

{y.cood}{height of image}

{nlevels}{Graduations of colour scale.}

{name}{main name of image.}

{subname}{sub name of image}

{rem.outliers}{Choose "only" to show only images after removing outliers, T for both with and without outliers}

{summary}{T/F to show extra information}

{title}{show titles, T/F"}

{offset}{number of columns of text preceding data}


#############

massPix Step-by-Step
------------------------

For troubleshooting and/or adapting the source code.

```

############################  massPix step by step ##################################################


######## setting up directories #####

## Set your working directory and ensure this has three folders:
## 1. "libraries" with the library files inside
## 2. "data" with your data file to be processed inside: this should be two files (.imzML and .ibd)
## 3. "imzMLConverter" with the imzML converter subfolders and files inside
## These can all be downloaded from github

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
    

 
```

