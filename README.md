massPix
==========
This is the initial README file for the massPix repository. 

massPix processes high resolution mass spectrometry imaging data, performs multivariate statistics (PCA, clustering) and lipid identification.

massPix is an R package which can be called by a single function (massPix) using "processing_script.R". Currently massPix processes one file at a time. massPix supports .imzML files - if your files are currently in this format proceed to read the "processing" section. For instructions on on converting raw data files to imzML format, the "raw_to_imzML_quickstart.pptx" file has further instructions.

Detailed step-by-step instructions on get started with using massPix are available in the "massPix_quickstart.pptx" file.

An R script has also been provided which shows the composite functions of massPix being used one by one ("massPix_step_by_step.R"). This may be useful for troubleshooting and/or adapting the source code. See "massPix Step-byStep" section below).


Brief introduction
-------------------
massPix is run from the R scripting interface (function: massPix), however a detailed knowledge of R is not required to install and use the software. Those with advanced knowledge of R programming can adapt the source code for their own needs. massPix outputs high quality images, a dataframe of the final normalised and annotated image which can be further manipulated in R, and csv files for spectra corresponding to cluster centers, PCA loadings, and lipid annotations. The massPix R package ("massPix_1.2.tar.gz"), all R scripts (wrapper, and individual functions), library files (for creating the lipid database) and the imzML Converter are available here. A more detailed step-by-step presentation on software use ("massPix_quickstart.ppt") and instructions on file conversion ("raw_to_imzML_quickstart.ppt") are also available here. Test data will shortly be available on the MetaboLights repository (accession number TBC). 

The overall data processing workflow consists of initial data pre-processing, filtering, image subsetting, deisotoping, annotation, normalisation, scaling, image “slicing” and multivariate statistics. Data in imzML format is parsed to R. Ions with intensities greater than a threshold, from each spectra, are extracted and grouped to user-adjustable mass bins (function: mzextractor). Spectral features are defined by the median m/z value in each bin, and only features detected above a threshold proportion of spectra are retained (function: peakpicker.bin). Average intensities for all features from a random subset of pixels are computed and used to perform deisotoping (function: subsetImage). The deisotoping algorithm identifies the molecular ion (M) and removes isotopes at m/z (M+1) and (M+2) which are within a calculated proportion of the intensity of M (function: filtered and deisotope). 

Putative lipid annotation by accurate mass is achieved by searching deisotoped ions against a generated library of lipid m/z ratios computed for all combinations of common fatty acids, lipid head-groups and anticipated adducts in each ionisation mode (functions: makelibrary and annotate). The criteria for a match can be adjusted according to different MS performance capabilities. Lipid classes searched in positive ion mode are diacylglycerides (DAG), triacylglycerides (TAG), phosphatidylcholines (PC), phosphatidylethanolamines (PE), phosphatidylserines (PS), LysoPC, cholesteryl esters (CE), sphingomyelins (SM) and ceramides (Cer). In negative ion mode, lipid classes searched are PC, phosphatidic acid (PA), PE, PS, phosphatidylglycerols (PG), phosphatidylinositols (PI), and free fatty acids (FFA). Whilst this list is not exhaustive, it does cover the most common lipid classes. Possible adducts considered are [M+K]+, [M+H]+, [M+Na]+, [M+NH4]+ in positive ion mode and [M–H]-, [M+Cl]-, [M+OAc]- in negative ion mode. It is important to point out that a database hit based on accurate mass should only be considered the first step in metabolite identification, and confirmation carried out using MS/MS is required, where this appropriate. 

massPix has the further capability to perform difference matching on deisotoped features to search for mass differences associated with measurement-introduced alternation (e.g. loss of headgroup, loss of water) or biological modifications (e.g. oxidation). The final image is constructed, based on all pixels (function: constructImage). Ion intensities are then normalised either to the median or total ion count, or to the average intensity of a set of standard ions (function: normalise). Single ion images can be produced (function: imageSlice) , or normalised intensities used to create multivariate statistical images based on k-means clustering (function: cluster) or PCA following centering and Pareto scaling (function: centreScale and imagePca). The analysis can be readily customised by replacing default parameters for filtering, normalisation and scaling, library composition, lipid assignment and image reporting.



Notes
-------------------
Download files from Github. 

The latest version of imzMLConverter can also be downloaded from http://www.cs.bham.ac.uk/~ibs/imzMLConverter/# (Race et al, J Proteomics, 2012)  - ensure this is saved in massPix master folder, where the "data/" and "libraries/" folders are located. Ensure folder name is "imzMLConverter/".

Test data is comprised of two files: test_POS.idb and test_POS.imzML file. The imzML fie is available to download on github, however test_POS.ibd file (~160 MB) must be downloaded from the dropbox link https://www.dropbox.com/s/5ff7yi3z1irqda0/test_POS.ibd?dl=0 - ensure download completes. These files will both be shortly available on MetaboLights repository.

Keep both .imzML and .ibd data files together in the "data/" folder. 



processing
----------

Use the "processing_script.R" to set your working directory and set processing parameters.

See below for parameter descriptions.


 ``` 
 ######## setting up directories #####

## Set your home directory and ensure this has three folders:
## 1. "libraries" with the library files inside
## 2. "data" with your data file to be processed inside: this should be two files (.imzML and .ibd)
## 3. "imzMLConverter" with the imzML converter subfolders and files inside


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
 
 
Process - T to process data from imzML file. F to skip processing. To process data for the first time choose process=T. You can use process=F, when wishing to re-process data with different settings for PCA, clustering and/or ion slicing. The final output is image.norm.csv (full normalised image data matrix - lipid annotations, m/z and intensities) and image.norm_short.csv (containing lipid annotations and m/z only).

pca - T to conduct PCA image analysis. If you want to perform PCA analysis choose pca=T; output are R plots of PC scores and loadings which can be view or exported as pdf. In addition the PC loadings are stored as loadingsPC.csv files.

slice - T to generate two dimensional image of one ion. The ion to be "sliced" is defined in parameter "row" (see below). Output is R plot which may be viewed/exported as pdf.

cluster.k - T to perform clustering. Output are R plots of clustered pixels and spectra corresponding to the average of each cluster. In addition spectra are written to Cluster.csv files.

ionisation_mode - Choose "positive" or "negative"; this will determine which lipid classes are used when building database

thres.int - defines intensity threshold, above which ions are retained

thres.low - defines the minimum m/z threshold, above which ions will be retained

thres.high - Defines the maximum m/z threshold, below which ions will be retained

bin.ppm - Mass accuracy for binning between spectra 

thres.filter - Defines threshold for proportion of missing values - in steps ranging from 1 to 21 where 0 is no ions retained (ion must be present in every pixel)and 21 is all ions retained (ion must be present in at least one pixel); suggest using values between 11 and 15.

ppm.annotate - Defines ppm threshold for which |observed m/z - theoretical m/z| must be less than for annotation to be retained. 

res.spatial - spatial resolution of image; units are micrometers

PCnum - number of principal components to calculate during PCA

row - indices of row which corresponds to spectral data to plot in image slice. The row number of the m/z to slice can be found by viewing the output file "image.norm_short.csv".

cluster.type - Currently only supports "kmeans" clustering

clusters - defines the number of clusters to use for clustering

ppm - Tolerance (ppm) within which mass of isotope must be within. 

no_isotopes - Number of isotopes to consider (1 or 2)

prop.1 - Proportion of monoisotope intensity the 1st isotope intensity must not exceed

prop.2 - Proportion of monoisotope intensity the 2nd isotope intensity must not exceed

search.mod - Search modifications T/F.

mod - modifications to search eg. c(NL=T, label=F, oxidised=T,desat=T). This step searches for modifications (defined in library file "lib_modifiction.csv") including neutral loss (loss of water, choline, phosphocholine, etc), label (13C palmitate) (oxidation and desaturations.

lookup_mod - A dataframe defining modifications - this is read from the library file "lib_modification.csv" and need not be defined by user

adducts - vector of adducts to be searched in the library of lipid masses. Pre-defined in annotate function. If ionisation mode is positive, adducts <- c(H = T, NH4 = F, Na = T, K = T, dH = F, Cl = F, OAc = F). If ionisation mode is negative, adducts<- c(H = F, NH4 = F, Na = F, K = F, dH = T, Cl = T, OAc = F)

sel.class -  A vector defining classes of lipids to be included in the library. Pre-defined in makelibrary function. If ionisation mode is positive, lipid classes searched are: TG, DG, PC, PE, PS, LysoPC, DG-H20, CE, SM, Cer. If ionisation mode is negative, lipid classes searched are:  PC, PA, PS, PE, PG, PI, PIP, PIP2, PIP3, FFA 

fixed - Defines if one of the SN positions is fixed, default is F.

fixed_FA - Defines the name of the fixed FA if fixed is T, e.g. 16, 16.1, 18.2.

lookup_lipid_class - A dataframe defining lipid classes - this is read from the library file "lib_class.csv" and need not be defined by user

lookup_FA - A dataframe defining fatty acids - this is read from the library file "lib_FA.csv" and need not be defined by user

lookup_element - A dataframe defining elements - this is read from the library file "lib_element.csv" and need not be defined by the user

files - a vector of file names; currently only supports processing one file at a time. Does not need to be defined by user.

spectra_dir - file path to the spectral files

imzMLparse - file path to imzMLConverter

percentage.deiso - Defines the proportion of total pixels to select, at random from the first file to produce a subset of the image. This step speeds up the image processing, by performing deisotoping on a subset of the image and then applying it to the whole image.   

steps - Sequence of values between 0 and 1 that define the thresholds of missing values to test

imagedata.in - Dataframe containing image - does not need to be defined by user.

norm.type - "TIC" or "median" or "standards" or "none" are the options. These are different ways to normalise the data, suggest to use "TIC". "standards" may be used when using an internal standard(s)

standards - default is NULL, vector of row indices corresponding to variables that are standards. When wishing to normalise to an internal standard(s), choose norm.type = "standards". The row number for the m/z of the standard(s) must be provided, e.g. standards = c(45, 60). Row indices can be found by looking at the image_norm_short.csv output file. 

scale.type - options are "cs" center and pareto scale, "c" center only, "pareto" pareto only and "none"}

transform - T/F to perform log transform

scale - range of scale that intensity values will be scaled to, typically 100.

x.cood - width of image in pixels - read from data file and no need to be defined by user

y.cood - height of image in pixels - read from data file and no need to be defined by user

nlevels - Graduations of colour scale.

name - main name of image.

subname - sub name of image

rem.outliers - Choose "only" to show only images after removing outliers, T for both with and without outliers

summary - T/F to show extra information

title - show titles, T/F

offset - number of columns of text preceding data


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
## 

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

Scripts for all the composite functions above are found in the R directory here.
