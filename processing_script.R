setwd("D:/MALDI_data/")

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


# running the main function and setting the options

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
