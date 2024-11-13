### Load libraries
################################################################
library(raster)
library(dismo)
library(sf)
library(sp)
library(rgeos) 
library(biomod2) # For this study we used version 3.5.1
library(ggplot2)
library(RColorBrewer) 
library(dplyr); library(tidyr)
library(sp)


### Loop por especie
setwd(sdmoutput)


# Load future data

time <- c("2021-2040","2081-2100")
time_name <- c("t21","t81")
scenario <- c("ssp126","ssp370","ssp585")
bio <- c("bio2","bio4","bio8","bio9","bio10","bio15","bio18","bio19")
namesbio <- c("bio_2","bio_4","bio_8","bio_9","bio_10","bio_15","bio_18","bio_19")
mod <- "mod_11"

### Variables ambientales (predictores)

# soil <-stack(list.files(path="../data_sources/environmental_variables/selection_10km/",pattern='tif',full.names=TRUE))
# soil <- soil[[9:14]]
# 
# files <- NULL
# for(t in 1:length(time)){
#   for(s in 1:length(scenario)){
#     for(b in 1:length(bio)){
#       file_tif <- paste("../data_sources/environmental_variables/future_climate/WORLDCLIM/",time[t],"/STACKS/",scenario[s],"/",bio[b],"/mean.tif", sep="")
#       files <- c(files,file_tif)
#     }
#     bioclim<-stack(files)
#     names(bioclim) <- namesbio
#     myExpl <- stack(bioclim,soil)
#     writeRaster(myExpl,filename = paste("../data_sources/environmental_variables/future_climate/stack/",time[t],scenario[s],".grd",sep=""))
#     files <- NULL
#     }
# }



for(i in ini:fin) {
  for(t in 1:length(time)){
      for(s in 1:length(scenario)){
        
  myRespName <- species_all[i]
  cat("\n######## \n Ensemble modelling for", myRespName, "\n Time:", date())  
  DataSpecies <- read.csv(paste("../data_outputs/GBIF_coords_corrected_spatial_bias/sp_correct_bias",myRespName,".csv",sep = ""))
  
  myRespName <- stringr::str_replace_all(myRespName, " ", "") #eliminamos espacio para evitar problemas al guardar archvos
  myRespName <- stringr::str_sub(myRespName, start=1,end=12) 


  # Environmental data
  
  exp_fut <- stack(paste("../data_sources/environmental_variables/future_climate/stack/",time[t],scenario[s],".grd",sep=""))
  crs(exp_fut) <- CRS("+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs")
  exp_fut <- subset(exp_fut,c("bio_10","bio_15","bio_18","bio_19","bio_2","bio_4","bio_8","bio_9","slope_def","soilcfvo_def","soilclay_def","soilpH_def","soilsilt_def"))

  
  ### 4.- Calibracion de los modelos ###
  ######################################
  
  ### Extraemos objeto de la carpeta
  # comprobamos si el modelo se ejecutó correctamente
  if(file.exists(paste(sdmoutput,"/" ,myRespName,"/",myRespName,".",mod,".models.out",sep=""))==F){next}
  load(paste(sdmoutput,"/" ,myRespName,"/",myRespName,".",mod,".models.out",sep=""))
  myBiomodModelOut <- get(paste(myRespName,".",mod,".models.out",sep=""))
  rm(list=ls(pattern=paste(myRespName,".",mod,".models.out",sep="")))

 
  ### 6.- Ensemble modeling ###
  #############################
  
  ### Extraemos objeto de la carpeta
  if(file.exists(paste(sdmoutput,"/" ,myRespName,"/",myRespName,".",mod,".ensemble.models.out",sep=""))==F){next}
  load(paste(sdmoutput ,myRespName,"/",myRespName,".",mod,".ensemble.models.out",sep=""))
  myBiomodEM <- get(paste(myRespName,".",mod,".ensemble.models.out",sep=""))
  rm(list=ls(pattern=paste(myRespName,".",mod,".ensemble.models.out",sep="")))
  

  ### 7.- Proyeccion de modelos "actuales" (al espacio geografico)###
  ##################################################################
  
  ### Proyeccion sobre el area de estudio bajo condiciones futura de cada modelo individual
  myBiomodProj <- BIOMOD_Projection(
    bm.mod = myBiomodModelOut, # Resultados de los modelos
    new.env = exp_fut, # Variables ambientales para toda Europa
    proj.name = paste(mod,time_name[t],scenario[s],sep = ""), # Nombre de las proyecciones
    models.chosen = 'all', # models to be projected
    metric.binary = 'all',
    metric.filter= 'all',
    compress = TRUE,
    output.format = '.img', # Formato de archivos GIS (tambien *.img)
  )
  


  ### 8.- Proyeccion de modelos de conjunto ###
  #############################################
  
  myBiomodEF <- BIOMOD_EnsembleForecasting( 
    bm.em = myBiomodEM,
    bm.proj = myBiomodProj,
    compress = 'gzip',
    output.format = '.img', # Formato de archivos GIS (tambien *.img)
    do.stack=FALSE,
    build.clamping.mask = TRUE,
    metric.binary = 'all',
    metric.filter= 'all',
    seed.val = seed,
    nb.cpu=cpu,
    do.stack=FALSE)
  
 
  
      }  
  }  
}   
    
