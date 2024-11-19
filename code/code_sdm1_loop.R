### Load libraries
################################################################
library(raster)
library(dismo)
library(sf)
library(sp)
library(rgeos) 
library(biomod2) # For this study we used version 4.1-2
library(ggplot2)
library(RColorBrewer) 
library(dplyr); library(tidyr)
library(sp)


sdmoutput <-  ""  # Work directory
setwd(sdmoutput)

### Environemntal variables. Different set depending on the bias correction
Exp <-stack(list.files(path="../data_sources/environmental_variables/selection_10km/",pattern='tif',full.names=TRUE))
crs(Exp) <- CRS("+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs")
# plot(Exp)

# ExpCou <-stack(list.files(path="../data_sources/environmental_variables/variables_biascorrect_country/",pattern='tif',full.names=TRUE))
# crs(ExpCou) <- CRS("+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs")
# # plot(ExpCou)
# 
# ExpPix <-stack(list.files(path="../data_sources/environmental_variables/variables_biascorrect_pixel/",pattern='tif',full.names=TRUE))
# crs(ExpPix) <- CRS("+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs")
# # plot(ExpPix)

# Species
plant_keys <- read.csv("../data_outputs/plant_keys.csv")
# plant_keys <- read.csv("C:/Users/claud/Desktop/UCADOCTORADO/REPOSITORIOS/chapter_3/data_outputs/plant_keys.csv")
species_all <- plant_keys[,"plant_sp"]

### Opciones de calibracion ###
###################################
seed <- 8739465
set.seed(seed)   
cpu=2
maxent_file <- "C:/Users/USUARIO/Documents/R/maxent/maxent.jar"
pa.nb.rep = 5
NbRunEval = 10

### GAM: 
myBiomodOption <- BIOMOD_ModelingOptions(GAM = list(k = 4),
                                         MAXENT.Phillips = list(path_to_maxent.jar=maxent_file))


print(myBiomodOption)

### Loop per species
setwd(sdmoutput)

for(i in ini:fin) {
 
  # Get species occurrence data
   myRespName <- species_all[i]
  cat("\n######## \n Ensemble modelling for", myRespName, "\n Time:", date())  
  DataSpecies <- read.csv(paste("../data_outputs/GBIF_coords_corrected_spatial_bias/sp_correct_bias",myRespName,".csv",sep = ""))
  myRespName <- stringr::str_replace_all(myRespName, " ", "") 
  myRespName <- stringr::str_sub(myRespName, start=1,end=12) 


  spPointData <- SpatialPoints(DataSpecies[,c("coords.x1","coords.x2")], 
                               proj4string = CRS("+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"))
  
  myResp <- spPointData
  Npresences <- nrow(DataSpecies)
  

  # Combinations of parameters to be tested
  myGrid <- data.frame(rbind(expand_grid(strategy = c('disk','random'),
                                         pa.dist.min = c(200000,""),
                                         pa.dist.max = "",
                                         PA.nb.absences = c(1000,Npresences),
                                         expl.var = c("Exp","ExpCou","ExpPix")))) %>%
    filter(!(strategy=="random" & (pa.dist.min!=""|pa.dist.max!=""))) %>%
    filter(!(strategy=="disk" & (pa.dist.min=="")))
  
  myGrid[is.na(myGrid$pa.dist.max)&!is.na(myGrid$pa.dist.min),"pa.dist.max"] <- ""
  myGrid$pa.dist.min <- as.numeric(myGrid$pa.dist.min)
  myGrid$pa.dist.max <- as.numeric(myGrid$pa.dist.max)
  
  write.csv(myGrid,paste0(sdmoutput,"/" ,myRespName,"/list_models_",myRespName,".csv",sep=""))                  

  for(j in 1:nrow(myGrid)){
 
     ### 2.- Data formatting ###
  ################################


    IDmodel <- paste0("mod_",j)
  
  
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                       expl.var = get(myGrid$expl.var[j]),
                       resp.name = myRespName,
                       eval.resp.var = NULL, 
                       eval.expl.var = NULL,
                       PA.nb.rep = pa.nb.rep, 
                       PA.nb.absences = myGrid$PA.nb.absences[j],
                       PA.strategy = myGrid$strategy[j],
                       PA.dist.min = ifelse(is.na(myGrid$pa.dist.min[j]),"",myGrid$pa.dist.min[j]),
                       PA.dist.max = ifelse(is.na(myGrid$pa.dist.max[j]),"",myGrid$pa.dist.max[j]),
                       PA.sre.quant = 0.025,
                       na.rm = TRUE)


  if(ncol(myBiomodData@PA.table)==1) {
    print(paste("There is no enough geographical space to draw all PA in: ",myRespName,IDmodel))
    
    #Save data format
    if(file.exists(myRespName)==F){dir.create(myRespName)}
    pdf(file = paste(myRespName,"/",myRespName,"_",IDmodel,"_Data.pdf",sep="") )  
        plot(myBiomodData)
    dev.off()
    next
    }
  
  
  ### Extraemos objeto de la carpeta
  # load(paste(sdmoutput,myRespName,"/",myRespName,".",IDmodel,".models.out",sep=""))
  # myBiomodModelOut <- get(paste(myRespName,".",IDmodel,".models.out",sep=""))
  # rm(list=ls(pattern=paste(myRespName,".",IDmodel,".models.out",sep="")))
  

  
  ### 4.- Model calibration ###
  ######################################
  

  myBiomodModelOut <- BIOMOD_Modeling(
    bm.format=myBiomodData, 
    models = c('GLM', 'GBM','GAM','ANN',
               'RF','MAXENT.Phillips.2'), 
    bm = myBiomodOption,
    nb.rep=NbRunEval, 
    data.split.perc=70,
    prevalence=0.5,
    var.import=3, 
    metric.eval = c('TSS','ROC','KAPPA'), 
    save.output = TRUE, 
    do.full.models = FALSE, 
    modeling.id = IDmodel,
    seed.val = seed,
    nb.cpu=cpu) 
 
  
  pdf(file = paste(myRespName,"/",myRespName,"_",IDmodel,"_Data.pdf",sep="") ) 
    plot(myBiomodData)
  dev.off()
  
  

  
   ### Extraemos objeto de la carpeta
  # load(paste(sdmoutput,myRespName,"/",myRespName,".",IDmodel,".models.out",sep=""))
  # myBiomodModelOut <- get(paste(myRespName,".",IDmodel,".models.out",sep=""))
  # rm(list=ls(pattern=paste(myRespName,".",IDmodel,".models.out",sep="")))



  ### 5.- Model asssessment ###
  #####################################
  
  options(max.print=999999) # Dar mas memoria RAM
  
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)

  print(fivenum(as.numeric(myBiomodModelEval["ROC","Testing.data",,,])))
  
  print(fivenum(as.numeric(myBiomodModelEval["TSS","Testing.data",,,])))
  
  get_variables_importance(myBiomodModelOut)
  
  var.imp.mod.out <- get_variables_importance(myBiomodModelOut, as.data.frame=T)
  head(var.imp.mod.out, n = 5)
  
  ### Save importance of variables
  capture.output(get_variables_importance(myBiomodModelOut),
                 file=file.path(myRespName, 
                                paste(myRespName,"_",IDmodel,"_variables_importance.csv", sep="")))
  
  ### Save ROC and TSS metrics
  evalDF.ROC <- as.data.frame(myBiomodModelEval["ROC","Testing.data",,,])
  evalDF.TSS <- as.data.frame(myBiomodModelEval["TSS","Testing.data",,,])
  
  write.csv(evalDF.ROC, file = paste(myRespName,"/",myRespName,"_",IDmodel,"_evalDF_ROC.csv",sep=""))
  write.csv(evalDF.TSS, file = paste(myRespName,"/",myRespName,"_",IDmodel,"_evalDF_TSS.csv",sep=""))


  ### 6.- Ensemble modeling ###
  #############################
  

  myBiomodEM <- BIOMOD_EnsembleModeling( 
    bm.mod = myBiomodModelOut, 
    models.chosen = 'all', 
    em.by='all', 
    metric.select =  c('TSS','ROC'), 
    metric.select.thresh  = c(0.7,0.7),
    metric.eval = c('TSS'), 
    prob.mean = TRUE, 
    prob.cv = TRUE, 
    prob.ci = FALSE, 
    prob.ci.alpha = 0.05, 
    prob.median = FALSE, 
    prob.mean.weight = TRUE, 
    prob.mean.weight.decay = 'proportional',  
    seed.val = seed,
    nb.cpu=cpu)

  
  ### Extraemos objeto de la carpeta
  # load(paste(sdmoutput ,myRespName,"/",myRespName,".",IDmodel,".ensemble.models.out",sep=""))
  # myBiomodEM <- get(paste(myRespName,".",IDmodel,".ensemble.models.out",sep=""))
  # rm(list=ls(pattern=paste(myRespName,".",IDmodel,".ensemble.models.out",sep="")))


  ### 7.- Current model projection ###
  ##################################################################
  
myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut, 
  new.env = Exp, 
  proj.name = paste0(IDmodel,'current'),
  models.chosen = 'all', 
  metric.binary = 'all',
  metric.filter= 'all',
  compress = TRUE,
  output.format = '.img', 
  )
  
  
  # load(paste(sdmoutput ,myRespName,"/proj_", IDmodel, "current/",myRespName,".",IDmodel,"current.projection.out",sep=""))
  # myBiomodProj <- get(paste(myRespName,".",IDmodel,"current.projection.out",sep=""))
  # rm(list=ls(pattern=paste(myRespName,".",IDmodel,"current.projection.out",sep="")))#
  
  ### 8.- Project ensemble model ###
  #############################################
  
  myBiomodEF <- BIOMOD_EnsembleForecasting( 
    bm.em = myBiomodEM,
    bm.proj = myBiomodProj,
    compress = 'gzip',
    output.format = '.img', 
    do.stack=FALSE,
    build.clamping.mask = TRUE,
    metric.binary = 'all',
    metric.filter= 'all',
    seed.val = seed,
    nb.cpu=cpu,
    do.stack=FALSE)
  
  ### Print summary
  pdf(file = paste(myRespName,"/",myRespName,"_",IDmodel,"_ensemble.pdf",sep="") ) 
    print(myBiomodEF)
    plot(myBiomodEF)
  dev.off()


 
  }
}  
