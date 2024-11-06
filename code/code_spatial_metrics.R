# In this section, we calculate two spatial metrics: the Congruence Index (COI) 
# and the Occurrence in Unsampled Cells (OUC). Both metrics rely on predictions 
# from Species Distribution Models (SDMs) for the target species, along with 
# independent data on its distribution.
 
# - Congruence Index (COI): A newly proposed index that measures the alignment, 
#   or congruence, of SDM predictions with an independent expert-defined 
#   distribution. CI provides a quantitative assessment of how closely the 
#   model’s predictions match expert insights, adding an extra layer of validation.
 
# - Occurrence in Unsampled Cells (OUC): Evaluates the presence of the species 
#   in areas predicted by the SDM that have not been sampled. This metric 
#   helps determine the model’s effectiveness in identifying suitable habitats 
#   outside of the areas with observed occurrences, providing insights into the 
#   impact of spatial sampling biases.

# Together, these metrics provide a comprehensive evaluation of SDM performance.


### 0. LOAD PACKAGES ---------------------------------------------------------------

# Define required packages
required_packages <- c(
          "tidyverse", "rgdal","sf","sp","rnaturalearth",
          "rnaturalearthdata", "raster", "maptools",
          "gdalUtilities", "tidyr", "ggspatial", "purrr",
          "ggpubr", "readxl")

# Install and load packages if not already loaded
invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# Check working directory
getwd()


### 1. CALCULATE COI METRIC ---------------------------------------------------------------

# Open the list of study species
spp_keys <- read_excel("data_inputs/spp_keys.xlsx")

# Identify incomplete species
incomplete_species <- c(
  "Sambucusnigr", "Ilexaquifoli", "Hederahelix", "Cornussangui", 
  "Juniperuscom", "Vacciniummyr", "Morusnigra", "Crataegusmon", 
  "Prunusspinos", "Rosacanina", "Solanumnigru", "Taxusbaccata", 
  "Vitisvinifer"
)

# Filter out the incomplete species from spp_keys
spp_keys <- subset(spp_keys, !plant_sp %in% incomplete_species)

# Load the proper coordinate reference system
euro_albers = CRS("+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs")

### COI metric

# Create an empty dataframe where we will upload data on the COI
# metric, the treatment combination, and the study species
combination <- NULL
COI <- NULL
COI_metric <- data.frame(combination, COI)

### Generate a for loop to repeat across species
for(j in 1:nrow(spp_keys)) { 
  
  # First, download independent data from the Atlas Florae Europaeae 
  # (AFE) regarding the distribution of study species, and save in
  # a proper folder
  # https://www.luomus.fi/en/new-grid-system-atlas-florae-europaeae
  
  # When downloaded, open shapefile with centroids from the AFE
  if(file.exists(paste("data_inputs/AFE/",
                       spp_keys[j,1],".shp",sep = ""))){
    # data from AFE could not available for all the study species
    
    # Read shapefile as a spatial feature
    points <- read_sf(paste("data_inputs/AFE/",
                            spp_keys[j,1],".shp",sep = ""))
    
    # Transform to the adequate coordinate reference system
    points <- st_transform(points, crs(euro_albers))
    
    # Set a rectangular buffer of 50x50 km around centroids
    buffer <- sf::st_buffer(points,
                            dist = 50000, # 50 km rectangular buffer
                            endCapStyle = "SQUARE")
    
    # Dissolve it and generate a spatial feature object
    buffer_dissol <- st_union(buffer)
    buffer_dissol <- st_sf(buffer_dissol)
    
    # Loop to get the value of the CI metric per each treatment 
    # (12 treatment combinations per species)
    comb <- as.character(c(1:12))
    for (i in 1:12) {
      if(file.exists(paste("data_inputs/sdm_prediction/",
                           spp_keys[j,1],"/proj_mod_",
                           comb[i],"EM_noTH.tif",sep = ""))){
        
        # Open the correspondent SDM prediction
        sdm <- raster(paste("data_inputs/sdm_prediction/",
                            spp_keys[j,1],"/proj_mod_",
                            comb[i],"EM_noTH.tif",sep = ""))
        
        # Capture the SDM prediction within the masked raster
        masked <- raster::mask(x = sdm, mask = buffer_dissol)
        
        # Calculate the mean probability of occurrence within the 
        # masked raster
        val <- getValues(masked) # get raster values
        mi <- mean(val,na.rm=T) # calculate the mean probability
        
        # Capture the SDM prediction outside the masked raster
        masked_inverse <- raster::mask(x = sdm, mask = buffer_dissol, 
                               inverse = TRUE)
        
        # Calculate the mean probability of occurrence outside 
        # the masked raster (= outside the real - AFE - 
        # distribution)
        val_inv <- getValues(masked_inverse) # get raster values
        mo <- mean(val_inv,na.rm=T) # calculate the mean probability
        
        # Calculate the COI metric
        # COI is defined within the open (0, 1) interval since mi>mo, 
        # where mi is the mean probability of occurrence within the 
        # real distribution of the species and mo is the mean 
        # probability of occurrence outside the real distribution of 
        # the species
        COI <- 1 - (mo/mi)
        
        # Larger COI values indicate better congruence of the models 
        # (spatial coincidence between AFE range and model 
        # prediction)
        
        # Add information on treatment combination and species
        combination <- i
        species <- paste(spp_keys[j,1])
        add <- data.frame(combination, species, COI)
        
        # Add to main dataframe
        COI_metric <- rbind(COI_metric, add)
      }
    }
  }
}



### 2. CALCULATE THE OUC METRIC ---------------------------------------------------------------

# Open the list of study species
spp_keys <- read_excel("data_inputs/spp_keys.xlsx")

# Open sampling effort raster
samplingeffort <- raster("data_inputs/samplingeffort_10km.tif")

# Reclassify to assign a value to NA
samplingeffort <- reclassify(samplingeffort, cbind(NA, 0))
plot(samplingeffort)

# Generate base dataframe
combination <- NULL
species <- NULL
mean <- NULL
OUC_metric <- data.frame(combination, species, mean)

# Loop to get the value of the OUC metric per each treatment 
# (12 treatment combinations per species)
comb <- as.character(c(1:12))
for(j in 1:nrow(spp_keys)){
  for(i in 1:length(comb)){
    if(file.exists(paste("data_inputs/sdm_prediction/",
                         spp_keys[j,1],"/proj_mod_",
                         comb[i],"EM_noTH.tif",sep = ""))){
      
      # Open the correspondent SDM prediction
      dm <- raster(paste("data_inputs/sdm_prediction/",
                         spp_keys[j,1],"/proj_mod_",
                         comb[i],"EM_noTH.tif",sep = ""))
    
      # Generate a raster stack with the sampling effort in the 
      # study area and the SDM prediction
      s <- stack(dm, samplingeffort)
    
      # Transform object to dataframe
      v <- data.frame(na.omit(values(s)))
    
      # Exclude sampling efforts over 0
      v <- v[v$samplingeffort_10km <= 0, ]
    
      # Rename columns accordingly
      colnames(v) <- c("value", "samplingeffort")
    
      # Obtain the mean probability of occurrence in unsampled cells
      mean <- mean(v$value) # Reminder: we are interested on this 
                          # mean probability in order to check the 
                          # effect of treatments in areas without
                          # sampling effort
    
      # Get 0 if NAs appear
      mean <- ifelse(is.na(mean),0,mean)
    
      # Add information on treatment combination and species
      combination <- i
      species <- paste(spp_keys[j,1])
      add <- data.frame(combination, species, mean)
    
      # Add to main dataframe
      OUC_metric <- rbind(OUC_metric, add)
    }
  }
}

