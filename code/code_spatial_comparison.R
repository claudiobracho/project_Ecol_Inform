# In this section, we compute spatial comparisons between treatments
# implemented in our Species Distribution Models (SDMs). These comparisons
# assess similarity and correlation across treatment scenarios using multiple
# metrics, each providing a unique perspective on spatial alignment:

# - Mean Pearson Correlation: Measures the linear correlation between
#   spatial predictions across treatments.
 
# - Schoener’s D: Quantifies the overlap between two spatial distributions,
#   ranging from 0 (no overlap) to 1 (complete overlap).

# - Hellinger’s I: Provides another measure of overlap or similarity in 
#   probability distributions, sensitive to differences across spatial predictions.
 
# - Spearman Rank Correlation: A non-parametric measure of rank correlation
#   that evaluates the relationship between two distributions.

# These metrics together offer a comprehensive comparison framework for 
# analyzing spatial patterns between SDM treatments.


### 0. LOAD PACKAGES ---------------------------------------------------------------

# Define required packages
required_packages <- c(
  "tidyverse", "rgdal", "sf", "sp", "rnaturalearth", 
  "rnaturalearthdata", "raster", "maptools", "Hmisc",
  "gdalUtilities", "tidyr", "ggspatial", "purrr",
  "ggpubr", "readxl", "devtools", "xlsx", "ENMTools",
  "heatmaply", "isocat", "terra", "Matrix", "corrplot"
)

# Install and load packages if not already loaded
invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# Check working directory
getwd()


### 1. LOAD DATA -------------------------------------------------------

# Load the list of study species
spp_keys <- read_excel("data_inputs/spp_keys.xlsx")

# Define incomplete species to remove
incomplete_species <- c(
  "Sambucusnigr", "Ilexaquifoli", "Hederahelix", "Cornussangui", 
  "Juniperuscom", "Vacciniummyr", "Morusnigra", "Crataegusmon", 
  "Prunusspinos", "Rosacanina", "Solanumnigru", "Taxusbaccata", 
  "Vitisvinifer"
)

# Filter out incomplete species
spp_keys <- subset(spp_keys, !plant_sp %in% incomplete_species)

# Define treatment combinations
comb <- as.character(1:12)


### 2. MEAN PEARSON CORRELATION --------------------------

# Load raster stack of predictions for the first species
Rstack <- stack(list.files(
  path = "data_inputs/sdm_prediction/Arbutusunedo",
  pattern = 'tif', full.names = TRUE
))

# Plot the raster stack to visualize layers
plot(Rstack)

# Calculate Pearson correlation coefficient between layers
Rstats <- layerStats(Rstack, stat = "pearson", asSample = TRUE, na.rm = TRUE)
correlation.x <- Rstats$`pearson correlation coefficient`


# Run a for loop
for(j in 2:nrow(spp_keys)){
  
  # Open predictions for further species
  Rstack <-stack(list.files(path=paste0("data_inputs/sdm_prediction/",
                                       spp_keys[j,1]),
                            pattern='tif',full.names=TRUE))
  
  # Calculate the Pearson correlation
  Rstats <- layerStats(Rstack, "pearson", asSample=TRUE, na.rm=TRUE)
  correlation_2 <- Rstats$`pearson correlation coefficient`
  correlation.x <- (correlation.x + correlation_2)
}

# Obtain the Mean Pearson correlation
correlation_v2 <- correlation.x / 18 # 18 = N of species with all combinations

# Obtain the Mean Pearson correlation matrix and plot it
r <- cor(correlation_v2)
mtcars.rcorr <- rcorr(as.matrix(correlation_v2))
p <- mtcars.rcorr$P
heatmaply_cor(
  r,
  node_type = "scatter",
  point_size_mat = -log10(p), 
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation"),
  symm = FALSE,
  col = viridis(n = 5),
  dendrogram = "both"
)


### 3. SCHOENER’S D, HELLINGER’S I, AND SPEARMAN RANK CORRELATION ---------------------------------------------------------------
 
# Generate empty dataframe
df <- data.frame(D = numeric(), I = numeric(), R = numeric(), 
                 sp = character(), combbase = character(),
                 combcompar = character())

# Run a for loop
for(j in 1:nrow(spp_keys)) {
  for (i in 1:12) {
    sdm_base <- raster(paste("data_inputs/sdm_prediction/",
                             spp_keys[j,1],"/proj_mod_",
                             comb[i],"EM_noTH.tif",sep = ""))
    sdm_base <- rast(sdm_base)
    
    # Compare with all the treatment combinations of the species
    for(s in 1:12) {
      sdm_compare <- raster(paste("data_inputs/sdm_prediction/",
                                  spp_keys[j,1],"/proj_mod_",
                                  comb[s],"EM_noTH.tif",sep = ""))
      sdm_compare <- rast(sdm_compare)
      
      # Get index values
      index_values <- ENMTools::raster.overlap(sdm_base,
                                               sdm_compare)
      D <- index_values$D # D index
      I <- index_values$I # I index
      R <- index_values$rank.cor # Spearman rank correlation
      
      # Add information on treatment combination and species
      combbase <- i
      combcompar <- s
      species <- paste(spp_keys[j,1])
      add <- data.frame(D, I, R, species, combbase, combcompar)
      df <- rbind(df, add)
    }
  }
}

# Aggregate by columns
D_agg <- aggregate(df$D, by = list(df$combcompar, df$combbase),
                   FUN = mean)
I_agg <- aggregate(df$I, by = list(df$combcompar, df$combbase),
                   FUN = mean)
R_agg <- aggregate(df$R, by = list(df$combcompar, df$combbase),
                   FUN = mean)

### Build the correlation matrix for the D index
inds <- cbind(as.integer(D_agg$Group.1), as.integer(D_agg$Group.2))
inds <- t(apply(inds, 1, sort))
res <- sparseMatrix(i = inds[,1], 
                    j = inds[,2], 
                    x = D_agg$x,
                    symmetric = TRUE)
res <- as.matrix(res)

# The values (excluding the diagonal) are +1 inflated. Correct this
res <- res - 1

# Transform zeros to ones
res_def <- res + diag(1, nrow(res))

# Save matrix for the D index
res_def1 <- round(as.data.frame(res_def), 2)

### Build the correlation matrix for the I index
inds <- cbind(as.integer(I_agg$Group.1), as.integer(I_agg$Group.2))
inds <- t(apply(inds, 1, sort))
res <- sparseMatrix(i = inds[,1], 
                    j = inds[,2], 
                    x = I_agg$x,
                    symmetric = TRUE)
res <- as.matrix(res)

# The values (excluding the diagonal) are +1 inflated. Correct this
res <- res - 1

# Transform zeros to ones
res_def <- res + diag(1, nrow(res))

# Save matrix for the I index
res_def2 <- round(as.data.frame(res_def), 2)

### The same for the R index
inds <- cbind(as.integer(R_agg$Group.1), as.integer(R_agg$Group.2))
inds <- t(apply(inds, 1, sort))
res <- sparseMatrix(i = inds[,1], 
                    j = inds[,2], 
                    x = R_agg$x,
                    symmetric = TRUE)
res <- as.matrix(res)

# The values (excluding the diagonal) are +1 inflated. Correct this
res <- res - 1

# Transform zero to one(s)
res_def <- res + diag(1, nrow(res))

# Save matrix for the R index
res_def3 <- round(as.data.frame(res_def), 2)
