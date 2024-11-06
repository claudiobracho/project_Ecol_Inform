# In this section we load all necessary packages required for the subsequent
# scripts.

# The renv package is used here to ensure reproducibility of the project
# environment. By restoring dependencies from a lockfile, renv can precisely 
# replicate the package versions originally used, which is crucial for consistency 
# across different systems working on this project.

# Note: To avoid compatibility issues, it is recommended to use the original 
# version of R specified in the lockfile when running this code, 
# as newer versions may introduce changes that could affect functionality or 
# package compatibility.

# Ensure the renv package is installed for environment management
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
library(renv)

# Restore project dependencies from the lockfile for reproducibility
renv::restore()

# Define a list of necessary packages
packages <- c(
  "tidyverse", "rgdal", "sf", "sp", "rnaturalearth", "rnaturalearthdata",
  "raster", "maptools", "gdalUtilities", "tidyr", "ggspatial", "purrr",
  "ggpubr", "readxl", "dplyr", "ggplot2", "lsmeans", "stringr", "multcomp",
  "glmmTMB", "ggrepel", "performance", "emmeans", "effects", "DHARMa", 
  "gridExtra", "ggeffects", "sjPlot", "sjmisc", "Hmisc", "devtools", "xlsx",
  "ENMTools", "heatmaply", "isocat", "terra", "Matrix", "corrplot", 
  "biomod2", "dismo", "rstudioapi", "rgeos", "RColorBrewer"
)

# Function to safely load or install a package if missing
load_or_install <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) install.packages(package)
  library(package, character.only = TRUE)
}

# Load or install required packages
invisible(lapply(packages, load_or_install))

# Load additional packages individually in case it is needed to address
# specific issues
individual_packages <- c("heatmaply", "terra", "Matrix", "Hmisc", "readxl",
                         "raster", "sf", "plotly", "ecospat", "TMB", "glmmTMB")
invisible(lapply(individual_packages, load_or_install))

# Take a snapshot of the current environment to save exact package versions
renv::settings$snapshot.type("all")
renv::snapshot()
