# In this section, we execute all Generalized Linear Mixed Models (GLMMs) 
# included in our analysis. Each GLMM uses a different performance 
# metric as the response variable, allowing us to assess model performance 
# across various criteria.
 
# The GLMMs specifically evaluate:
 
# - AUC (Area Under the Curve): Measures the model's ability to distinguish 
#   between presence and absence locations, with higher values indicating 
#   better predictive accuracy.
 
# - TSS (True Skill Statistic): Balances sensitivity and specificity in 
#   model predictions, helping assess predictive reliability.
 
# - OUC (Occurrence in Unsampled Cells): Evaluates the model’s extrapolative 
#   power and impact of sampling biases by measuring the probability of
#   occurrence in unsampled cells.

# - Boyce Index: Assesses the model's ability to spatially rank the probability
#   of occurrence of species. Higher values of the Boyce Index indicate that the
#   model effectively ranks suitable areas in alignment with observed species
#   distributions.
 
# - COI (Congruence Index): Quantifies the agreement between model predictions 
#   and an independent, expert-defined distribution, validating alignment 
#   with expert knowledge.
 
# Each metric is treated as the response variable in its respective GLMM, 
# enabling a comprehensive evaluation of model performance across multiple 
# criteria.


### 0. LOAD PACKAGES ---------------------------------------------------------------

# Define required packages
required_packages <- c(
       "tidyr", "dplyr", "ggplot2", "lsmeans","stringr","multcomp",
       "glmmTMB", "ggrepel", "ggpubr", "performance", 
       "emmeans", "effects", "DHARMa", "gridExtra", 
       "purrr", "readxl", "ggeffects", "sjPlot", "sjmisc")

# Install and load packages if not already loaded
invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# Check working directory
getwd()


### 1. CONDUCT GLMMS ---------------------------------------------------------------

# Open dataframe
df<- read.csv("data_inputs/metrics_SDM.csv")

## Conduct GLMM for AUC
model_AUC <- glmmTMB(AUC ~ filter + methodPA + nPA + (1|species), 
                     data = df,
                     family = beta_family(link = "logit"))
# Get its summary
summary(model_AUC)

# Get its confidence intervals
confint(model_AUC)

# Get a post-hoc test regarding filters methods to generate 
# pseudo-absences and numbers of pseudo-absences tested in our 
# GLMMs
lsm <- lsmeans(model_AUC, c("filter"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_AUC, c("methodPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_AUC, c("nPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")

# Include interactive effects of latitude, longitude and the 
# number of cleaned occurrences
model_AUC <- glmmTMB(AUC ~ filter * cleaned + filter * y + 
                       filter * x + methodPA * cleaned  + 
                       methodPA * y + methodPA * x + 
                       nPA * cleaned + nPA * y + nPA * x +
                       (1|species), 
                     data = df,
                     family = beta_family(link = "logit"))

# Get its summary
summary(model_AUC)

# Get its confidence intervals
confint(model_AUC)


## Conduct GLMM for TSS
model_TSS <- glmmTMB(TSS ~ filter + methodPA + nPA + (1|species), 
                     data = df,
                     family = beta_family(link = "logit"))

# Get its summary
summary(model_TSS)

# Get its confidence intervals
confint(model_TSS)

# Get a post-hoc test regarding filters methods to generate 
# pseudo-absences and numbers of pseudo-absences tested in our 
# GLMMs
lsm <- lsmeans(model_TSS, c("filter"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_TSS, c("methodPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_TSS, c("nPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")

# Include interactive effects of latitude, longitude and the 
# number of cleaned occurrences
model_TSS <- glmmTMB(TSS ~ filter * cleaned + filter * y + 
                       filter * x + methodPA * cleaned  + 
                       methodPA * y + methodPA * x + nPA *
                       cleaned + nPA * y + nPA * x + (1|species), 
                     data = df,
                     family = beta_family(link = "logit"))

# Get its summary
summary(model_TSS)

# Get its confidence intervals
confint(model_TSS)


## Conduct GLMM for OUC
model_OUC <- glmmTMB(OUC ~ filter + methodPA + nPA + 
                              (1|species), 
                            data = df,
                            family = beta_family(link = "logit"))

# Get its summary
summary(model_OUC)

# Get its confidence intervals
confint(model_OUC)

# Get a post-hoc test regarding filters methods to generate 
# pseudo-absences and numbers of pseudo-absences tested in our 
# GLMMs
lsm <- lsmeans(model_OUC, c("filter"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_OUC, c("methodPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_OUC, c("nPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")

# Include interactive effects of latitude, longitude and the 
# number of cleaned occurrences
model_OUC <- glmmTMB(OUC ~ filter * cleaned + 
                              filter * y + filter * x +
                              methodPA * cleaned  + methodPA * y + 
                              methodPA * x + nPA * cleaned + 
                              nPA * y + nPA * x + (1|species), 
                            data = df,
                            family = beta_family(link = "logit"))

# Get its summary
summary(model_OUC)

# Get its confidence intervals
confint(model_OUC)


## Conduct GLMM for boyce
model_boyce <- glmmTMB(boyce ~ filter + methodPA + nPA + 
                       (1|species), 
                     data = df,
                     family = beta_family(link = "logit"))

# Get its summary
summary(model_boyce)

# Get its confidence intervals
confint(model_boyce)

# Get a post-hoc test regarding filters methods to generate 
# pseudo-absences and numbers of pseudo-absences tested in our 
# GLMMs
lsm <- lsmeans(model_boyce, c("filter"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_boyce, c("methodPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_boyce, c("nPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")

# Include interactive effects of latitude, longitude and the 
# number of cleaned occurrences
model_boyce <- glmmTMB(boyce ~ filter * cleaned + 
                       filter * y + filter * x +
                       methodPA * cleaned  + methodPA * y + 
                       methodPA * x + nPA * cleaned + 
                       nPA * y + nPA * x + (1|species), 
                     data = df,
                     family = beta_family(link = "logit"))

# Get its summary
summary(model_boyce)

# Get its confidence intervals
confint(model_boyce)



## Conduct GLMM for COI

# Identify incomplete species
incomplete_species <- c(
  "Sambucusnigr", "Ilexaquifoli", "Hederahelix", "Cornussangui", 
  "Juniperuscom", "Vacciniummyr", "Morusnigra", "Crataegusmon", 
  "Prunusspinos", "Rosacanina", "Solanumnigru", "Taxusbaccata", 
  "Vitisvinifer"
)

# Filter out the incomplete species and remove the column 'X'
congruence_index <- subset(df, !species %in% incomplete_species)
congruence_index$X <- NULL

# Then, conduct the GLMM
model_COI <- glmmTMB(COI ~ filter + methodPA + nPA + (1|species), 
                    data = congruence_index,
                    family = beta_family(link = "logit"))

# Get its summary
summary(model_COI)

# Get its confidence intervals
confint(model_COI)

# Get a post-hoc test regarding filters methods to generate 
# pseudo-absences and numbers of pseudo-absences tested in our 
# GLMMs
lsm <- lsmeans(model_COI, c("filter"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_COI, c("methodPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")
lsm <- lsmeans(model_COI, c("nPA"))
lsmeans::contrast(lsm, Letters=letters, method = "pairwise")

# Include interactive effects of latitude, longitude and the 
# number of cleaned occurrences
model_COI <- glmmTMB(COI ~ filter * cleaned + filter * y + 
                      filter * x + methodPA * cleaned  + 
                      methodPA * y + methodPA * x + nPA * cleaned +
                      nPA * y + nPA * x + (1|species), 
                    data = congruence_index,
                    family = beta_family(link = "logit"))

# Get its summary
summary(model_COI)

# Get its confidence intervals
confint(model_COI)


## Run again AUC, TSS, and OUC metrics, but only with the six 
# species available for the COI metric

# First, subset our dataframe accordingly
congruence_index_six <- subset(congruence_index, 
                               congruence_index$COI > 0)

# GLMM for AUC
model_AUC <- glmmTMB(AUC ~ filter * cleaned + filter * y + 
                       filter * x + methodPA * cleaned  + 
                       methodPA * y + methodPA * x + nPA *
                       cleaned + nPA * y + nPA * x + (1|species), 
                     data = congruence_index_six,
                     family = beta_family(link = "logit"))

options(scipen = 999) # to delete scientific notation

# Get its summary
summary(model_AUC)

# Get its confidence intervals
confint(model_AUC)

# GLMM for TSS
model_TSS <- glmmTMB(TSS ~ filter * cleaned + filter * y + 
                       filter * x + methodPA * cleaned  + 
                       methodPA * y + methodPA * x + nPA * cleaned +
                       nPA * y + nPA * x + (1|species), 
                     data = congruence_index_six,
                     family = beta_family(link = "logit"))

# Get its summary
summary(model_TSS)

# Get its confidence intervals
confint(model_TSS)

# GLMM for OUC
model_OUC <- glmmTMB(OUC ~ filter * cleaned + 
                              filter * y + filter * x +
                              methodPA * cleaned  + methodPA * y + 
                              methodPA * x + nPA * cleaned + 
                              nPA * y + nPA * x + (1|species), 
                            data = congruence_index_six,
                            family = beta_family(link = "logit"))

# Get its summary
summary(model_OUC)

# Get its confidence intervals
confint(model_OUC)



# GLMM for boyce
model_boyce <- glmmTMB(boyce ~ filter * cleaned + 
                       filter * y + filter * x +
                       methodPA * cleaned  + methodPA * y + 
                       methodPA * x + nPA * cleaned + 
                       nPA * y + nPA * x + (1|species), 
                     data = congruence_index_six,
                     family = beta_family(link = "logit"))

# Get its summary
summary(model_boyce)

# Get its confidence intervals
confint(model_boyce)



## Get mean AUC, TSS, OUC, and COI per plant sp.
df_mean_global <- df %>%
  group_by(species) %>%
  summarise(AUC=mean(AUC,na.rm=T),TSS=mean(TSS,na.rm=T),
            OUC=mean(OUC,na.rm=T), COI=mean(COI,na.rm=T))


### 2. PLOT EFFECTS ---------------------------------------------------------------

# Change names to adjust sizes in the final plot
df$nPA <- gsub("Npresences", "nP", df$nPA)
df <- rename(df, numberPA = nPA)

# First example with the AUC model that includes interactive effects
plot(ggpredict(model_AUC, terms = c("filter", "methodPA", 
                                    "nPA")))

# Plot interactions per each explanatory variable
plot_model(model_AUC, type = "int")

# Change plot theme aesthetics
set_theme(
  geom.outline.color = "antiquewhite4", 
  geom.outline.size = 1, 
  geom.label.size = 2,
  geom.label.color = "grey50",
  title.color = "black", 
  title.size = 1.5, 
  axis.angle.x = 45, 
  axis.textcolor = "black", 
  legend.size = 1,
  legend.title.size = 1,
  legend.title.color = "black",
  legend.color = "black",
  base = theme_bw()
)

# Plot interactions per each explanatory variable and metric, in
# order to obtain a grid plot
model_AUC <- glmmTMB(AUC ~ filter * cleaned + filter * y + 
                       filter * x + methodPA * cleaned  + 
                       methodPA * y + methodPA * x + numberPA *
                       cleaned + numberPA * y + numberPA * x + 
                       (1|species), 
                     data = df,
                     family = beta_family(link = "logit"))

# Define parameters for each plot
plot_params <- list(
  list(terms = c("y", "methodPA"), axis_title = c("", "AUC")),
  list(terms = c("x", "methodPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "methodPA"), axis_title = c("", "")),
  list(terms = c("y", "numberPA"), axis_title = c("", "AUC")),
  list(terms = c("x", "numberPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "numberPA"), axis_title = c("", "")),
  list(terms = c("y", "filter"), axis_title = c("(scaled) Latitude", "AUC")),
  list(terms = c("x", "filter"), axis_title = c("(scaled) Longitude", "")),
  list(terms = c("cleaned", "filter"), axis_title = c("(scaled) N Ps", ""))
)

# Generate plots in a loop and store in a list
plots <- lapply(plot_params, function(params) {
  plot_model(
    model_AUC, type = "pred", terms = params$terms,
    show.values = TRUE, show.p = TRUE, show.legend = FALSE,
    axis.title = params$axis_title, title = "", axis.lim = c(0.70, 1)
  )
})
plots

# Arrange previous objects
iauc <- grid.arrange(grobs = plots, ncol = 3)

# save plot
~#ggsave("figures/Figure_IAUC.jpeg", width = 22, height = 18, 
#       units = "cm", iauc)


# Same with TSS
model_TSS <- glmmTMB(TSS ~ filter * cleaned + filter * y + 
                       filter * x + methodPA * cleaned  + 
                       methodPA * y + methodPA * x + 
                       numberPA * cleaned + numberPA * y + 
                       numberPA * x + (1|species), 
                     data = df,
                     family = beta_family(link = "logit"))


# Define parameters for each plot
plot_params <- list(
  list(terms = c("y", "methodPA"), axis_title = c("", "TSS")),
  list(terms = c("x", "methodPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "methodPA"), axis_title = c("", "")),
  list(terms = c("y", "numberPA"), axis_title = c("", "TSS")),
  list(terms = c("x", "numberPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "numberPA"), axis_title = c("", "")),
  list(terms = c("y", "filter"), axis_title = c("(scaled) Latitude", "TSS")),
  list(terms = c("x", "filter"), axis_title = c("(scaled) Longitude", "")),
  list(terms = c("cleaned", "filter"), axis_title = c("(scaled) N Ps", ""))
)

# Generate plots in a loop and store in a list
plots <- lapply(plot_params, function(params) {
  plot_model(
    model_TSS, type = "pred", terms = params$terms,
    show.values = TRUE, show.p = TRUE, show.legend = FALSE,
    axis.title = params$axis_title, title = "", axis.lim = c(0.40, 1)
  )
})
plots

# Arrange previous objects
iTSS <- grid.arrange(grobs = plots, ncol = 3)

# save plot
#ggsave("figures/Figure_ITSS.jpeg", width = 22, height = 18, 
#       units = "cm", itss)


# Same with OUC
model_OUC <- glmmTMB(OUC ~ filter * cleaned + 
                              filter * y + filter * x + 
                              methodPA * cleaned  + methodPA * y + 
                              methodPA * x + numberPA * cleaned +
                              numberPA * y + numberPA * x + 
                              (1|species), 
                            data = df,
                            family = beta_family(link = "logit"))


# Define parameters for each plot
plot_params <- list(
  list(terms = c("y", "methodPA"), axis_title = c("", "OUC")),
  list(terms = c("x", "methodPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "methodPA"), axis_title = c("", "")),
  list(terms = c("y", "numberPA"), axis_title = c("", "OUC")),
  list(terms = c("x", "numberPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "numberPA"), axis_title = c("", "")),
  list(terms = c("y", "filter"), axis_title = c("(scaled) Latitude", "OUC")),
  list(terms = c("x", "filter"), axis_title = c("(scaled) Longitude", "")),
  list(terms = c("cleaned", "filter"), axis_title = c("(scaled) N Ps", ""))
)

# Generate plots in a loop and store in a list
plots <- lapply(plot_params, function(params) {
  plot_model(
    model_OUC, type = "pred", terms = params$terms,
    show.values = TRUE, show.p = TRUE, show.legend = FALSE,
    axis.title = params$axis_title, title = "", axis.lim = c(0, 1)
  )
})
plots

# Arrange previous objects
iOUC <- grid.arrange(grobs = plots, ncol = 3)

# save plot
#ggsave("figures/Figure_IOUC.jpeg", width = 22, height = 18, 
#       units = "cm", iouc)


# Same with boyce
model_boyce <- glmmTMB(boyce ~ filter * cleaned + 
                       filter * y + filter * x + 
                       methodPA * cleaned  + methodPA * y + 
                       methodPA * x + numberPA * cleaned +
                       numberPA * y + numberPA * x + 
                       (1|species), 
                     data = df,
                     family = beta_family(link = "logit"))

# Define parameters for each plot
plot_params <- list(
  list(terms = c("y", "methodPA"), axis_title = c("", "boyce")),
  list(terms = c("x", "methodPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "methodPA"), axis_title = c("", "")),
  list(terms = c("y", "numberPA"), axis_title = c("", "boyce")),
  list(terms = c("x", "numberPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "numberPA"), axis_title = c("", "")),
  list(terms = c("y", "filter"), axis_title = c("(scaled) Latitude", "boyce")),
  list(terms = c("x", "filter"), axis_title = c("(scaled) Longitude", "")),
  list(terms = c("cleaned", "filter"), axis_title = c("(scaled) N Ps", ""))
)

# Generate plots in a loop and store in a list
plots <- lapply(plot_params, function(params) {
  plot_model(
    model_boyce, type = "pred", terms = params$terms,
    show.values = TRUE, show.p = TRUE, show.legend = FALSE,
    axis.title = params$axis_title, title = "", axis.lim = c(0, 1)
  )
})
plots

# Arrange previous objects
iboyce <- grid.arrange(grobs = plots, ncol = 3)

# save plot
#ggsave("figures/Figure_Iboyce.jpeg", width = 22, height = 18, 
#       units = "cm", iboyce)


# Same with COI

# Identify incomplete species
incomplete_species <- c(
  "Sambucusnigr", "Ilexaquifoli", "Hederahelix", "Cornussangui", 
  "Juniperuscom", "Vacciniummyr", "Morusnigra", "Crataegusmon", 
  "Prunusspinos", "Rosacanina", "Solanumnigru", "Taxusbaccata", 
  "Vitisvinifer"
)

# Filter out the incomplete species and remove the 'X' column
df <- subset(df, !species %in% incomplete_species)

# Get model and individual objects
model_COI <- glmmTMB(COI ~ filter * cleaned + filter * y + 
                      filter * x + methodPA * cleaned  + 
                      methodPA * y + methodPA * x + 
                      numberPA * cleaned + numberPA * y + 
                      numberPA * x + (1|species), 
                    data = df,
                    family = beta_family(link = "logit"))

# Define parameters for each plot
plot_params <- list(
  list(terms = c("y", "methodPA"), axis_title = c("", "COI")),
  list(terms = c("x", "methodPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "methodPA"), axis_title = c("", "")),
  list(terms = c("y", "numberPA"), axis_title = c("", "COI")),
  list(terms = c("x", "numberPA"), axis_title = c("", "")),
  list(terms = c("cleaned", "numberPA"), axis_title = c("", "")),
  list(terms = c("y", "filter"), axis_title = c("(scaled) Latitude", "COI")),
  list(terms = c("x", "filter"), axis_title = c("(scaled) Longitude", "")),
  list(terms = c("cleaned", "filter"), axis_title = c("(scaled) N Ps", ""))
)

# Generate plots in a loop and store in a list
plots <- lapply(plot_params, function(params) {
  plot_model(
    model_COI, type = "pred", terms = params$terms,
    show.values = TRUE, show.p = TRUE, show.legend = FALSE,
    axis.title = params$axis_title, title = "", axis.lim = c(0.50, 1)
  )
})
plots

# Arrange previous objects
iCOI <- grid.arrange(grobs = plots, ncol = 3)

# save plot
#ggsave("figures/Figure_ICOI.jpeg", width = 22, height = 18, 
#       units = "cm", iCOI)

