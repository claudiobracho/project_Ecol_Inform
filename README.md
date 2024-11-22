# project_Ecol_Inform

# Spatially explicit metrics improve the evaluation of species distribution models facing sampling biases

Authors: Bracho-Estévanez, C. A., Arenas-Castro, S., González-Varo, J. P., & González-Moreno, P.

Citation: Bracho-Estévanez, C. A., Arenas-Castro, S., González-Varo, J. P., & González-Moreno, P. 2024. Spatially explicit metrics improve the evaluation of species distribution models facing sampling biases. Ecological Informatics, XX(X). XXX-XXX.

## Brief description

Here we provide easy and open access both to our code and inputs, which reproduce most relevant results detailed in the main text of our article. Importantly, this includes the script to obtain OUC and COI metrics, which allows to implement them in further works.

## Data structure

The folder "project_Ecol_Inform" includes three folders (code, data_inputs, and renv) and three additional objects (project_Ecol_Inform, README, and renv.lock). Note that renv folder and large data inputs are only available via figshare at <https://doi.org/10.6084/m9.figshare.27623874> .

## 1- code folder

> This folder includes six objects ("code_packages.R", "code_spatial_metrics.R", "code_spatial_comparison.R", "code_GLMMs.R", "code_sdm1_loop.R", "code_sdm2_metrics.R").
>
> -   "code_packages" permits (via the renv package) to refresh or load all necessary packages. Thus, it should be run first.
>
> -   "code_spatial_metrics.R" calculates OUC and CI metrics given SDM predictions (requires spatial outputs).
>   
> -   "code_spatial_comparison.R" calculates the mean Pearson correlation, the Schoener’s D, the Hellinger’s I, and the Spearman rank correlation between SDM predictions (requires spatial outputs).
>
> -   "code_GLMMs.R" runs all generalized mixed models reported in our main text.
>
> -   "code_sdm1_loop.R", uses biomod2 package to generate SDM per species. Requires occurrence and environmental data available via figshare at <https://doi.org/10.6084/m9.figshare.27623874>.
>
> -   "code_sdm2_metrics.R integrates the metrics results of th SDM in a single file for all species

## 2- data_inputs folder

> This folder includes four subfolders ("environmental variables", "cleaned records", "AFE" and "SDM prediction", only available via figshare at <https://doi.org/10.6084/m9.figshare.27623874>), and four additional objects ("metrics_SDM.csv", "samplingeffort_10km.tif", "samplingeffort_10km.xml", and "spp_keys.csv") necessary to run provided codes.

## 3- renv folder

> This folder is automatiquely generated when implementing the renv package. It includes a subfolder ("library") with information on the packages used by our R project.

## 4- Further objects

> project_Ecol_Inform R project renv.lock (here readers can consult the versions of each package loaded in our R project)
