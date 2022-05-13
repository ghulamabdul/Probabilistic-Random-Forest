# RAMEFI

This repository provides up-to-date code and data accompanying the paper

> Eisenstein, L., Schulz, B., Qadir, G. A., Pinto, J. G. and Knippertz, P. (2022). 
> Objective identification of high-wind features within extratropical cyclones using a probabilistic random forest (RAMEFI). Part I: Method and illustrative case studies.
> Weather Clim. Dynam. Discuss. \[preprint], https://doi.org/10.5194/wcd-2022-29, in review, 2022.

In particular, code for the implementation of the RAMEFI method, the data that was used in the study and videos of the predictions are available. The data, code and videos at the time of the submission are available at Zenodo (click here for [code and data](https://doi.org/10.5281/zenodo.6541303) and here for the [videos](https://doi.org/10.5281/zenodo.6541277)).


## Code

### Interactive Labelling Tool

We developed an interactive visualisation using the open-source data visualisation package bokeh for Python (https://bokeh.org). The tool allows us to explore several parameters and to select an area by mouse to set labels representing the mesoscale high-wind features, which can then be saved in new files. We supply the code for the tool and an example data set with all needed parameters[^1].

### Random Forests

We applied the random forests using the statistical software R (https://www.r-project.org/). We supply code for the preprocessing of the raw data, the application of the random forest method, evalaution and the postprocessing of the forecasts, which can be found in the `random_forests/` subfolder. 

| File | Description |
| ---- | ----------- | 
| `station_preprocessing` | Preprocessing of NCDF files for station data (transformation to R data). |
| `station_rf_total` | Generation of a random forest based on all case studies. |
| `station_rf_cv` | Application of the random forests for the station data in the cross-validation setting (including an evaluation of forecast performance). |
| `station_pdp` | Generation of the PDP objects for the station-based random forests. |
| `station_postprocessing` | Transformation of the station-based forecasts to csv files. |
| `station_kriging` | Application of the Kriging method on the staion-based random forest predictions. |
| `cosmo_preprocessing` | Preprocessing of NCDF files for COSMO-REA data (transformation to R data). |
| `cosmo_rf_prediction` | Application of the station-based random forests on the COSMO-REA data in the cross-validation setting (including a comparison to climatology). |
| `cosmo_postprocessing` | Transformation of the COSMO-REA-based forecasts to csv files. |

### Figures from the paper

Code for the generation of the figures in the paper can be found in the subfolder `plotting_scripts/`.

| File | Description |
| ---- | ----------- | 
| `figures_paper` | R code for the generation of the reliability diagrams and the predictor importance plots. |
| `plots_paper` | Python code for the generation of the storm and probability maps. |


## Data

Here, we provide a short overview of the data accompanying the study. For more information, we refer to the publication.

### Surface observations

The surface observations were supplied by the German Weather Service (DWD, [Deutscher Wetterdienst](https://www.dwd.de/EN/)) and are proprietary. The wind gust observations were obtained from the Climate Data Center (CDC, [link](https://www.dwd.de/EN/climate_environment/cdc/cdc_node_en.html)) of the DWD. 
Instead of the original data that was used to train the random forests, we supply a simulated data set sampled from the original data such that the interested reader can apply the code using data of the desired format[^1].  All other data such as the (Kriging) predictions and the random forest objects are based on the original data and can be used.

[^1]: Note that the simulated data set is not physically consistent and the predictions will be of poor quality.

The following data is provided in this repository (and/or Zenodo)

| File | Description |
| ---- | ----------- | 
| `example_obs` | Fake surface observations with random values for the interactive labelling tool. |
| `station_fake_data` | Fake surface observations and feature labels generated based on the original data. |
| `station_data` | Original surface observations and feature labels (NOT SUPPLIED). |
| `station_preds` | Random forest predictions, permuted predictions for feature importances and scores (in R and csv format). |
| `station_preds_for_kriging` | Random forest predictions preprocessed for Kriging (csv format). |
| `station_pdp` | PDP curves for the individual storms and accumulated. |
| `station_fi` | Brier score permutation importances calculated for the predictor variables. |
| `random_forest` | Random forest object trained using all storm cases (**ONLY ON ZENODO**). |
| <code>rf/<em>name</em></code> | Random forest object generated in the (cross-validation) fold of *name* (used for prediction of *name*) (**ONLY ON ZENODO**). |
| <code>kriging/<em>name</em>_<em>t</em></code> | Kriging predictions for *name* at time *t*. |
| <code>pdp/<em>name</em></code> | PDP objects corresponding to the random forest object generated in the (cross-validation) fold of *name*. |

### COSMO-REA6

The COSMO-REA6 data is available under https://reanalysis.meteo.uni-bonn.de (Hans-Ertel-Centre for Weather Research, 2019).

The following data is provided in this repository

| File | Description |
| ---- | ----------- | 
| `cosmo_clim` | Brier scores of the random forest and climatological predictions for the individual storms. |
| `cosmo_preds` | Random forest predictions (in R and csv format). |
| <code>cosmo/cosmo_data_<em>name</em></code> | Predictor variables and feature labels for *name* (used for prediction of *name*). |
| <code>cosmo/cosmo_preds_<em>name</em></code> | Predictions for the (cross-validation) fold of *name* generated with the corresponding station-based random forest obejct (in R and csv format). |

## Videos

Videos showing the temporal evolution of the Kriging predictions for the wind features can be found in the subfolder `videos/`.

##
