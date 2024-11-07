---
title: "Read me"
author: "F Tavernier, V. Segura"
---

Corresponding authors: <flora-tavernier@hotmail.fr>, <vincent.segura@inrae.fr>

This document describes the scripts and datasets used to produce the results presented in our article entitled "Near infrared real time non destructive monitoring of sugar accumulation reveals that single berries ripen two fold faster than previously documented on standard asynchronous samples".

# Repository Structure
## `data/`: raw and processed datasets
Datasets are available at <https://doi.org/10.57745/YGKPZA>

- `NIRSBaies_2021et2022_HPLC.csv`: This table corresponds to the HPLC measurements and NIRS from 2021 and 2022.

- `Training_dataset_2021-2022_CalibNIRS.csv`: This table was used to train the PLSR models for predicting the traits.

- `Validation_dataset_2021-2022_CalibNIRS.csv`: This table was used to validate the models.

- `Suivi_2021-2022.csv`: This table contains the monitoring data from 2021 and 2022.

- `Env_tmps_acc_sucres_baies.RData`: This file contains the PLSDA model, the monitoring list with preprocessed spectra (2021 and 2022) and the traits predictions for the monitoring data.

- `Data_for_calcul_time.csv`: This table was used for calculating the time differences. It contains only the monitoring of 2022 and without any notes. The sugar_cal values were also limited from 0 to 1500 mM (values under 0 were replaced by 0).

- `Data_for_calcul_time_rsquare_poly2.csv`: This table contains the same data as `Data_for_calcul_time.csv` but with also the coefficients and the R² values for the polynomial fit.

- `coef_logi_rsq.csv`: Tis table contains the coefficients and the R² values calculated for the sigmoid fit (logistic regression).

## `scripts/`: analysis scripts and result generation (with R)
These scripts were made in R Studio v 28.3.1 and running R v. 4.3.2 (R Core Team (2023). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>).

- `calibSPIR.R`: This little script contains to all the functions used to build the PLSR model.

- `Prediction_model_preparation.Rmd`: This script generate the training and validation sets using `NIRSBaies_2021et2022_HPLC.RData`, which contains `NIRSBaies_2021et2022_HPLC.csv`. It also create the PLSR models for each trait using `Valid-Train_Dataset_CalibNIRS_2021-2022.Rdata` (`Training_dataset_2021-2022_CalibNIRS.csv` + `Validation_dataset_2021-2022_CalibNIRS.csv`), `outliers.RData` created when the training and validation sets are generated and `calibSPIR.R`.

- `HPLC_traits_PLSDAmodel_and_predict_monit_berries.Rmd`: This script explore `NIRSBaies_2021et2022_HPLC_OK.RData` (univariate analysis, PCA). It builds the PLSDA model for the classification of developmental stages. It generates the PLSR model figures for each traits using `Valid-Train_Dataset_CalibNIRS_2021-2022.Rdata`, `output_cal_pred_MicroNIR_KS_noout.Rdata` (model and outliers generated in `Prediction_model_preparation.Rmd`) and `calibSPIR.R`.

- `Estimation_sugar_acc_time_berries.Rmd`: This script use the PLSDA model and the sugar_cal data saved from `HPLC_traits_PLSDAmodel_and_predict_monit_berries.Rmd` to fit a sigmoid curve (logistic regression) on each monitored berry's kinetic. It generates `Data_for_calcul_time.csv`, `Data_for_calcul_time_rsquare_poly2.csv` and `coef_logi_rsq.csv`. It also fit a 2nd order polynomial, it was done to test another method than the sigmoid.

## `notebooks/`: Google Colab notebook used for the estimation of time differences (with Python)
This notebook was made using Google Colab (Python v. 3.11)

- `Timediff_9.ipynb`: This notebook generates the polynomial and the sigmoid fit plots for each berry and calculate the sugar accumulation time differences between each variety in 2022.

# Utilisation
To reproduce our results, you need to first run `Prediction_model_preparation.Rmd` in order to generate the PLSR models outputs. 
Then, run `HPLC_traits_PLSDAmodel_and_predict_monit_berries.Rmd` to generate the PLSDA model and the sugar_cal data. 

Finally, run `Estimation_sugar_acc_time_berries.Rmd` to generate the 3 tables used in the Google Colab notebook `Timediff_9.ipynb` (Python).
