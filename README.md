# Machine learning on electronic health records

Repository of code developed for *Machine learning models in electronic health records can outperform conventional survival models for predicting patient mortality in coronary artery disease*.

## Introduction

The R scripts provided in this repository were developed to perform survival modelling on 100,000 patients’ electronic health records.

## Usage

To use these scripts, download them into a subfolder ``code`` of a folder which also contains subfolders ``data`` (which should contain input data) and ``output`` (where output files are placed by the scripts).

Most scripts [can use ``spin`` to generate reports](http://deanattali.com/2015/03/24/knitrs-best-hidden-gem-spin/), and the simplest way to do this is to execute them using the ``spinme.R`` script provided. For example, to run the discretised Cox model script and generate a report, navigate to the ``code/cox-ph`` folder of your project, and run

```
Rscript ../scripts/spinme.R cox-discretised.R
```

This will generate a report ``cox-discretised.html`` in the working directory, as well as additional files in the ``output`` folder.

## Project inventory

This is a brief description of the important files in this project, broken down by folder.

### age-only

#### age-only.R

This calculates the C-index and calibration score for a model where age is the only variable. For the C-index, it is assumed that the older patient will die first, and for calibration the Kaplan–Meier estimator for patients of a given age is assumed to be the risk estimate for all patients of that age.

### cox-ph

Various Cox proportional hazards models. Those prefixed ``caliber-replicate-`` are based on the Cox model developed in [Rapsomaniki et al. 2014](https://academic.oup.com/eurheartj/article-lookup/doi/10.1093/eurheartj/eht533) (DOI: [10.1093/eurheartj/eht533](https://dx.doi.org/10.1093/eurheartj/eht533)).

#### caliber-replicate-with-imputation.R

This model is as close to identical to Rapsomaniki et al. as possible, using a five-fold multiply-imputed dataset, with continuous variables scaled as in that paper.

#### caliber-replicate-with-missing.R

This model uses the same scaling as Rapsomaniki et al., but is conducted on a single dataset with missing values represented by missing indicator variables rather than imputed.

#### caliber-scale.R

Functions to scale data for Cox modelling.

#### cox-discrete-elasticnet.R

Discrete elastic net Cox model for the data-driven modelling, which cross-validates to find the optimal _α_ and then bootstraps to establish distributions for the other fitted parameters.

#### cox-discrete-varsellogrank.R

Discrete Cox model for the data-driven modelling, which cross-validates over number of variables used, drawing from a list ranked by univariate logrank tests.

#### cox-discretised.R

This model uses the expert-selected dataset with discretised versions of continuous variables to allow missing values to be incorporated, and cross-validates to determine the discretisation scheme.

#### cox-discretised-imputed.R

This model uses the imputed version of the expert-selected dataset with discretised versions of continuous variables, following the same method as above.

#### rapsomaniki-cox-values-from-paper.csv

Values for Cox coefficients transcribed from Rapsomaniki et al., used to check for consistency between that model and these.

### lib

Shared libraries for functions and common routines.

#### all-cv-bootstrap.R

Script to cross-validate discretisation schemes, followed by bootstrapping the selected optimal model. Works for both Cox modelling and random forests with either ``randomForestSRC`` or ``ranger``. Discretised random forests were not used in the final analysis, as there was no appreciable performance gain.

#### handy.R

Shortcuts and wrappers, from [Andrew’s](https://github.com/ajsteele/) handy [handy.R](https://github.com/ajsteele/handy.R) script.

#### handymedical.R

Useful functions and wrappers for preparing and manipulating data, and making use of Cox models and random forests with either ``randomForestSRC`` or ``ranger``, including bootstrapping, as transparent and consistent as possible. These functions are hopefully of general use for other survival modelling projects; dataset-specific functions are defined in ``shared.R``.

#### rfsrc-cv-mtry-nsplit-logical.R

Script to cross-validate ``randomForestSRC`` survival forests using the large dataset, optimising the ``mtry`` and ``nsplit`` hyperparameters.

#### rfsrc-cv-nsplit-bootstrap.R

Script to cross-validate ``randomForestSRC`` survival forests using the expert-selected dataset, optimising ``nsplit``.

#### shared.R

This script is run at the start of most model scripts, and defines a random seed, plus variables and functions which will be useful. The data-parsing functions here are specific to the scheme of this particular dataset and so were excluded from ``handymedical.R``.

### overview

Various scripts for exploring the dataset and retrieving and plotting results for publication.

#### all-models.R

Produces a graph of the C-index and calibration score from all models. The basis of Fig. 1 in the paper.

#### bigdata-mtry-nsplit.R

Plots a line graph showing C-index performance of random forests depending on ``mtry`` and ``nsplit`` in the large dataset.

#### calibration-plots.R

Plots two example calibration curves to show how the calibration score is calculated. The basis of Fig. 2 in the paper.

#### cohort-tables.R

Prints a number of summary statistics used for Table 2 in the paper.

#### explore-dataset.R

A number of quick exploratory graphs and comparisons to explore the expert-selected dataset, with a particular focus on degrees and distribution of missing data.

#### missing-values-risk.R

Compares coefficients for Cox models. First, continuous imputed vs continuous with missing indicators and discrete; second, ranges of continuous values’ associated risks with those associated with a value being missing; finally, survival curves for patients with a particular value missing vs present. The basis of Fig. 3 in the paper.

#### performance-differences.R

Pairwise differences with uncertainty in C-index and calibration between all models tested, ascertained by finding the distribution of differences between bootstrap replicates for each model-pair.

#### variable-effects.R

Plots of variable effects for continuous and discrete Cox models, and random forests. The basis of Fig. 4 in the paper.

#### variable-importances.R

Plots permutation variable importances calculated for the final data-driven models, post variable selection. The basis of Fig. 5 in the paper.

### random-forest

#### rf-age.R

Building a random forest with fewer variables (including just age) to experiment with predictive power.

#### rf-classification.R

Classification forest for death at 5 years, in an attempt to improve calibration score of the resulting model.

#### rf-imputed.R

Random forest on the imputed dataset as an empirical test of whether imputation provides an advantage.

#### rfsrc-cv.R

Random forest on the expert-selected dataset, which uses ``rfsrc-cv-nsplit-bootstrap.R`` from ``lib`` (see above) to fit its forest.

#### rf-varsellogrank.R

Random forest for the data-driven modelling, which cross-validates over number of variables used, drawing from a list ranked by univariate logrank tests.

#### rf-varselmiss.R

Random forest for the data-driven modelling, which cross-validates over number of variables used, drawing from a list ranked by decreasing missingness.

#### rf-varselrf-eqv.R

Random forest for the data-driven modelling, which cross-validates over number of variables used, drawing from a list ranked by the variable importance of a large random forest fitted to all the data. Modelled after [varSelRF](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-3).

### scripts

This folder is for a few miscellaneous short scripts.

#### export-bigdata.R

This script was used to export anonymised data for the data-driven modelling with ~600 variables. Its file paths are not correct because they are localised for the secure environment where the raw data are stored.

#### spinme.R

This wrapper makes it easy to spin a script into an HTML report from the command line (see the example command at the top of this readme).

## Notes

This repository has been tidied up so that only scripts relevant to the final publication are preserved. Various initial and exploratory analysis scripts have been removed for clarity. If for any reason these are of interest, they are present in commit 08934808c497a0f094c71a731cb9cb2564e4cc0f, the final commit before the tidy-up began.