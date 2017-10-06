# Machine learning on electronic health records

Repository of code developed for TODO: paper title and reference.

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

### scripts

This folder is for a few short, useful scripts.

#### export-bigdata.R

This script was used to export anonymised data for the data-driven modelling with ~600 variables. Its file paths are not correct because they are localised for the secure environment where the raw data are stored.