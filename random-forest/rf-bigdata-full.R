#' # Summarising big data
#' 
#' Having extracted a huge number of variables, let's find out what we got...
#' 
#+ data_setup

output.filename.base <- '../../output/rf-bigdata-try4'

data.filename.big <- '../../data/cohort-datadriven-02.csv'
model.type <- 'rfsrc'
n.data <- NA

surv.predict.old <- c('age', 'smokstatus', 'imd_score', 'gender')
untransformed.vars <- c('time_death', 'endpoint_death', 'exclude')
exclude.vars <-
  c(
    # Entity type 4 is smoking status, which we already have
    "clinical.values.4_data1", "clinical.values.4_data5",
    "clinical.values.4_data6",
    # Entity 13 data2 is the patient's weight centile, and not a single one is
    # entered, but they come out as 0 so the algorithm, looking for NAs, thinks
    # it's a useful column
    "clinical.values.13_data2"
  )

source('../lib/shared.R')

COHORT <- fread(data.filename.big)

percentMissing <- function(x) {
  sum(is.na(x))/length(x) * 100
}

missingness <- sapply(COHORT, percentMissing)

bigdata.prefixes <-
  c(
    'hes.icd.',
    'hes.opcs.',
    'tests.enttype.',
    'clinical.history.',
    'clinical.values.',
    'bnf.'
  )

bigdata.columns <-
  which(
    # Does is start with one of the data column names?
    startsWithAny(names(COHORT), bigdata.prefixes) &
      # And it's not one of the columns we want to exclude?
      !(names(COHORT) %in% exclude.vars)
  )

top.bigdata <- sort(missingness[bigdata.columns])[1:100]

COHORT.bigdata <-
  COHORT[, c(
    untransformed.vars, surv.predict.old, names(top.bigdata)
    ),
    with = FALSE
  ]

# Deal appropriately with missing data
# Most of the variables are number of days since the first record of that type
time.based.vars <-
  names(COHORT.bigdata)[
    startsWithAny(
      names(COHORT.bigdata),
      c('hes.icd.', 'hes.opcs.', 'clinical.history.')
    )
  ]
# Missing values can therefore be replaced with -1, because, if this is to be
# detected, it will be in the future
for (j in time.based.vars) {
  set(COHORT.bigdata, j = j, value = NA2val(COHORT.bigdata[[j]], -1))
}

# The drug data are number of packs prescribed recently, so missing data can be
# replaced with 0, because none were prescribed
prescriptions.vars <- names(COHORT.bigdata)[startsWith(names(COHORT.bigdata), 'bnf.')]
for (j in prescriptions.vars) {
  set(COHORT.bigdata, j = j, value = NA2val(COHORT.bigdata[[j]], 0))
}

# This leaves tests and clinical.values, which are test results and should be
# imputed.

# Manually fix clinical values items...
#
# "clinical.values.1_data1"  "clinical.values.1_data2"
# These are just blood pressure values...fine to impute
#
# "clinical.values.13_data1" "clinical.values.13_data3"
# These are weight and BMI...also fine to impute
#
# Entity 5 is alcohol consumption status, 1 = Yes, 2 = No, 3 = Ex, so should be
# a factor, and NA can be a factor level
COHORT.bigdata$clinical.values.5_data1 <-
  factorNAfix(factor(COHORT.bigdata$clinical.values.5_data1), NAval = 'missing')

# Both gender and smokstatus are factors...fix that
COHORT.bigdata$gender <- factor(COHORT.bigdata$gender)
COHORT.bigdata$smokstatus <-
  factorNAfix(factor(COHORT.bigdata$smokstatus), NAval = 'missing')

# Exclude invalid patients
COHORT.bigdata <- COHORT.bigdata[!COHORT.bigdata$exclude]
COHORT.bigdata$exclude <- NULL

COHORT.bigdata <-
  prepSurvCol(data.frame(COHORT.bigdata), 'time_death', 'endpoint_death', 'Death')

test.set <- testSetIndices(COHORT.bigdata, random.seed = 78361)

surv.predict <- c(surv.predict.old, names(top.bigdata))


source('../lib/rfsrc-cv-mtry-nsplit-logical.R')


#' # Results
#' 
#' 
#' ## Performance
#' 
#' ### Discrimination
#' 
#' C-index is **`r round(surv.model.fit.coeffs['c.test', 'val'], 3)` +/-
#' `r round(surv.model.fit.coeffs['c.test', 'err'], 3)`** on the held-out test
#' set.
#' 
#' ### Calibration
#' 
#' Does the model predict realistic probabilities of an event?
#' 
#+ calibration_plot

calibration.table <-
  calibrationTable(
    # Standard calibration options
    surv.model.fit, COHORT.prep[test.set,],
    # Always need to specify NA imputation for rfsrc
    na.action = 'na.impute'
  )

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' The area between the calibration curve and the diagonal is 
#' **`r round(calibration.score['area'], 3)`** +/-
#' **`r round(calibration.score['se'], 3)`**.
#' 
#' ## Model fit
#' 
#+ resulting_fit

print(surv.model.fit)

#' ## Variable importance
#' 

print(importance(surv.model.fit))