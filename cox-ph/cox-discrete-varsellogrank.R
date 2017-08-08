#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- TRUE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Variable selection in data-driven health records
#' 
#' Having extracted around 600 variables which occur most frequently in patient
#' records, let's try to narrow these down using a methodology based on varSelRf
#' combined with survival modelling. We'll find the predictability of variables
#' as defined by the p-value of a logrank test on survival curves of different
#' categories within that variable, and then iteratively throw out unimportant
#' variables, cross-validating for optimum performance.
#' 
#' ## User variables
#' 
#+ user_variables

output.filename.base <- '../../output/cox-bigdata-varsellogrank-01'

cv.n.folds <- 3
vars.drop.frac <- 0.2 # Fraction of variables to drop at each iteration
bootstraps <- 200

n.data <- NA # This is after any variables being excluded in prep

n.threads <- 20

#' ## Data set-up
#' 
#+ data_setup

data.filename.big <- '../../data/cohort-datadriven-02.csv'

surv.predict.old <- c('age', 'smokstatus', 'imd_score', 'gender')
untransformed.vars <- c('time_death', 'endpoint_death', 'exclude')

source('../lib/shared.R')
require(xtable)

# Define these after shared.R or they will be overwritten!
exclude.vars <-
  c(
    # Entity type 4 is smoking status, which we already have
    "clinical.values.4_data1", "clinical.values.4_data5",
    "clinical.values.4_data6",
    # Entity 13 data2 is the patient's weight centile, and not a single one is
    # entered, but they come out as 0 so the algorithm, looking for NAs, thinks
    # it's a useful column
    "clinical.values.13_data2",
    # Entities 148 and 149 are to do with death certification. I'm not sure how 
    # it made it into the dataset, but since all the datapoints in this are
    # looking back in time, they're all NA. This causes rfsrc to fail.
    "clinical.values.148_data1", "clinical.values.148_data2",
    "clinical.values.148_data3", "clinical.values.148_data4",
    "clinical.values.148_data5",
    "clinical.values.149_data1", "clinical.values.149_data2"
  )

COHORT <- fread(data.filename.big)

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
  colnames(COHORT)[
    which(
      # Does is start with one of the data column names?
      startsWithAny(names(COHORT), bigdata.prefixes) &
        # And it's not one of the columns we want to exclude?
        !(colnames(COHORT) %in% exclude.vars)
    )
    ]

COHORT.bigdata <-
  COHORT[, c(
    untransformed.vars, surv.predict.old, bigdata.columns
  ),
  with = FALSE
  ]

# Get the missingness before we start removing missing values
missingness <- sort(sapply(COHORT.bigdata, percentMissing))
# Remove values for the 'untransformed.vars' above, which are the survival
# values plus exclude column
missingness <- missingness[!(names(missingness) %in% untransformed.vars)]

# Deal appropriately with missing data
# Most of the variables are number of days since the first record of that type
time.based.vars <-
  names(COHORT.bigdata)[
    startsWithAny(
      names(COHORT.bigdata),
      c('hes.icd.', 'hes.opcs.', 'clinical.history.')
    )
    ]
# We're dealing with this as a logical, so we want non-NA values to be TRUE,
# is there is something in the history
for (j in time.based.vars) {
  set(COHORT.bigdata, j = j, value = !is.na(COHORT.bigdata[[j]]))
}

# Again, taking this as a logical, set any non-NA value to TRUE.
prescriptions.vars <- names(COHORT.bigdata)[startsWith(names(COHORT.bigdata), 'bnf.')]
for (j in prescriptions.vars) {
  set(COHORT.bigdata, j = j, value = !is.na(COHORT.bigdata[[j]]))
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

# If n.data was specified, trim the data table down to size
if(!is.na(n.data)) {
  COHORT.bigdata <- sample.df(COHORT.bigdata, n.data)
}

# Define test set
test.set <- testSetIndices(COHORT.bigdata, random.seed = 78361)

# Start by predicting survival with all the variables provided
surv.predict <- c(surv.predict.old, bigdata.columns)

# Set up a csv file to store calibration data, or retrieve previous data
calibration.filename <- paste0(output.filename.base, '-varselcalibration.csv')

varLogrankTest <- function(df, var) {
  # If there's only one category, this is a single-valued variable so you can't
  # do a logrank test on different values of it...
  if(length(unique(NArm(df[, var]))) == 1) {
    return(NA)
  }
  
  # If it's a logical, make an extra column for consistency of later code
  if(class(df[, var]) == 'logical') {
    df$groups <- factor(ifelse(df[, var], 'A', 'B'))
  # If it's numeric, split it into four quartiles
  } else if(class(df[, var]) == 'numeric') {
    # First, discard all rows where the value is missing
    df <- df[!is.na(df[, var]), ]
    # Then, assign quartiles
    df$groups <- 
      factor(
        findInterval(
          df[, var],
          quantile(df[, var], probs=c(0, 0.25, .5, .75, 1))
        )
      )
    
  } else {
    # Otherwise, it's a factor, so leave it as-is
    df$groups <- df[, var]
  }
  
  # Perform a logrank test on the data
  lr.test <- 
    survdiff(
      as.formula(paste0('Surv(surv_time, surv_event) ~ groups')),
      df
    )
  # Return the p-value of the logrank test
  pchisq(lr.test$chisq, length(lr.test$n)-1, lower.tail = FALSE)
}

# Don't use the output variables in our list
vars.to.check <-
  names(COHORT.bigdata)[!(names(COHORT.bigdata) %in% c('surv_time', 'surv_event'))]

var.logrank.p <-
  sapply(
    X = vars.to.check, FUN = varLogrankTest,
    df = COHORT.bigdata[-test.set, ]
  )

# Sort them, in ascending order because small p-values indicate differing
# survival curves
var.logrank.p <- sort(var.logrank.p, na.last = TRUE)

# Create process settings

# Variables to leave alone, including those whose logrank p-value is NA because
# that means there is only one value in the column and so it can't be discretised
# properly anyway
vars.noprocess <- c('surv_time', 'surv_event', names(var.logrank.p)[is.na(var.logrank.p)])
process.settings <-
  list(
    var        = vars.noprocess,
    method     = rep(NA, length(vars.noprocess)),
    settings   = rep(list(NA), length(vars.noprocess))
  )
# Find continuous variables which will need discretising
continuous.vars <- names(COHORT.bigdata)[sapply(COHORT.bigdata, class) %in% c('integer', 'numeric')]
# Remove those variables already explicitly excluded, mainly for those whose
# logrank score was NA
continuous.vars <- continuous.vars[!(continuous.vars %in% process.settings$var)]
process.settings$var <- c(process.settings$var, continuous.vars)
process.settings$method <-
  c(process.settings$method,
    rep('binByQuantile', length(continuous.vars))
  )
process.settings$settings <-
  c(
    process.settings$settings,
    rep(
      list(
        seq(
          # Quantiles are obviously between 0 and 1
          0, 1,
          # Choose a random number of bins (and for n bins, you need n + 1 breaks)
          length.out = 10
        )
      ),
      length(continuous.vars)
    )
  )

COHORT.prep <-
  prepData(
    # Data for cross-validation excludes test set
    COHORT.bigdata,
    names(COHORT.bigdata),
    process.settings,
    'surv_time', 'surv_event',
    TRUE
  )

# Kludge...remove surv_time.1 and rename surv_event.1
COHORT.prep$surv_time.1 <- NULL
names(COHORT.prep)[names(COHORT.prep) == 'surv_event.1'] <- 'surv_event'

# Further kludge...remove negative survival times
COHORT.prep <- subset(COHORT.prep, surv_time > 0)

#' ## Run random forest calibration
#' 
#' If there's not already a calibration file, we run the rfVarSel methodology:
#'   1. Fit a big forest to the whole dataset to obtain variable importances.
#'   2. Cross-validate as number of most important variables kept is reduced.
#' 
#' (If there is already a calibration file, just load the previous work.)
#' 
#+ rf_var_sel_calibration

# If we've not already done a calibration, then do one
if(!file.exists(calibration.filename)) {
  # Create an empty data frame to aggregate stats per fold
  cv.performance <- data.frame()
  
  # Cross-validate over number of variables to try
  cv.vars <-
    getVarNums(
      length(var.logrank.p),
      # no point going lower than the point at which all the p-values are 0,
      # because the order is alphabetical and therefore meaningless below this!
      min =  sum(var.logrank.p == 0, na.rm = TRUE)
    )
  
  COHORT.cv <- COHORT.prep[-test.set, ]
  
  # Run crossvalidations. No need to parallelise because rfsrc is parallelised
  for(i in 1:length(cv.vars)) {
    # Get the subset of most important variables to use
    surv.predict.partial <- names(var.logrank.p)[1:cv.vars[i]]
    
    # Get folds for cross-validation
    cv.folds <- cvFolds(nrow(COHORT.cv), cv.n.folds)
    
    cv.fold.performance <- data.frame()
    
    for(j in 1:cv.n.folds) {
      time.start <- handyTimer()
      # Fit model to the training set
      surv.model.fit <-
        survivalFit(
          surv.predict.partial,
          COHORT.cv[-cv.folds[[j]],],
          model.type = 'survreg',
          n.threads = n.threads
        )
      time.learn <- handyTimer(time.start)
      
      time.start <- handyTimer()
      # Get C-index on validation set
      c.index.val <-
        cIndex(
          surv.model.fit, COHORT.cv[cv.folds[[j]],]
        )
      time.c.index <- handyTimer(time.start)
      
      time.start <- handyTimer()
      # Get calibration score validation set
      calibration.score <-
        calibrationScore(
          calibrationTable(
            surv.model.fit, COHORT.cv[cv.folds[[j]],]
          )
        )
      time.calibration <- handyTimer(time.start)
      
      # Append the stats we've obtained from this fold
      cv.fold.performance <-
        rbind(
          cv.fold.performance,
          data.frame(
            calibration = i,
            cv.fold = j,
            n.vars = cv.vars[i],
            c.index.val,
            calibration.score,
            time.learn,
            time.c.index,
            time.calibration
          )
        )
      
    } # End cross-validation loop (j)
    
    
    # rbind the performance by fold
    cv.performance <-
      rbind(
        cv.performance,
        cv.fold.performance
      )
    
    # Save output at the end of each loop
    write.csv(cv.performance, calibration.filename)
    
  } # End calibration loop (i)
} else {
  cv.performance <- read.csv(calibration.filename)
}

#' ## Find the best model from the calibrations
#' 
#' ### Plot model performance
#' 
#+ model_performance

# Find the best calibration...
# First, average performance across cross-validation folds
cv.performance.average <-
  aggregate(
    c.index.val ~ n.vars,
    data = cv.performance,
    mean
  )

cv.calibration.average <-
  aggregate(
    area ~ n.vars,
    data = cv.performance,
    mean
  )

ggplot(cv.performance.average, aes(x = n.vars, y = c.index.val)) +
  geom_line() +
  geom_point(data = cv.performance) +
  ggtitle(label = 'C-index by n.vars')

ggplot(cv.calibration.average, aes(x = n.vars, y = area)) +
  geom_line() +
  geom_point(data = cv.performance) +
  ggtitle(label = 'Calibration performance by n.vars')

# Find the highest value
n.vars <-
  cv.performance.average$n.vars[
    which.max(cv.performance.average$c.index.val)
    ]

# Fit a full model with the variables provided
surv.predict.partial <- names(var.logrank.p)[1:n.vars]

#' ## Best model
#' 
#' The best model contained `r n.vars` variables. Let's see what those were...
#' 
#+ variables_used

vars.df <-
  data.frame(
    vars = surv.predict.partial
  )

vars.df$descriptions <- lookUpDescriptions(surv.predict.partial)

vars.df$missingness <- missingness[surv.predict.partial]

#+ variables_table, results='asis'

print(
  xtable(vars.df),
  type = 'html',
  include.rownames = FALSE
)

#' ## Perform the final fit
#' 
#' Having found the best number of variables by cross-validation, let's perform
#' the final fit with the full training set and `r n.trees.final` trees.
#' 
#+ final_fit

time.start <- handyTimer()
surv.model.fit.final <-
  survivalFit(
    surv.predict.partial,
    COHORT.prep[-test.set,],
    model.type = 'survreg'
  )
time.fit.final <- handyTimer(time.start)

saveRDS(surv.model.fit.final, paste0(output.filename.base, '-finalmodel.rds'))

#' Final model of `r n.trees.final` trees fitted in `r round(time.fit.final)`
#' seconds! 
#' 
#' Also bootstrap this final fitting stage. A fully proper bootstrap would
#' iterate over the whole model-building process including variable selection,
#' but that would be prohibitive in terms of computational time.
#' 
#+ bootstrap_final

surv.model.fit.boot <-
  survivalBootstrap(
    surv.predict.partial,
    COHORT.bigdata[-test.set,], # Training set
    COHORT.bigdata[test.set,],  # Test set
    model.type = 'survreg',
    bootstraps = bootstraps
  )

# Get coefficients and variable importances from bootstrap fits
surv.model.fit.coeffs <- bootStats(surv.model.fit.boot, uncertainty = '95ci')

#' ## Performance
#' 
#' ### C-index
#' 
#' C-indices are **`r round(surv.model.fit.coeffs['c.train', 'val'], 3)`
#' (`r round(surv.model.fit.coeffs['c.train', 'lower'], 3)` - 
#' `r round(surv.model.fit.coeffs['c.train', 'upper'], 3)`)**
#' on the training set and
#' **`r round(surv.model.fit.coeffs['c.test', 'val'], 3)`
#' (`r round(surv.model.fit.coeffs['c.test', 'lower'], 3)` - 
#' `r round(surv.model.fit.coeffs['c.test', 'upper'], 3)`)** on the test set.
#' 
#'
#' ### Calibration
#' 
#' The bootstrapped calibration score is
#' **`r round(surv.model.fit.coeffs['calibration.score', 'val'], 3)`
#' (`r round(surv.model.fit.coeffs['calibration.score', 'lower'], 3)` - 
#' `r round(surv.model.fit.coeffs['calibration.score', 'upper'], 3)`)**.
#' 
#' Let's draw a representative curve from the unbootstrapped fit... (It would be
#' better to draw all the curves from the bootstrap fit to get an idea of
#' variability, but I've not implemented this yet.)
#' 
#+ calibration_plot

calibration.table <-
  calibrationTable(
    # Standard calibration options
    surv.model.fit.final, COHORT.bigdata[test.set,],
    # Always need to specify NA imputation for rfsrc
    na.action = 'na.impute'
  )

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' The area between the calibration curve and the diagonal is 
#' **`r round(calibration.score[['area']], 3)`** +/-
#' **`r round(calibration.score[['se']], 3)`**.
#'
#+ save_results

# Save performance results
varsToTable(
  data.frame(
    model = 'cox-logrank',
    imputation = FALSE,
    discretised = TRUE,
    c.index = surv.model.fit.coeffs['c.test', 'val'],
    c.index.lower = surv.model.fit.coeffs['c.test', 'lower'],
    c.index.upper = surv.model.fit.coeffs['c.test', 'upper'],
    calibration.score = surv.model.fit.coeffs['calibration.score', 'val'],
    calibration.score.lower =
      surv.model.fit.coeffs['calibration.score', 'lower'],
    calibration.score.upper =
      surv.model.fit.coeffs['calibration.score', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)