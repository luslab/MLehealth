#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Discrete Cox model with imputed data
#' 
#' The discrete Cox model performs very similarly to the normal Cox model, even
#' without performing imputation first. Let's try it with imputation, and see
#' if the performance is boosted.
#' 
#' We're going to use the same parameters for discretising as were found for the
#' data with missing values. On the one hand, this is slightly unfair because
#' these values might be suboptimal now that the data are imputed but, on the
#' other, it does allow for direct comparison. If we cross-validated again, we
#' would be very likely to find different parameters (there are far more than
#' can plausibly be tried), which may lead performance to be better or worse
#' entirely at random.

# Calibration from cross-validation performed on data before imputation
calibration.filename <- '../../output/survreg-crossvalidation-try4.csv'
# The first part of the filename for any output
output.filename.base <- '../../output/survreg-discrete-imputed-try1'
imputed.data.filename <- '../../data/COHORT_complete.rds'
boot.filename <- paste0(output.filename.base, '-boot.rds')

n.threads <- 20
bootstraps <- 50

model.type <- 'survreg'

#' ## Discretise data
#' 
#' First, let's load the results from the calibrations and find the parameters
#' for discretisation.

source('../lib/shared.R')

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )

# read file containing calibrations
cv.performance <- read.csv(calibration.filename)

# Find the best calibration...
# First, average performance across cross-validation folds
cv.performance.average <-
  aggregate(
    c.index.val ~ calibration,
    data = cv.performance,
    mean
  )
# Find the highest value
best.calibration <-
  cv.performance.average$calibration[
    which.max(cv.performance.average$c.index.val)
    ]
# And finally, find the first row of that calibration to get the n.bins values
best.calibration.row1 <-
  min(which(cv.performance$calibration == best.calibration))

# Get its parameters
n.bins <-
  t(
    cv.performance[best.calibration.row1, continuous.vars]
  )


# Reset process settings with the base setings
process.settings <-
  list(
    var        = c('anonpatid', 'time_death', 'imd_score', 'exclude'),
    method     = c(NA, NA, NA, NA),
    settings   = list(NA, NA, NA, NA)
  )
for(j in 1:length(continuous.vars)) {
  process.settings$var <- c(process.settings$var, continuous.vars[j])
  process.settings$method <- c(process.settings$method, 'binByQuantile')
  process.settings$settings <-
    c(
      process.settings$settings,
      list(
        seq(
          # Quantiles are obviously between 0 and 1
          0, 1,
          # Choose a random number of bins (and for n bins, you need n + 1 breaks)
          length.out = n.bins[j]
        )
      )
    )
}

#' ## Load imputed data
#' 
#' Load the data and prepare it with the settings above

# Load the data from its RDS file
imputed.data <- readRDS(imputed.data.filename)

# Remove rows with death time of 0 to avoid fitting errors, and get the survival
# columns ready
for(i in 1:length(imputed.data)) {
  imputed.data[[i]] <- imputed.data[[i]][imputed.data[[i]][, surv.time] > 0, ]
  # Put in a fake exclude column for the next function (exclusions are already
  # excluded in the imputed dataset)
  imputed.data[[i]]$exclude <- FALSE
  imputed.data[[i]]$imd_score <- as.numeric(imputed.data[[i]]$imd_score)
  imputed.data[[i]] <-
    prepData(
      imputed.data[[i]],
      cols.keep,
      process.settings,
      'time_death', 'endpoint_death', 1,
      extra.fun = caliberExtraPrep
    )
}

# Define test set
test.set <- testSetIndices(imputed.data[[1]], random.seed = 78361)

# Do a quick and dirty fit on a single imputed dataset, to draw calibration
# curve from
time.start <- handyTimer()
surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.optimised[-test.set,], # Training set
    model.type = model.type,
    n.threads = n.threads
  )
time.fit <- handyTimer(time.start)

# Perform bootstrap fitting for every multiply imputed dataset
surv.fit.boot <- list()
time.start <- handyTimer()
for(i in 1:length(imputed.data)) {
  surv.fit.boot[[i]] <-
    survivalBootstrap(
      surv.predict,
      imputed.data[[i]][-test.set, ],
      imputed.data[[i]][test.set, ],
      model.type = model.type,
      bootstraps = bootstraps,
      n.threads = n.threads
    )
}
time.boot <- handyTimer(time.start)

# Save the fits, because it might've taken a while!
saveRDS(surv.fit.boot, boot.filename)

#' Model fitted in `r round(time.fit)` seconds, and `r bootstraps` fits
#' performed on `r length(imputed.data)` imputed datasets in
#' `r round(time.boot)` seconds.

# Unpackage the uncertainties from the bootstrapped data
surv.fit.boot.ests <-  bootMIStats(surv.fit.boot)

# Save bootstrapped performance values
varsToTable(
  data.frame(
    model = 'cox',
    imputation = TRUE,
    discretised = TRUE,
    c.index = surv.fit.boot.ests['c.test', 'val'],
    c.index.lower = surv.fit.boot.ests['c.test', 'lower'],
    c.index.upper = surv.fit.boot.ests['c.test', 'upper'],
    calibration.score = surv.fit.boot.ests['calibration.score', 'val'],
    calibration.score.lower = surv.fit.boot.ests['calibration.score', 'lower'],
    calibration.score.upper = surv.fit.boot.ests['calibration.score', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)


#' # Results
#' 
#' ## Performance
#' 
#' Having fitted the Cox model, how did we do? The c-indices were calculated as
#' part of the bootstrapping, so we just need to take a look at those...
#' 
#' C-indices are **`r round(surv.fit.boot.ests['c.train', 'val'], 3)`
#' (`r round(surv.fit.boot.ests['c.train', 'lower'], 3)` - 
#' `r round(surv.fit.boot.ests['c.train', 'upper'], 3)`)** on the training set and
#' **`r round(surv.fit.boot.ests['c.test', 'val'], 3)`
#' (`r round(surv.fit.boot.ests['c.test', 'lower'], 3)` - 
#' `r round(surv.fit.boot.ests['c.test', 'upper'], 3)`)** on the test set.
#' 
#' 
#' ### Calibration
#' 
#' The bootstrapped calibration score is
#' **`r round(surv.fit.boot.ests['calibration.score', 'val'], 3)`
#' (`r round(surv.fit.boot.ests['calibration.score', 'lower'], 3)` - 
#' `r round(surv.fit.boot.ests['calibration.score', 'upper'], 3)`)**.
#' 
#' Let's draw a representative curve from the unbootstrapped fit... (It would be
#' better to draw all the curves from the bootstrap fit to get an idea of
#' variability, but I've not implemented this yet.)
#' 
#+ calibration_plot

calibration.table <-
  calibrationTable(surv.model.fit, imputed.data[[1]][test.set, ])

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' The area between the calibration curve and the diagonal is 
#' **`r round(calibration.score[['area']], 3)`** +/-
#' **`r round(calibration.score[['se']], 3)`**.