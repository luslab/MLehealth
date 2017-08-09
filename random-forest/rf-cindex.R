#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Random survival forests split on C-index
#' 
#' Most of our random forests are split by maximising the difference between
#' branches as measured by the logrank test. In theory, since one of our two
#' performance measures is C-index, splitting on C-index instead should allow
#' a better performance. Let's try it...

#+ user_variables, message=FALSE

data.filename <- '../../data/cohort-sanitised.csv'
endpoint <- 'death'

n.trees <- 500
n.split <- 10
n.impute <- 3
n.threads <- 20

bootstraps <- 50

output.filename.base <- '../../output/rf-cindex-try1'
boot.filename <- paste0(output.filename.base, '-boot.rds')

# What to do with missing data
continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )
untransformed.vars <- c('anonpatid', 'surv_time', 'imd_score', 'exclude')
continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )
untransformed.vars <- c('anonpatid', 'surv_time', 'imd_score', 'exclude')

source('../lib/shared.R')
require(ggrepel)

missingReplace <- NA # ie, don't change the values!

#' ## Fit random forest model

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- data.frame(fread(data.filename))

# Define process settings; at first, we're going to do nothing to anything
# continuous, because we may want to make _missing columns which this function
# can't do...
process.settings <-
  list(
    var        = c(untransformed.vars, continuous.vars),
    method     = rep(NA, length(untransformed.vars) + length(continuous.vars)),
    settings   = rep(NA, length(untransformed.vars) + length(continuous.vars))
  )

COHORT.prep <-
  prepData(
    # Data for cross-validation excludes test set
    COHORT.full,
    cols.keep,
    process.settings,
    surv.time, surv.event,
    surv.event.yes,
    extra.fun = caliberExtraPrep
  )
n.data <- nrow(COHORT.prep)

# Define indices of test set
test.set <- testSetIndices(COHORT.prep, random.seed = 78361)

# Now, process the data such that we can fit a model to it...
COHORT.use <- prepCoxMissing(COHORT.prep, missingReplace = missingReplace)

#' ## Fit random forest model

# Do a quick and dirty fit on a single imputed dataset, to draw calibration
# curve from
time.start <- handyTimer()
surv.model.fit <-
  survivalFit(
    surv.predict,
    imputed.data[[1]][-test.set,],
    model.type = 'rfsrc',
    n.trees = n.trees,
    split.rule = split.rule,
    n.threads = n.threads,
    nsplit = n.split,
    na.action = 'na.impute',
    nimpute = n.impute
  )
time.fit <- handyTimer(time.start)

fit.rf.boot <- list()

# Perform bootstrap fitting for every multiply imputed dataset
time.start <- handyTimer()
for(i in 1:length(imputed.data)) {
  fit.rf.boot[[i]] <-
    survivalBootstrap(
      surv.predict,
      imputed.data[[i]][-test.set, ],
      imputed.data[[i]][test.set, ],
      model.type = 'rfsrc',
      bootstraps = bootstraps,
      n.threads = n.threads
    )
}
time.boot <- handyTimer(time.start)

# Save the fits, because it might've taken a while!
saveRDS(fit.rf.boot, boot.filename)

#' Model fitted in `r round(time.fit)` seconds, and `r bootstraps` fits
#' performed on `r length(imputed.data)` imputed datasets in
#' `r round(time.boot)` seconds.

# Unpackage the uncertainties from the bootstrapped data
fit.rf.boot.ests <-  bootMIStats(fit.rf.boot)

# Save bootstrapped performance values
varsToTable(
  data.frame(
    model = 'rfsrc',
    imputation = TRUE,
    discretised = FALSE,
    c.index = fit.rf.boot.ests['c.test', 'val'],
    c.index.lower = fit.rf.boot.ests['c.test', 'lower'],
    c.index.upper = fit.rf.boot.ests['c.test', 'upper'],
    calibration.score = fit.rf.boot.ests['calibration.score', 'val'],
    calibration.score.lower = fit.rf.boot.ests['calibration.score', 'lower'],
    calibration.score.upper = fit.rf.boot.ests['calibration.score', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)

#' ## Performance
#' 
#' Having fitted the Cox model, how did we do? The c-indices were calculated as
#' part of the bootstrapping, so we just need to take a look at those...
#' 
#' C-indices are **`r round(fit.rf.boot.ests['c.train', 'val'], 3)`
#' (`r round(fit.rf.boot.ests['c.train', 'lower'], 3)` - 
#' `r round(fit.rf.boot.ests['c.train', 'upper'], 3)`)** on the training set and
#' **`r round(fit.rf.boot.ests['c.test', 'val'], 3)`
#' (`r round(fit.rf.boot.ests['c.test', 'lower'], 3)` - 
#' `r round(fit.rf.boot.ests['c.test', 'upper'], 3)`)** on the test set.
#' 
#' 
#' ### Calibration
#' 
#' The bootstrapped calibration score is
#' **`r round(fit.rf.boot.ests['calibration.score', 'val'], 3)`
#' (`r round(fit.rf.boot.ests['calibration.score', 'lower'], 3)` - 
#' `r round(fit.rf.boot.ests['calibration.score', 'upper'], 3)`)**.
#' 
#' Let's draw a representative curve from the unbootstrapped fit... (It would be
#' better to draw all the curves from the bootstrap fit to get an idea of
#' variability, but I've not implemented this yet.)
#' 
#+ calibration_plot

calibration.table <-
  calibrationTable(surv.model.fit, imputed.data[[i]][test.set, ])

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' The area between the calibration curve and the diagonal is 
#' **`r round(calibration.score[['area']], 3)`** +/-
#' **`r round(calibration.score[['se']], 3)`**.