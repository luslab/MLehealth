#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Random survival forests on pre-imputed data
#' 
#' The imputation algorithm used for the Cox modelling 'cheats' in a couple of
#' different ways. Firstly, the imputation model is fitted on the whole dataset,
#' rather than training a model on a training set and then using it on the test
#' set. Secondly, causality is violated in that future measurements and even the
#' outcome variable (ie whether a patient died) is included in the model. The
#' rationale for this is that you want to have the best and most complete
#' dataset possible, and if doctors are going to go on and use this in clinical
#' practice they will simply take the additional required measurements when
#' calculating a patient's risk score, rather than trying to do the modelling on
#' incomplete data.
#' 
#' However, this could allow the imputation to 'pass back' useful data to the
#' Cox model, and thus artificially inflate its performance statistics.
#' 
#' It is non-trivial to work around this, not least because the imputation
#' package we used does not expose the model to allow training on a subset of
#' the data. A quick and dirty method to check whether this may be an issue,
#' therefore, is to train a random forest model on the imputed dataset, and see
#' if it can outperform the Cox model. So, let's try it...

#+ user_variables, message=FALSE

imputed.data.filename <- '../../data/COHORT_complete.rds'
n.data <- NA # This is of full dataset...further rows may be excluded in prep
endpoint <- 'death.imputed'

n.trees <- 500
n.split <- 10
n.threads <- 40

bootstraps <- 50

output.filename.base <- '../../output/rf-imputed-try1'
boot.filename <- paste0(output.filename.base, '-boot.rds')

# What to do with missing data
continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )
untransformed.vars <- c('anonpatid', 'surv_time', 'imd_score', 'exclude')

source('../lib/shared.R')
require(ggrepel)

#' ## Read imputed data

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
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
    caliberExtraPrep(
      prepSurvCol(
        imputed.data[[i]], 'time_death', 'endpoint_death', 1
      )
    )
}

# Define test set
test.set <- testSetIndices(imputed.data[[1]], random.seed = 78361)


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
    nsplit = n.split
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
    model = 'rf',
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
#' Not too bad!
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