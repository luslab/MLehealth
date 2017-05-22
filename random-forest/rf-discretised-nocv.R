#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Cross-validating discretisation of input variables in a survival model
#' 
#' In difference to previous attempts at cross-validation, this uses between 10
#' and 20 bins, not between 2 and 20, in an attempt to avoid throwing away data.

# The first part of the filename for any output
output.filename.base <- '../../output/rfsrc-discrete-nocv-try1'

# What kind of model to fit to...currently 'cph' (Cox model), 'ranger' or
# 'rfsrc' (two implementations of random survival forests)
model.type <- 'rfsrc'

n.data <- NA
split.rule <- 'logrank'
n.trees <- 2000
n.threads <- 19

# If surv.vars is defined as a character vector here, the model only uses those
# variables specified, eg c('age') would build a model purely based on age. If
# not specified (ie commented out), it will use the defaults.
# surv.predict <- c('age')

#' ## Fit the model
#' 
#' Now, let's fit the model, but without cross-validating the number of factor
#' levels! The issue is that, if we're allowing factor levels to be grouped into
#' two branches arbitrarily, there are 2^n - 2 combinations, which rapidly
#' becomes a huge number. Thus, cross-validating, especially with large numbers
#' of factor levels, is very impractical.
#' 
#' We'll also leave age as a pure number: we know that it's both a very
#' significant variable, and also that it makes sense to treat it as though it's
#' ordered, because risk should increase monotonically with it.
#' 
#+ rf_discretised, cache=cacheoption

source('../lib/shared.R')

continuous.vars <-
  c(
    'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )

process.settings <-
  list(
    var        = c(
      'anonpatid', 'time_death', 'imd_score', 'exclude', 'age',
      continuous.vars
    ),
    method     = c(
      NA, NA, NA, NA, NA,
      rep('binByQuantile', length(continuous.vars))
      ),
    settings   = 
      c(
        # Variables to leave alone have no process settings
        list(NA, NA, NA, NA, NA),
        # The other process settings are passed as a list...
        rep(
          list(
              seq(
                0, 1, # Quantiles are obviously between 0 and 1
                length.out = 9 # Let's have 8 bins, so 9 breaks
              )
            ),
          # Need an item in this list per continuous variable to be binned
          length(continuous.vars)
        )
      )
  )

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- data.frame(fread(data.filename))

# If n.data was specified...
if(!is.na(n.data)){
  # Take a subset n.data in size
  COHORT.use <- sample.df(COHORT.full, n.data)
  rm(COHORT.full)
} else {
  # Use all the data
  COHORT.use <- COHORT.full
  rm(COHORT.full)
}

# We now need a quick null preparation of the data to get its length (some rows
# may be excluded during preparation)
COHORT.prep <-
  prepData(
    COHORT.use,
    cols.keep,
    process.settings,
    surv.time, surv.event,
    surv.event.yes,
    extra.fun = caliberExtraPrep
  )
n.data <- nrow(COHORT.prep)

# Define indices of test set
test.set <- sample(1:n.data, (1/3)*n.data)

surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.prep[-test.set,], # Training set
    model.type = model.type,
    n.trees = n.trees,
    split.rule = split.rule,
    n.threads = n.threads
  )

surv.model.fit.boot <-
  survivalBootstrap(
    surv.predict,
    COHORT.optimised[-test.set,], # Training set
    COHORT.optimised[test.set,],  # Test set
    model.type = model.type,
    n.trees = n.trees,
    split.rule = split.rule,
    n.threads = n.threads,
    bootstraps = bootstraps
  )


# Save the fit object
saveRDS(
  surv.model.fit.boot,
  paste0(output.filename.base, '-surv-model-bootstraps.rds')
)

# Get C-indices for training and test sets
surv.model.fit.coeffs <- bootStats(surv.model.fit.boot, uncertainty = '95ci')

# Save them to the all-models comparison table
varsToTable(
  data.frame(
    model = 'rfsrc',
    imputation = FALSE,
    discretised = FALSE,
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

#' # Results
#' 
#' 
#' #' ## Performance
#' 
#' C-indices are **`r round(surv.model.fit.coeffs['c.train', 'val'], 3)` +/-
#' `r round(surv.model.fit.coeffs['c.train', 'err'], 3)`** on the training set and
#' **`r round(surv.model.fit.coeffs['c.test', 'val'], 3)` +/-
#' `r round(surv.model.fit.coeffs['c.test', 'err'], 3)`** on the held-out test set.
#' 
#' ## Model fit
#' 
#+ resulting_fit

print(surv.model.fit)

#' ## Variable importance

# First, load data from Cox modelling for comparison
cox.var.imp <- read.csv(comparison.filename)

# Then, get the variable importance from the model just fitted
var.imp <-
  data.frame(
    var.imp = importance(surv.model.fit)/max(importance(surv.model.fit))
  )
var.imp$quantity <- rownames(var.imp)

var.imp <- merge(var.imp, cox.var.imp)

# Save the results as a CSV
write.csv(var.imp, paste0(output.filename, '-var-imp.csv'))

#' ## Variable importance vs Cox model replication variable importance

print(
  ggplot(var.imp, aes(x = our_range, y = var.imp)) +
    geom_point() +
    geom_text_repel(aes(label = quantity)) +
    # Log both...old coefficients for linearity, importance to shrink range!
    scale_x_log10() +
    scale_y_log10()
)

print(cor(var.imp[, c('var.imp', 'our_range')], method = 'spearman'))