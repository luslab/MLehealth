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
output.filename.base <- '../../output/rfsrc-classification-try1'

risk.time <- 5

n.data <- NA
split.rule <- 'logrank'
n.trees <- 2000
n.threads <- 19

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )

untransformed.vars <- c('anonpatid', 'surv_time', 'imd_score', 'exclude')

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

# Process settings: don't touch anything!!
process.settings <-
  list(
    var       = c(untransformed.vars, continuous.vars),
    method    = rep(NA, length(untransformed.vars) + length(continuous.vars)),
    settings  = rep(NA, length(untransformed.vars) + length(continuous.vars))
  )

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

# Create column for whether or not the patient had an event before risk.time
COHORT.prep$event <- NA
# Event before risk.time
COHORT.prep$event[
  COHORT.prep$surv_event & COHORT.prep$surv_time <= risk.time
  ] <- TRUE
# Event after, whether censorship or not, means no event by risk.time
COHORT.prep$event[COHORT.prep$surv_time > risk.time] <- FALSE
# Otherwise, censored before risk.time, let's remove the row
COHORT.prep <- COHORT.prep[!is.na(COHORT.prep$event), ]

surv.model.fit <-
  rfsrc(
    as.formula(
      paste('event ~', paste(surv.predict, collapse = '+'))
    ),
    COHORT.prep[-test.set,], # Training set
    ntree = n.trees,
    splitrule = 'gini',
    n.threads = n.threads,
    na.action = 'na.impute',
    nimpute = 3
  )