#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- TRUE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Trying out continuous variables in random survival forests
#' 
#' So far, I've been looking at binning to account for missing data in random
#' forests. Let's try continuing to treat them continuously, and shunting
#' missing values to one or other end of the data range.

#+ user_variables, message=FALSE

data.filename <- '../../data/cohort-sanitised.csv'

n.trees <- 500

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )
untransformed.vars <- c('anonpatid', 'time_death', 'imd_score', 'exclude')

source('shared.R')
require(ggrepel)

#' ## missingToBig
#' 
#' The continuous variables with missing values will have them shunted to the
#' large end of the range by `missingToBig`. This turns NAs into a value larger
#' than any present in the data, and rounded up nicely so as to be
#' human-readable. For example, if all the data are negative, it uses 0; if the
#' data include 0, it goes with 100; and if the data are positive, it uses 1
#' followed by enough zeros to be at least ten times larger than the largest
#' data value.

print(missingToBig(c(-3.1, NA)))
print(missingToBig(c(0, NA)))
print(missingToBig(c(4.7, NA)))

#' ## Model time!

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- data.frame(fread(data.filename))

# Define process settings; nothing for those to not transform, and missingToBig
# for the continuous ones...
process.settings <-
  list(
    var        = c(untransformed.vars, continuous.vars),
    method     =
      c(
        rep(NA, length(untransformed.vars)),
        rep('missingToBig', length(continuous.vars))
      ),
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
test.set <- sample(1:n.data, (1/3)*n.data)

# Fit random forest
surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.prep[-test.set,],
    model.type = 'ranger',
    n.trees = n.trees,
    split.rule = 'logrank',
    n.threads = 8,
    respect.unordered.factors = 'partition'
  )

# Get C-indices for training and test sets
c.index.train <-
  cIndex(surv.model.fit, COHORT.prep[-test.set, ], model.type = 'ranger')
c.index.test <- 
  cIndex(surv.model.fit, COHORT.prep[test.set, ], model.type = 'ranger')

#' # Results
#' 
#' ## Performance
#' 
#' The C-index on the full training set is **`r round(c.index.train, 3)`**.
#' 
#' The C-index on the held-out test set is **`r round(c.index.test, 3)`**.
#' 
#' ## Model fit
#' 
#+ resulting_fit

print(surv.model.fit)