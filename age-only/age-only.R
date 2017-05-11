#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Performance of an age-only model
#' 
#+ setup, message=FALSE

bootstraps <- 200

data.filename <- '../../data/cohort-sanitised.csv'
n.threads <- 8

source('../lib/shared.R')
require(rms)

COHORT.full <- data.frame(fread(data.filename))

COHORT.use <- subset(COHORT.full, !exclude)

n.data <- nrow(COHORT.use)

# Define indices of test set
test.set <- sample(1:n.data, (1/3)*n.data)

COHORT.use <- prepSurvCol(COHORT.use, surv.time, surv.event,surv.event.yes)

#' ## Concordance index
#' 
#' The c-index can be calculated taking age itself as a risk score. Since it's
#' purely rank-based, it doesn't matter that age is nonlinearly related to true
#' risk.
#' 
#' We bootstrap this based purely on the test set to make the bootstrap
#' variability commensurate with other models tested on the test set. You can
#' imagine that this model was 'trained' on the training set, even though it's
#' so simple that we actually did no such thing...
#' 
#+ c_index

# Define a trivial function used for bootstrapping which simply returns the
# c-index of a 'model' based purely on age.

ageOnly <- function(df, indices) {
  c(
    c.index = 
      as.numeric(
        survConcordance(
          Surv(surv_time, surv_event) ~ age,
          df[indices,]
        )$concordance
      )
  )
}

age.only.c.index <-
  boot(
    data = COHORT.use[test.set,],
    statistic = ageOnly,
    R = bootstraps,
    parallel = 'multicore',
    ncpus = n.threads
  )

age.only.c.index.ci <- bootStats(age.only.c.index, uncertainty = '95ci')

#' The c-index is `r age.only.c.index$val` +/- `r age.only.c.index$err`.
#' 

varsToTable(
  data.frame(
    model = 'age',
    imputation = FALSE,
    discretised = FALSE,
    c.index = age.only.c.index.ci['c.index', 'val'],
    c.index.lower = age.only.c.index.ci['c.index', 'lower'],
    c.index.upper = age.only.c.index.ci['c.index', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)
