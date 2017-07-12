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
n.threads <- 16

source('../lib/shared.R')
require(rms)

COHORT.full <- data.frame(fread(data.filename))

COHORT.use <- subset(COHORT.full, !exclude)

n.data <- nrow(COHORT.use)

# Define indices of test set
test.set <- testSetIndices(COHORT.use, random.seed = 78361)

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

ageOnly <- function(df, indices, df.test) {
  # Create a Kaplan-Meier curve from the bootstrap sample
  km.by.age <-
    survfit(
      Surv(surv_time, surv_event) ~ age,
      data = df[indices, ],
      conf.type = "log-log"
    )

  # Return the calibration score
  c(
    calibration.score =
      calibrationScore(
        calibrationTable(km.by.age, df.test)
      )$area
  )
}

# C-index is by definition non-variable because we do not bootstrap or otherwise
# permute the test set, and there's no 'training' phase
age.c.index <- 
  as.numeric(
    survConcordance(
      Surv(surv_time, surv_event) ~ age,
      COHORT.use[test.set,]
    )$concordance
  )

# Bootstrap to establish variability of calibration
age.only.boot <-
  boot(
    data = COHORT.use[-test.set,],
    statistic = ageOnly,
    R = bootstraps,
    parallel = 'multicore',
    ncpus = n.threads,
    df.test =  COHORT.use[test.set,]
  )

age.only.boot.stats <- bootStats(age.only.boot, uncertainty = '95ci')

#' C-index is
#' **`r round(age.c.index, 3)`** 
#' on the held-out test set (not that it really matters, the model isn't
#' 'trained' as such for the discrimination test...it's just oldest patient dies
#' first).
#' 
#' 
#' ## Calibration
#' 
#' 
#+ calibration

km.by.age <-
  survfit(
    Surv(surv_time, surv_event) ~ age,
    data = COHORT.use[-test.set,],
    conf.type = "log-log"
  )

calibration.table <- calibrationTable(km.by.age, COHORT.use[test.set, ])

print(calibrationScore(calibration.table))

calibrationPlot(calibration.table)

#' Calibration score is
#' **`r round(age.only.boot.stats['calibration.score', 'val'], 3)`
#' (`r round(age.only.boot.stats['calibration.score', 'lower'], 3)` -
#' `r round(age.only.boot.stats['calibration.score', 'upper'], 3)`)** 
#' on the held-out test set.

varsToTable(
  data.frame(
    model = 'age',
    imputation = FALSE,
    discretised = FALSE,
    c.index = age.c.index,
    c.index.lower = NA,
    c.index.upper = NA,
    calibration.score = age.only.boot.stats['calibration.score', 'val'],
    calibration.score.lower = age.only.boot.stats['calibration.score', 'lower'],
    calibration.score.upper = age.only.boot.stats['calibration.score', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)
