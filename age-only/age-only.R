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

#' The c-index is `r age.only.c.index$val`
#' (`r age.only.c.index$lower`-`r age.only.c.index$upper`).
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

km.df <- data.frame(
  age = rep(
    # Chop off first 4 characters ('age=') and turn age into a number
    as.numeric(substring(names(km.by.age$strata), 5)),
    # Repeat each number as many times as there are patients that age
    times = km.by.age$strata
  ),
  time = km.by.age$time,
  surv = km.by.age$surv
)

risk.time <- 5
surv.by.age <-
  data.frame(
    age = unique(km.df$age),
    surv.at.t = NA
  )

for(age in unique(km.df$age)) {
  # If anyone in that age bracket lived long enough for us to make a prediction...
  if(max(km.df$time[km.df$age == age]) > risk.time) {
    # Find the first event after that point, which gives us the survival
    surv.by.age$surv.at.t[surv.by.age$age == age] <-
      km.df$surv[
        # The datapoint needs to be for the correct age of patient
        km.df$age == age &
        # And pick the time which is the smallest value greater than the time
        # in which we're interested.
        km.df$time ==
          minGt(km.df$time[km.df$age == age], risk.time)
      ]
  }
}

# This is basically a kludged version of the calibrationTable() function for
# this setting, which is so simple as not to merit functionalising
calibration.table.precursor <-
  merge(
    COHORT.use[, c('age', 'surv_time', 'surv_event'), drop = FALSE],
    surv.by.age[, c('age', 'surv.at.t')],
    sort = FALSE, all.x = TRUE
  )
calibration.table <-
  data.frame(
    risk = 1 - calibration.table.precursor$surv.at.t,
    event = NA
  )
# Event before risk.time
calibration.table$event[calibration.table.precursor$surv_event & calibration.table.precursor$surv_time <= risk.time] <- TRUE
# Event after risk.time, whether censored or not
calibration.table$event[calibration.table.precursor$surv_time > risk.time] <- FALSE
# Otherwise, censored before risk.time, leave as NA

calibrationScore(sample.df(calibration.table, 30000))

calibrationPlot(sample.df(calibration.table, 30000))

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
