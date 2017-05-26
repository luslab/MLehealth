#' # Trying with bigger data
#' This ignores the hand-picked variables used in the original study, and fits
#' only to my big-data-style most common diagnoses/measurements etc.

data.filename.big <- '../../data/cohort-datadriven.csv'
model.type <- 'rfsrc'
n.data <- NA

source('../lib/shared.R')

COHORT.full <- fread(data.filename.big)

if(!is.na(n.data)){
  # Take a subset n.data in size
  COHORT.use <- sample.df(COHORT.full, n.data)
  rm(COHORT.full)
} else {
  # Use all the data
  COHORT.use <- COHORT.full
  rm(COHORT.full)
}

COHORT.prep <-
  prepData(
    as.data.frame(COHORT.use),
    cols.keep, discretise.settings, surv.time, surv.event,
    surv.event.yes, extra.fun = caliberExtraPrep, n.keep = n.data
  )
n.data <- nrow(COHORT.prep)

# Define indices of test set
test.set <- sample(1:n.data, (1/3)*n.data)

# Columns which could contain missing data and we'd like to leave it missing
# because they're continuous
continuous.vars <-
  c(
    # A few of the original cohort variables, to keep things fair
    c('age', 'smokstatus', 'imd_score'),
    # Those starting tests.enttype. are tests with continuous results
    names(COHORT.use)[substring(names(COHORT.use), 1, nchar('tests.enttype.')) == 'tests.enttype.']
  )

# Columns containing time since x diagnosis type variables, where we can
# legitimately replace all missing values with a number less than zero
time.vars <-
  c(
     names(COHORT.use)[substring(names(COHORT.use), 1, nchar('hes.icd.')) == 'hes.icd.'],
     names(COHORT.use)[substring(names(COHORT.use), 1, nchar('hes.opcs.')) == 'hes.opcs.']
  )

# Those variables which we don't want to touch
untransformed.vars <- c('surv_time', 'imd_score', 'exclude')

# Aggregate them together, using unique because we've probably double-specified but just to be safe
cols.keep <- unique(c('gender', cols.keep, continuous.vars, time.vars))
surv.predict <- c('gender', continuous.vars, time.vars)

process.settings <-
  list(
    var        = c(untransformed.vars, continuous.vars, time.vars),
    method     =
      # Untransformed vars and vars which are continuous and thus legitimately missing, leave alone
      c(
      	rep(NA, length(untransformed.vars) + length(continuous.vars)),
        # time-like vars, we want to replace missing values with -1
        rep('NA2val', length(time.vars))
      ),
    settings   =
    c(
      rep(NA, length(untransformed.vars) + length(continuous.vars)),
      rep(-1, length(time.vars))
     )
  )

COHORT.prep <-
  prepData(
    # Data for cross-validation excludes test set
    as.data.frame(COHORT.use),
    cols.keep,
    process.settings,
    surv.time, surv.event,
    surv.event.yes,
    extra.fun = caliberExclude
  )
n.data <- nrow(COHORT.prep)

# I have no idea why gender is being preserved as a character vector which then
# causes rfsrc to fail, so let's kludge it for now to get a C-index estimate...
COHORT.prep$gender <- factor(COHORT.prep$gender)
COHORT.prep$smokstatus <- factor(COHORT.prep$smokstatus)
if(!file.exists('../../output/rf-bigdata-only-try1-fit.rds')) {
  time.start <- handyTimer()
  surv.model.fit <-
    survivalFit(
      surv.predict,
      COHORT.prep[-test.set,],
      model.type = model.type,
      n.trees = 2000,
      n.threads = 16,
      split.rule = 'logrank',
      na.action = 'na.impute',
      nimpute = 3,
      nsplit = 8
    )
  time.fit <- handyTimer(time.start)
  
  saveRDS(surv.model.fit, '../../output/rf-bigdata-only-try1-fit.rds')
} else {
  surv.model.fit <- readRDS('../../output/rf-bigdata-only-try1-fit.rds')
  time.fit <- NA
}

time.start <- handyTimer()
c.index <-
  cIndex(surv.model.fit, COHORT.prep[test.set,], na.action = 'na.impute')
time.c.index <- handyTimer(time.start)

#' ## Results
#'
#' Here's the fit:

print(surv.model.fit)

#' C-index is `r round(c.index, 4)`.
#' 
#' It took `r round(time.fit/60, 1)` minutes to fit the model and
#' `r round(time.c.index/60, 1)` minutes to evaluate the C-index on the test set.
#' 
#' ## Variable importance
#' 
#' Can this approach recover the most important variables? Let's see!

variable.importance <-
  generalVarImp(
    surv.model.fit,
    COHORT.prep[test.set,],
    risk.time = 5,
    tod.round = 0.1,
    na.action = 'na.impute'
  )

write.csv(variable.importance, '../../output/rf-bigdata-try1-varimp.csv')

require(xtable)

#+ varimp_table, results='asis'

print(
  xtable(
    variable.importance[
      order(variable.importance$var.imp, decreasing = TRUE)[1:20],
      ],
    digits = c(0,3)
  ),
  type = 'html',
  include.rownames = FALSE
)

#' ## Calibration
#' 
#+ calibration

calibration.table <- calibrationTable(surv.model.fit, COHORT.prep[test.set, ])

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' Calibration score is is `r round(calibration.score$area, 4)` +/-
#' `r round(calibration.score$se, 4)`.