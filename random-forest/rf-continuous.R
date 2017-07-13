#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Trying out continuous variables in random survival forests
#' 
#' So far, I've been looking at binning to account for missing data in random
#' forests. Let's try using RFSRC's imputation function, together with a flag
#' to indicate missingness.

#+ user_variables, message=FALSE

data.filename <- '../../data/cohort-sanitised.csv'
output.base <- '../../output/rf-continuous-try7'

endpoint <- 'death' # Change to MI to look for MI...anything else uses death

n.trees <- 5000
n.imputations <- 5
n.splits <- 10
n.threads <- 40
# What to do with missing data
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

surv.predict <- c(surv.predict, names(COHORT.use)[grepl('_missing', names(COHORT.use))])

# Fit random forest
surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.use[-test.set,],
    model.type = 'rfsrc',
    n.trees = n.trees,
    nsplit = n.splits,
    nimpute = n.imputations,
    n.threads = n.threads,
    split.rule = 'logrank',
    na.action = 'na.impute'
  )

saveRDS(surv.model.fit, paste0(output.base, '-surv-forest.rds'))

# Get C-indices for training and test sets
c.index.train <- cIndex(surv.model.fit, COHORT.use[-test.set, ],
                        na.action = 'na.impute')
c.index.test <- cIndex(surv.model.fit, COHORT.use[test.set, ],
                       na.action = 'na.impute')

#' # Results
#' 
#' ## Performance
#' 
#' The C-index on the full training set is **`r round(c.index.train, 3)`**.
#' 
#' The C-index on the held-out test set is **`r round(c.index.test, 3)`**.
#' 
#' 
calibration.table <-
  calibrationTable(fit.exp, COHORT.scaled.demissed[test.set, ])

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' The area between the calibration curve and the diagonal is 
#' **`r round(calibration.score[['area']], 3)`** +/-
#' **`r round(calibration.score[['se']], 3)`**.
#' 

# Save the model performance
varsToTable(
  data.frame(
    model = 'rfsrc',
    imputation = FALSE,
    discretised = FALSE,
    c.index = c.index.test,
    c.index.lower = NA, # bootstrapping not yet implemented
    c.index.upper = NA,
    calibration.score = calibration.score[['area']],
    calibration.score.lower = NA, # bootstrapping not yet implemented
    calibration.score.upper = NA
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)


#' ## Model fit
#' 
#+ resulting_fit

print(surv.model.fit)

#' ## Effects of variables
#' 
#' Whilst random forests are often considered something of a black box, it is
#' possible to get some insight into how each variable affects the output simply
#' by varying it and seeing what happens. That's what these functions do: take,
#' say, 1000 random patients, and run each of them through the model many times
#' with different values of the variable of interest, and see how the resulting
#' risk varies. Let's try it...
#' 
#+ rf_variable_effects

for(variable in continuous.vars) {
  risk.by.variable <- generalEffectDf(surv.model.fit, COHORT.prep, variable)
  
  # Get the mean of the normalised risk for every value of the variable
  risk.aggregated <-
    aggregate(
      as.formula(paste0('risk.normalised ~ ', variable)),
      risk.by.variable, mean
    )
  
  # work out the limits on the x-axis by taking the 1st and 99th percentiles
  x.axis.limits <-
    quantile(COHORT.full[, variable], c(0.01, 0.99), na.rm = TRUE)
  
  print(
    ggplot(risk.by.variable, aes_string(x = variable, y = 'risk.normalised')) +
      geom_line(alpha=0.01, aes(group = id)) +
      geom_line(data = risk.aggregated, colour = 'blue') +
      coord_cartesian(xlim = c(x.axis.limits))
  )
}
