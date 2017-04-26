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
#' forests. Let's try continuing to treat them continuously, and shunting
#' missing values to one or other end of the data range.

#+ user_variables, message=FALSE

data.filename <- '../../data/cohort-sanitised.csv'
output.base <- '../output/rf-continuous-try2'

endpoint <- 'death' # Change to MI to look for MI...anything else uses death

n.trees <- 500
n.threads <- 8

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )
untransformed.vars <- c('anonpatid', 'surv_time', 'imd_score', 'exclude')

source('../lib/shared.R')
require(ggrepel)

#' ## Fit random forest model

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
    n.threads = n.threads
  )

saveRDS(surv.model.fit, paste0(output.base, '-surv-forest.rds'))

# Get C-indices for training and test sets
c.index.train <- cIndex(surv.model.fit, COHORT.prep[-test.set, ])
c.index.test <- cIndex(surv.model.fit, COHORT.prep[test.set, ])

# Save the C-index on the test set as the model performance
varsToTable(
  data.frame(
    model = 'random_forest',
    imputation = FALSE,
    discretised = FALSE,
    c.index = c.index.test,
    c.index.lower = NA, # bootstrapping not yet implemented
    c.index.upper = NA
  ),
  performance.file
)

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
