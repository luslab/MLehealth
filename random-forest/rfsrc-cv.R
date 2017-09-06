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

output.filename.base <- '../../output/rfsrc-cv-nsplit-try3'
data.filename <- '../../data/cohort-sanitised.csv'

# If surv.vars is defined as a character vector here, the model only uses those
# variables specified, eg c('age') would build a model purely based on age. If
# not specified (ie commented out), it will use the defaults.
# surv.predict <- c('age')

#' ## Do the cross-validation
#' 
#' The steps for this are common regardless of model type, so run the script to
#' get a cross-validated model to further analyse...
#+ rf_discretised_cv, cache=cacheoption

source('../lib/rfsrc-cv-nsplit-bootstrap.R', chdir = TRUE)

#' # Results
#' 
#' 
#' ## Performance
#' 
#' ### Discrimination
#' 
#' C-indices are **`r round(surv.model.fit.coeffs['c.train', 'val'], 3)` +/-
#' `r round(surv.model.fit.coeffs['c.train', 'err'], 3)`** on the training set and
#' **`r round(surv.model.fit.coeffs['c.test', 'val'], 3)` +/-
#' `r round(surv.model.fit.coeffs['c.test', 'err'], 3)`** on the held-out test set.
#' 
#' ### Calibration
#' 
#' Does the model predict realistic probabilities of an event?
#' 
#+ calibration_plot

calibration.table <-
  calibrationTable(
    # Standard calibration options
    surv.model.fit, COHORT.prep[test.set,],
    # Always need to specify NA imputation for rfsrc
    na.action = 'na.impute'
  )

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' The area between the calibration curve and the diagonal is 
#' **`r round(calibration.score['area'], 3)`** +/-
#' **`r round(calibration.score['se'], 3)`**.
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

#' ## Variable effects
#' 
#+ variable_effects

risk.by.variables <- data.frame()

for(variable in continuous.vars) {
  # Create a partial effect table for this variable
  risk.by.variable <-
    partialEffectTable(
      surv.model.fit, COHORT.prep[-test.set,], variable, na.action = 'na.impute'
    )
  # Slight kludge...rename the column which above is given the variable name to
  # just val, to allow rbinding
  names(risk.by.variable)[2] <- 'val'
  # Append a column with the variable's name so we can distinguish this in
  # a long data frame
  risk.by.variable$var <- variable
  # Append it to our ongoing big data frame
  risk.by.variables <- rbind(risk.by.variables, risk.by.variable)
  # Save the risks as we go
  write.csv(risk.by.variables, paste0(output.filename.base, '-var-effects.csv'))
  
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