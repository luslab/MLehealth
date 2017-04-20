#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Cross-validating discretisation of input variables in a survival model
#' 
#' 

calibration.filename <- '../../output/survreg-crossvalidation-try1.csv'
comparison.filename <-
  '../../output/caliber-replicate-with-missing-var-imp-try1.csv'
# The first part of the filename for any output
output.filename.base <- '../../output/all-cv-survreg-boot-try1'


# What kind of model to fit to...currently 'cph' (Cox model), 'ranger' or
# 'rfsrc' (two implementations of random survival forests)
model.type <- 'survreg'

# If surv.vars is defined as a character vector here, the model only uses those
# variables specified, eg c('age') would build a model purely based on age. If
# not specified (ie commented out), it will use the defaults.
# surv.predict <- c('age')

#' ## Do the cross-validation
#' 
#' The steps for this are common regardless of model type, so run the script to
#' get a cross-validated model to further analyse...
#+ cox_discretised_cv, cache=cacheoption

source('../lib/all-cv-bootstrap.R', chdir = TRUE)

#' # Results
#' 
#' ## Performance
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

#' ## Cox coefficients
#'
#+ cox_coefficients_plot

cph.coeffs <-
  cphCoeffs(
    surv.model.fit.coeffs, COHORT.optimised, surv.predict,
    model.type = model.type
  )

# Create columns to hold minimum and maximum values of bins
cph.coeffs$bin.min <- NA
cph.coeffs$bin.max <- NA

for(variable in unique(cph.coeffs$var)) {
  # If it's a continuous variable, get the real centres of the bins
  if(variable %in% process.settings$var) {
    process.i <- which(variable == process.settings$var)
    
    if(process.settings$method[[process.i]] == 'binByQuantile') {
      
      variable.quantiles <-
        getQuantiles(
          COHORT.use[, variable],
          process.settings$settings[[process.i]]
        )
      # For those rows which relate to this variable, and whose level isn't
      # missing, put in the appropriate quantile boundaries for plotting
      cph.coeffs$bin.min[cph.coeffs$var == variable & 
                           cph.coeffs$level != 'missing'] <-
        variable.quantiles[1:(length(variable.quantiles) - 1)]
      cph.coeffs$bin.max[cph.coeffs$var == variable & 
                           cph.coeffs$level != 'missing'] <-
        variable.quantiles[2:length(variable.quantiles)]
      
      # Now, plot this variable as a stepped line plot using those quantile
      # boundaries
      print(
        ggplot(
          subset(cph.coeffs, var == variable),
          aes(x = bin.min, y = exp(-val))
        ) +   
          geom_step() +
          geom_step(aes(y = exp(-(val-err))), colour = 'grey') +
          geom_step(aes(y = exp(-(val+err))), colour = 'grey') +
          geom_text(aes(label = bin.min), angle = 90, hjust = 0) +
          # Missing value central estimate
          geom_hline(
            yintercept = exp(-cph.coeffs$val[cph.coeffs$var == variable & 
                                               cph.coeffs$level == 'missing']),
            colour = 'red'
          ) +
          # Missing value bounds
          geom_hline(
            yintercept = exp(-(cph.coeffs$val[cph.coeffs$var == variable & 
                                                cph.coeffs$level == 'missing'] +
                                 -cph.coeffs$err[cph.coeffs$var == variable & 
                                                   cph.coeffs$level == 'missing']
            )),
            colour = 'red'
          ) +
          geom_hline(
            yintercept = exp(-(cph.coeffs$val[cph.coeffs$var == variable & 
                                                cph.coeffs$level == 'missing'] -
                                 -cph.coeffs$err[cph.coeffs$var == variable & 
                                                   cph.coeffs$level == 'missing']
            )),
            colour = 'red'
          ) +
          ggtitle(variable)
      )
    }
  }  else {
    # If there were no instructions, it's probably a normal factor, so plot it
    # in the conventional way...
    print(
      ggplot(
        subset(cph.coeffs, var == variable),
        aes(x = factor(level), y = exp(-val))
      ) +   
        geom_bar(
          width=0.9,
          stat = "identity"
        ) +
        geom_errorbar(aes(ymin = exp(-(val + err)), ymax = exp(-(val - err)))) +
        geom_text(aes(label = level), angle = 90, hjust = 0) +
        ggtitle(variable)
    )
  }
}
