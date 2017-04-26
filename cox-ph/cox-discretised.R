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

calibration.filename <- '../../output/survreg-crossvalidation-try2.csv'
caliber.missing.coefficients.filename <-
  '../../output/caliber-replicate-with-missing-survreg-bootstrap-coeffs-1.csv'
comparison.filename <-
  '../../output/caliber-replicate-with-missing-var-imp-try2.csv'
# The first part of the filename for any output
output.filename.base <- '../../output/all-cv-survreg-boot-try2'


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

# Unpackage the uncertainties from the bootstrapped data
surv.boot.ests <- bootStats(surv.model.fit, uncertainty = '95ci')

# Save bootstrapped performance values
varsToTable(
  data.frame(
    model = 'cox',
    imputation = FALSE,
    discretised = TRUE,
    c.index = surv.boot.ests['c.test', 'val'],
    c.index.lower = surv.boot.ests['c.test', 'lower'],
    c.index.upper = surv.boot.ests['c.test', 'upper']
  ),
  performance.file
)

# Unpackage the uncertainties again, this time transformed to risks
surv.boot.ests <- 
  bootStats(surv.model.fit, uncertainty = '95ci', transform = negExp)

#' First, plot the factors and logicals as a scatter plot to compare with the
#' continuous Cox model...

# Pull coefficients from model with missing data
caliber.missing.coeffs <- read.csv(caliber.missing.coefficients.filename)

# Rename surv.boot.ests ready for merging
names(surv.boot.ests) <-
  c('cox_discrete_value', 'cox_discrete_lower', 'cox_discrete_upper')
surv.boot.ests$quantity.level <- rownames(surv.boot.ests)
# Convert variablemissing to variable_missingTRUE for compatibility
vars.with.missing <- endsWith(surv.boot.ests$quantity.level, 'missing')
surv.boot.ests$quantity.level[vars.with.missing] <-
  paste0(
    substr(
      surv.boot.ests$quantity.level[vars.with.missing],
      1,
      nchar(surv.boot.ests$quantity.level[vars.with.missing]) - nchar('missing')
    ),
    '_missingTRUE'
  )

# Create a data frame comparing them
compare.coefficients <- merge(caliber.missing.coeffs, surv.boot.ests)

ggplot(
    compare.coefficients,
    aes(x = our_value, y = cox_discrete_value, colour = unit == 'missing')
  ) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 1, colour = 'grey') +
  geom_vline(xintercept = 1, colour = 'grey') +
  geom_point() +
  geom_errorbar(aes(ymin = cox_discrete_lower, ymax = cox_discrete_upper)) +
  geom_errorbarh(aes(xmin = our_lower, xmax = our_upper)) +
  geom_text_repel(aes(label = long_name)) +
  theme_classic(base_size = 8)

# Unpack variable and level names
cph.coeffs <- cphCoeffs(
  bootStats(surv.model.fit, uncertainty = '95ci', transform = identity),
  COHORT.optimised, surv.predict, model.type = 'boot.survreg'
)

# Transform these after extraction, because the default value for the base
# level of a factor in bootStats is 0, and e^0 = 1
cph.coeffs$val <- negExp(cph.coeffs$val)
cph.coeffs$lower <- negExp(cph.coeffs$lower)
cph.coeffs$upper <- negExp(cph.coeffs$upper)

# We'll need the CALIBER scaling functions for plotting
source('caliber-scale.R')

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
      # Make the final bin the 99th percentile
      cph.coeffs$bin.max[cph.coeffs$var == variable & 
                           cph.coeffs$level != 'missing'][
                             length(variable.quantiles) - 1] <-
        quantile(COHORT.use[, variable], 0.99, na.rm = TRUE)
      
      # Add a fake data point at the highest value to finish the graph
      cph.coeffs <-
        rbind(
          cph.coeffs,
          cph.coeffs[cph.coeffs$var == variable & 
                      cph.coeffs$level != 'missing', ][
                        length(variable.quantiles) - 1, ]
          )
      # Change it so that bin.min is bin.max from the old one
      cph.coeffs$bin.min[nrow(cph.coeffs)] <-
        cph.coeffs$bin.max[cph.coeffs$var == variable & 
                            cph.coeffs$level != 'missing'][
                              length(variable.quantiles) - 1]
      
      # work out the limits on the x-axis by taking the 1st and 99th percentiles
      x.axis.limits <- quantile(COHORT.use[, variable], c(0.01, 0.99), na.rm = TRUE)
      # Use the max to provide a max value for the final bin
      
      
      # Finally, we need to scale this such that the baseline value is equal
      # to the value for the equivalent place in the Cox model, to make the
      # risks comparable...
      
      # First, we need to find the average value of this variable in the lowest
      # bin (which is always the baseline here)
      baseline.bin <- variable.quantiles[1:2]
      baseline.bin.avg <- 
        mean(
          # Take only those values of the variable which are in the range
          COHORT.use[
            inRange(COHORT.use[, variable], baseline.bin, na.false = TRUE),
            variable
          ]
        )
      # Then, scale it with the caliber scaling
      baseline.bin.val <- caliberScaleUnits(baseline.bin.avg, variable)
      # Finally, convert it to a risk
      baseline.bin.risk <- 
        caliber.missing.coeffs$our_value[
          caliber.missing.coeffs$quantity == variable
          ] ^ baseline.bin.val
      
      # And now, multiply all the discretised values by that value to make them
      # comparable...
      cph.coeffs[cph.coeffs$var == variable, c('val', 'lower', 'upper')] <-
        cph.coeffs[cph.coeffs$var == variable, c('val', 'lower', 'upper')] *
        baseline.bin.risk
      
      # Now, plot this variable as a stepped line plot using those quantile
      # boundaries
      cox.discrete.plot <-
        ggplot(
          subset(cph.coeffs, var == variable),
          aes(x = bin.min, y = val)
        ) +   
        geom_step() +
        geom_step(aes(y = lower), colour = 'grey') +
        geom_step(aes(y = upper), colour = 'grey') +
        #geom_text(aes(label = bin.min), angle = 90, hjust = 0) +
        # Missing values central estimate
        coord_cartesian(xlim = c(x.axis.limits)) +
        ggtitle(variable)
      
      # If there's a missing value risk, add it
      if(any(cph.coeffs$var == variable & cph.coeffs$level == 'missing')) {
        cox.discrete.plot <-
          cox.discrete.plot +
          geom_pointrange(
            data = cph.coeffs[cph.coeffs$var == variable & 
                                cph.coeffs$level == 'missing', ],
            aes(
              x = x.axis.limits[2], # Place at the right-hand edge of plot
              y = val,
              ymin = lower,
              ymax = upper
            ),
            colour = 'red'
          )
      }
      
      # Now, let's add the line from the continuous Cox model
      continuous.cox <-
        data.frame(
          # Use 1000 values for a smooth line
          var.x.values =
            seq(
              quantile(COHORT.use[, variable], 0.01, na.rm = TRUE),
              quantile(COHORT.use[, variable], 0.99, na.rm = TRUE),
              length.out = 1000
            )
        )
      # Scale the x-values
      continuous.cox$var.x.scaled <-
        caliberScaleUnits(continuous.cox$var.x.values, variable)
      # Use the risks to calculate risk per x for central estimate and errors
      continuous.cox$y <-
        caliber.missing.coeffs$our_value[
          caliber.missing.coeffs$quantity == variable
          ] ^ continuous.cox$var.x.scaled
      continuous.cox$upper <-
        caliber.missing.coeffs$our_upper[
          caliber.missing.coeffs$quantity == variable
          ] ^ continuous.cox$var.x.scaled
      continuous.cox$lower <-
        caliber.missing.coeffs$our_lower[
          caliber.missing.coeffs$quantity == variable
          ] ^ continuous.cox$var.x.scaled
      
      cox.discrete.plot <-
        cox.discrete.plot +
        geom_line(
          data = continuous.cox,
          aes(x = var.x.values, y = y),
          colour = 'blue'
        ) +
        geom_line(
          data = continuous.cox,
          aes(x = var.x.values, y = upper),
          colour = 'lightblue'
        ) +
        geom_line(
          data = continuous.cox,
          aes(x = var.x.values, y = lower),
          colour = 'lightblue'
        )
      
      # If there is one, add missing value risk from the continuous model
      if(any(caliber.missing.coeffs$quantity == paste0(variable, '_missing') &
             caliber.missing.coeffs$unit == 'missing')) {
        cox.discrete.plot <-
          cox.discrete.plot +
          geom_pointrange(
            data = caliber.missing.coeffs[
              caliber.missing.coeffs$quantity == paste0(variable, '_missing') &
              caliber.missing.coeffs$unit == 'missing',
            ],
            aes(
              x = x.axis.limits[2] - diff(x.axis.limits)/100, # Just before RHS
              y = our_value,
              ymin = our_lower,
              ymax = our_upper
            ),
            colour = 'blue'
          )
      }
      
      print(cox.discrete.plot)
    }
  }
}
