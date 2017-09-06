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

calibration.filename <- '../../output/survreg-crossvalidation-try5.csv'
caliber.missing.coefficients.filename <-
  '../../output/caliber-replicate-with-missing-survreg-bootstrap-coeffs-1.csv'
comparison.filename <-
  '../../output/caliber-replicate-with-missing-var-imp-try2.csv'
# The first part of the filename for any output
output.filename.base <- '../../output/all-cv-survreg-boot-try5'


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
#' ### C-index
#' 
#' C-index is **`r round(surv.model.fit.coeffs['c.index', 'val'], 3)`
#' (`r round(surv.model.fit.coeffs['c.index', 'lower'], 3)` - 
#' `r round(surv.model.fit.coeffs['c.index', 'upper'], 3)`)** on the held-out
#' test set.
#' 
#'
#' ### Calibration
#' 
#' The bootstrapped calibration score is
#' **`r round(surv.model.fit.coeffs['calibration.score', 'val'], 3)`
#' (`r round(surv.model.fit.coeffs['calibration.score', 'lower'], 3)` - 
#' `r round(surv.model.fit.coeffs['calibration.score', 'upper'], 3)`)**.
#' 
#' Let's draw a representative curve from the unbootstrapped fit... (It would be
#' better to draw all the curves from the bootstrap fit to get an idea of
#' variability, but I've not implemented this yet.)
#' 
#+ calibration_plot

calibration.table <-
  calibrationTable(surv.model.fit, COHORT.optimised[test.set, ])

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' 
#' ## Model fit
#'
#+ resulting_fit

print(surv.model.fit)

#' ## Cox coefficients
#'
#+ cox_coefficients_plot


# Save bootstrapped performance values
varsToTable(
  data.frame(
    model = 'cox',
    imputation = FALSE,
    discretised = TRUE,
    c.index = surv.model.fit.coeffs['c.index', 'val'],
    c.index.lower = surv.model.fit.coeffs['c.index', 'lower'],
    c.index.upper = surv.model.fit.coeffs['c.index', 'upper'],
    calibration.score = surv.model.fit.coeffs['calibration.score', 'val'],
    calibration.score.lower =
      surv.model.fit.coeffs['calibration.score', 'lower'],
    calibration.score.upper =
      surv.model.fit.coeffs['calibration.score', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)

# Unpackage the uncertainties again, this time transformed because survreg
# returns negative values
surv.boot.ests <- bootStatsDf(surv.model.params.boot, transform = `-`)

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
  bootStats(surv.model.fit.boot, uncertainty = '95ci', transform = `-`),
  COHORT.optimised, surv.predict, model.type = 'boot.survreg'
)

# We'll need the CALIBER scaling functions for plotting
source('../cox-ph/caliber-scale.R')

# set up list to store the plots
cox.discrete.plots <- list()
# Add dummy columns for x-position of missing values
cph.coeffs$missing.x.pos.cont <- NA
cph.coeffs$missing.x.pos.disc <- NA

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
      
      # Work out data range by taking the 1st and 99th percentiles
      # Use the max to provide a max value for the final bin
      # Also use for x-axis limits, unless there are missing values to
      # accommodate on the right-hand edge.
      x.data.range <-
        quantile(COHORT.use[, variable], c(0.01, 0.99), na.rm = TRUE)
      x.axis.limits <- x.data.range
      
      
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
      baseline.bin.val <-
        caliberScaleUnits(baseline.bin.avg, variable) * 
        caliber.missing.coeffs$our_value[
          caliber.missing.coeffs$quantity == variable
          ]
      
      # And now, add all the discretised values to that value to make them
      # comparable...
      cph.coeffs[cph.coeffs$var == variable, c('val', 'lower', 'upper')] <-
        cph.coeffs[cph.coeffs$var == variable, c('val', 'lower', 'upper')] -
        baseline.bin.val
      
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
        ggtitle(variable)
      
      # If there's a missing value risk, add it
      if(any(cph.coeffs$var == variable & cph.coeffs$level == 'missing')) {
        # Expand the x-axis to squeeze the missing values in
        x.axis.limits[2] <- 
          x.axis.limits[2] + diff(x.data.range) * missing.padding
        # Put this missing value a third of the way into the missing area
        cph.coeffs$missing.x.pos.disc[
          cph.coeffs$var == variable &
            cph.coeffs$level == 'missing'] <-
          x.axis.limits[2] + diff(x.data.range) * missing.padding / 3
        
        # Add the point to the graph (we'll set axis limits later)
        cox.discrete.plot <-
          cox.discrete.plot +
          geom_pointrange(
            data = cph.coeffs[cph.coeffs$var == variable & 
                                cph.coeffs$level == 'missing', ],
            aes(
              x = missing.x.pos.disc,
              y = val, ymin = lower,
              ymax = upper
            ),
            colour = 'red'
          )
      }
      
      # Now, let's add the line from the continuous Cox model. We only need two
      # points because the lines are straight!
      continuous.cox <-
        data.frame(
          var.x.values = x.data.range
        )
      # Scale the x-values
      continuous.cox$var.x.scaled <-
        caliberScaleUnits(continuous.cox$var.x.values, variable)
      # Use the risks to calculate risk per x for central estimate and errors
      continuous.cox$y <-
        -caliber.missing.coeffs$our_value[
          caliber.missing.coeffs$quantity == variable
          ] * continuous.cox$var.x.scaled
      continuous.cox$upper <-
        -caliber.missing.coeffs$our_upper[
          caliber.missing.coeffs$quantity == variable
          ] * continuous.cox$var.x.scaled
      continuous.cox$lower <-
        -caliber.missing.coeffs$our_lower[
          caliber.missing.coeffs$quantity == variable
          ] * continuous.cox$var.x.scaled
      
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
        # Expand the x-axis to squeeze the missing values in
        x.axis.limits[2] <- 
          x.axis.limits[2] + diff(x.data.range) * missing.padding
        # Put this missing value 2/3rds of the way into the missing area
        cph.coeffs$missing.x.pos.cont[
          cph.coeffs$var == variable &
            cph.coeffs$level == 'missing'] <-
          x.axis.limits[2] + diff(x.data.range) * missing.padding / 3
        x.axis.limits[2] + 2 * diff(x.data.range) * missing.padding / 3
        
        cox.discrete.plot <-
          cox.discrete.plot +
          geom_pointrange(
            data = cph.coeffs[
              cph.coeffs$var == variable &
                cph.coeffs$level == 'missing',
              ],
            aes(
              x = missing.x.pos.cont,
              y = val, ymin = lower, ymax = upper
            ),
            colour = 'blue'
          )
      }
      
      # Finally, set the x-axis limits; will just be the data range, or data
      # range plus a bit if there are missing values to squeeze in
      cox.discrete.plot <-
        cox.discrete.plot +
        coord_cartesian(xlim = x.axis.limits)
      
      cox.discrete.plots[[variable]] <- cox.discrete.plot
    }
  }
}

print(cox.discrete.plots)
