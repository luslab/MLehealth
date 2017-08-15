source('../lib/shared.R')

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )

# Load in the discretised Cox model for plotting
surv.model.fit.boot <- readRDS('../../output/all-cv-survreg-boot-try5-surv-model.rds')

# Pull coefficients from model with missing data
caliber.missing.coefficients.filename <-
  '../../output/caliber-replicate-with-missing-survreg-bootstrap-coeffs-1.csv'
caliber.missing.coeffs <- read.csv(caliber.missing.coefficients.filename)
# Log them to get them on the same scale as discrete model
caliber.missing.coeffs$our_value <- -log(caliber.missing.coeffs$our_value)
caliber.missing.coeffs$our_lower <- -log(caliber.missing.coeffs$our_lower)
caliber.missing.coeffs$our_upper <- -log(caliber.missing.coeffs$our_upper)

# Load the data
COHORT.use <- data.frame(fread(data.filename))

# Open the calibration to find the best binning scheme
calibration.filename <- '../../output/survreg-crossvalidation-try5.csv'
cv.performance <- read.csv(calibration.filename)

# Find the best calibration...
# First, average performance across cross-validation folds
cv.performance.average <-
  aggregate(
    c.index.val ~ calibration,
    data = cv.performance,
    mean
  )
# Find the highest value
best.calibration <-
  cv.performance.average$calibration[
    which.max(cv.performance.average$c.index.val)
    ]
# And finally, find the first row of that calibration to get the n.bins values
best.calibration.row1 <-
  min(which(cv.performance$calibration == best.calibration))

# Get its parameters
n.bins <-
  t(
    cv.performance[best.calibration.row1, continuous.vars]
  )

# Prepare the data with those settings...

# Reset process settings with the base setings
process.settings <-
  list(
    var        = c('anonpatid', 'time_death', 'imd_score', 'exclude'),
    method     = c(NA, NA, NA, NA),
    settings   = list(NA, NA, NA, NA)
  )
for(j in 1:length(continuous.vars)) {
  process.settings$var <- c(process.settings$var, continuous.vars[j])
  process.settings$method <- c(process.settings$method, 'binByQuantile')
  process.settings$settings <-
    c(
      process.settings$settings,
      list(
        seq(
          # Quantiles are obviously between 0 and 1
          0, 1,
          # Choose a random number of bins (and for n bins, you need n + 1 breaks)
          length.out = n.bins[j]
        )
      )
    )
}

# prep the data given the variables provided
COHORT.optimised <-
  prepData(
    # Data for cross-validation excludes test set
    COHORT.use,
    cols.keep,
    process.settings,
    surv.time, surv.event,
    surv.event.yes,
    extra.fun = caliberExtraPrep
  )

# Unpack variable and level names
cph.coeffs <- cphCoeffs(
  bootStats(surv.model.fit.boot, uncertainty = '95ci', transform = `-`),
  COHORT.optimised, surv.predict, model.type = 'boot.survreg'
)

# We'll need the CALIBER scaling functions for plotting
source('../cox-ph/caliber-scale.R')

# set up a list to store the plotted plots
cox.discrete.plots <- list()

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
      
      cox.discrete.plots[[variable]] <- cox.discrete.plot
    }
  }
}

ggplot(subset(cph.coeffs, var == 'age'), aes(x = bin.max, y = val)) +
  geom_step() +
  geom_line(data = df, aes(x, y/10))
