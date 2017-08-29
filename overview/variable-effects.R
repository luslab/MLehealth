cox.disc.filename <- '../../output/all-cv-survreg-boot-try5-surv-model.rds'
caliber.missing.coefficients.filename <-
  '../../output/caliber-replicate-with-missing-survreg-6-linear-age-coeffs-3.csv'
rf.filename <- '../../output/rfsrc-cv-nsplit-try3-var-effects.csv'

source('../lib/shared.R')
requirePlus('cowplot')

# Amount of padding at the right-hand side to make space for missing values
missing.padding <- 0.05

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )

# Load in the discretised Cox model for plotting
surv.model.fit.boot <- readRDS(cox.disc.filename)

# Pull coefficients from model with missing data
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

# set up list to store the plots
cox.discrete.plots <- list()
# Add dummy columns for x-position of missing values
caliber.missing.coeffs$missing.x.pos.cont <- NA
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
        labs(x = variable, y = 'Bx')
      
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
        
        # Put this missing value 2/3rds of the way into the missing area
        caliber.missing.coeffs$missing.x.pos.cont[
          caliber.missing.coeffs$quantity == paste0(variable, '_missing') &
            caliber.missing.coeffs$unit == 'missing'] <-
                x.axis.limits[2] + 2 * diff(x.data.range) * missing.padding / 3
        
        cox.discrete.plot <-
          cox.discrete.plot +
          geom_pointrange(
            data = caliber.missing.coeffs[
              caliber.missing.coeffs$quantity == paste0(variable, '_missing') &
                caliber.missing.coeffs$unit == 'missing',
              ],
            aes(
              x = missing.x.pos.cont,
              y = our_value, ymin = our_lower, ymax = our_upper
            ),
            colour = 'blue'
          )
      }
      
      # Finally, set the x-axis limits; will just be the data range, or data
      # range plus a bit if there are missing values to squeeze in
      cox.discrete.plot <-
        cox.discrete.plot +
        coord_cartesian(xlim = x.axis.limits) +
        theme(axis.title.y = element_blank()) +
        theme(plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "cm"))
      
      cox.discrete.plots[[variable]] <- cox.discrete.plot
    }
  }
}

# Load the random forest variable effects file
risk.by.variables <- read.csv(rf.filename)
rf.vareff.plots <- list()

for(variable in unique(risk.by.variables$var)) {
  # Get the mean of the normalised risk for every value of the variable
  risk.aggregated <-
    aggregate(
      as.formula(paste0('risk.normalised ~ val')),
      subset(risk.by.variables, var == variable), median
    )
  
  # work out the limits on the axes by taking the 1st and 99th percentiles
  x.axis.limits <-
    quantile(COHORT.use[, variable], c(0.01, 0.99), na.rm = TRUE)
  y.axis.limits <-
    quantile(subset(risk.by.variables, var == variable)$risk.normalised, c(0.05, 0.95), na.rm = TRUE)
  
  # If there's a missing value risk in the graph above, expand the axes so they
  # match
  if(any(cph.coeffs$var == variable & cph.coeffs$level == 'missing')) {
    x.axis.limits[2] <- 
      x.axis.limits[2] + diff(x.data.range) * missing.padding
  }
  
  rf.vareff.plots[[variable]] <-
    ggplot(
      subset(risk.by.variables, var == variable), 
      aes(x = val, y = log(risk.normalised))
    ) +
    geom_line(alpha=0.003, aes(group = id)) +
    geom_line(data = risk.aggregated, colour = 'blue') +
    coord_cartesian(xlim = x.axis.limits, ylim = log(y.axis.limits)) +
    labs(x = variable) +
    theme(
      plot.margin = unit(c(0.2, 0.1, 0.2, 0.1), "cm"),
      axis.title.y = element_blank()
    )
}


plot_grid(
  cox.discrete.plots[['age']],
  cox.discrete.plots[['haemoglobin_6mo']],
  cox.discrete.plots[['total_wbc_6mo']],
  cox.discrete.plots[['crea_6mo']],
  rf.vareff.plots[['age']],
  rf.vareff.plots[['haemoglobin_6mo']],
  rf.vareff.plots[['total_wbc_6mo']],
  rf.vareff.plots[['crea_6mo']],
  labels = c('A', rep('', 3), 'B', rep('', 3)),
  align = "v", ncol = 4
)

ggsave(
  '../../output/variable-effects.pdf',
  width = 16,
  height = 10,
  units = 'cm',
  useDingbats = FALSE
)
