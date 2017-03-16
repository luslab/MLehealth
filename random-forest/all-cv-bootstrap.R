#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- TRUE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Cross-validating discretisation of input variables in a survival model
#' 
#' 

data.filename <- '../../data/cohort-sanitised.csv'
calibration.filename <- '../../output/cph-crossvalidation-try1.csv'
comparison.filename <-
  '../../output/caliber-replicate-with-missing-var-imp-try1.csv'
# The first part of the filename for any output
output.filename.base <- '../../output/all-cv-cph-boot-try1'

# What kind of model to fit to...currently 'cph' (Cox model), 'ranger' or
# 'rfsrc' (two implementations of random survival forests)
model.type <- 'survreg'

# All model types are bootstrapped this many times
bootstraps <- 200
# n.trees is (obviously) only relevant for random forests
n.trees <- 500
# The following two variables are only relevant if the model.type is 'ranger'
split.rule <- 'logrank'
n.threads <- 8

# Cross-validation variables
input.n.bins <- 2:20
cv.n.folds <- 3
n.calibrations <- 100
n.data <- 30000 # This is of full dataset...further rows may be excluded in prep

# If surv.vars is defined as a character vector here, the model only uses those
# variables specified, eg c('age') would build a model purely based on age. If
# not specified (ie commented out), it will use the defaults.
# surv.predict <- c('age')

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )

source('shared.R')
require(ggrepel)

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- data.frame(fread(data.filename))

# If n.data was specified...
if(!is.na(n.data)){
  # Take a subset n.data in size
  COHORT.use <- sample.df(COHORT.full, n.data)
  rm(COHORT.full)
} else {
  # Use all the data
  COHORT.use <- COHORT.full
  rm(COHORT.full)
}

# We now need a quick null preparation of the data to get its length (some rows
# may be excluded during preparation)
COHORT.prep <-
  prepData(
    COHORT.use,
    cols.keep, discretise.settings, surv.time, surv.event,
    surv.event.yes, extra.fun = caliberExtraPrep, n.keep = n.data
  )
n.data <- nrow(COHORT.prep)

# Define indices of test set
test.set <- sample(1:n.data, (1/3)*n.data)

# If we've not already done a calibration, then do one
if(!file.exists(calibration.filename)) {
  # Create an empty data frame to aggregate stats per fold
  cv.performance <- data.frame()
  
  for(i in 1:n.calibrations) {
    cat(
      'Calibration', i, '...\n'
    )
    
    # Reset process settings with the base setings
    process.settings <-
      list(
        var        = c('anonpatid', 'time_death', 'imd_score', 'exclude'),
        method     = c(NA, NA, NA, NA),
        settings   = list(NA, NA, NA, NA)
      )
    # Generate some random numbers of bins (and for n bins, you need n + 1 breaks)
    n.bins <- sample(input.n.bins, length(continuous.vars), replace = TRUE) + 1
    names(n.bins) <- continuous.vars
    # Go through each variable setting it to bin by quantile with a random number of bins
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
    COHORT.cv <-
      prepData(
        # Data for cross-validation excludes test set
        COHORT.use[-test.set, ],
        cols.keep,
        process.settings,
        surv.time, surv.event,
        surv.event.yes,
        extra.fun = caliberExtraPrep
      )
    
    # Get folds for cross-validation
    cv.folds <- cvFolds(nrow(COHORT.cv), cv.n.folds)
  
    for(j in 1:cv.n.folds) {
      time.start <- handyTimer()
      # Fit model to the training set
      surv.model.fit <-
        survivalFit(
          surv.predict,
          COHORT.cv[-cv.folds[[j]],],
          model.type = model.type,
          n.trees = n.trees,
          split.rule = split.rule,
          n.threads = n.threads
        )
      time.learn <- handyTimer(time.start)
      
      time.start <- handyTimer()
      # Get C-indices for training and validation sets
      c.index.train <-
        cIndex(
          surv.model.fit, COHORT.cv[-cv.folds[[j]],], model.type = model.type
        )
      c.index.val <-
        cIndex(
          surv.model.fit, COHORT.cv[cv.folds[[j]],], model.type = model.type
        )
      time.predict <- handyTimer(time.start)
      
      # Append the stats we've obtained from this fold
      cv.performance <-
        rbind(
          cv.performance,
          data.frame(
            calibration = i,
            cv.fold = j,
            as.list(n.bins),
            c.index.train,
            c.index.val,
            time.learn,
            time.predict
          )
        )
      
      # Save output at each step
      write.csv(cv.performance, calibration.filename)
      
    } # End cross-validation loop (j)
    
  } # End calibration loop (i)

} else { # If we did previously calibrate, load it
  cv.performance <- read.csv(calibration.filename)
}

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

# Variable importance argument varies depending on the package being used
if(model.type == 'ranger'){
  var.imp.arg <- 'permutation'
} else if(model.type == 'rfsrc') {
  var.imp.arg <- 'permute'
} else {
  var.imp.arg <- 'NULL'
}

#' ## Fit the final model
#' 
#' This may take some time, so we'll cache it if possible...

#+ fit_final_model, cache=cacheoption

# Fit to whole training set, calculating variable importance if appropriate
surv.model.fit <-
  survivalBootstrap(
    surv.predict,
    COHORT.optimised[-test.set,], # Training set
    COHORT.optimised[test.set,],  # Test set
    model.type = model.type,
    n.trees = n.trees,
    split.rule = split.rule,
    n.threads = n.threads,
    bootstraps = bootstraps
  )

# Save the fit object
saveRDS(surv.model.fit, paste0(output.filename.base, '-surv-model.rds'))

# Get C-indices for training and test sets
surv.model.fit.coeffs <-  bootStats(surv.model.fit)

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

#' ## Individual variables
#' 
#' How we examine individual variables will vary somewhat depending on the model
#' type fitted. Let's take a look...

if(model.type == 'survreg') {
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
  

} else {
  # For random forests, take a look at the variable importance
  
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
}