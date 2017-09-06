# All model types are bootstrapped this many times
bootstraps <- 200
# n.trees is (obviously) only relevant for random forests
n.trees <- 500
# The following two variables are only relevant if the model.type is 'ranger'
split.rule <- 'logrank'
n.threads <- 10

# Cross-validation variables
input.n.bins <- 10:20
cv.n.folds <- 3
n.calibrations <- 1000
n.data <- NA # This is of full dataset...further rows may be excluded in prep

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
test.set <- testSetIndices(COHORT.prep, random.seed = 78361)

# If we've not already done a calibration, then do one
if(!file.exists(calibration.filename)) {
  # Create an empty data frame to aggregate stats per fold
  cv.performance <- data.frame()
  
  # We can parallelise this bit with foreach, so set that up
  initParallel(n.threads)
  
  # Run crossvalidations in parallel
  cv.performance <- 
    foreach(i = 1:n.calibrations, .combine = 'rbind') %dopar% {
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
    
    cv.fold.performance <- data.frame()
    
    for(j in 1:cv.n.folds) {
      time.start <- handyTimer()
      # Fit model to the training set
      surv.model.fit <-
        survivalFit(
          surv.predict,
          COHORT.cv[-cv.folds[[j]],],
          model.type = model.type,
          n.trees = n.trees,
          split.rule = split.rule
          # n.threads not used because this is run in parallel
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
      cv.fold.performance <-
        rbind(
          cv.fold.performance,
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
      
    } # End cross-validation loop (j)
    
    # rbind the performance by fold
    cv.fold.performance
  } # End calibration loop (i)
  
  # Save output at end of calibration
  write.csv(cv.performance, calibration.filename)

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

#' ## Fit the final model
#' 
#' This may take some time, so we'll cache it if possible...

#+ fit_final_model

# Fit to whole training set
surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.optimised[-test.set,], # Training set
    model.type = model.type,
    n.trees = n.trees,
    split.rule = split.rule,
    n.threads = n.threads
  )

cl <- initParallel(n.threads, backend = 'doParallel')

surv.model.params.boot <-
  foreach(
    i = 1:bootstraps,
    .combine = rbind,
    .packages = c('survival'),
    .verbose = TRUE
    ) %dopar% {
    
    # Bootstrap-sampled training set
    COHORT.boot <-
      sample.df(
        COHORT.optimised[-test.set,],
        nrow(COHORT.optimised[-test.set,]),
        replace = TRUE
      )
    
    surv.model.fit.i <-
      survivalFit(
        surv.predict,
        COHORT.boot,
        model.type = model.type,
        n.trees = n.trees,
        split.rule = split.rule,
        # 1 thread, because we're parallelising the bootstrapping
        n.threads = 1
      )
    
    # Work out other quantities of interest
    #var.imp.vector <- bootstrapVarImp(surv.model.fit.i, COHORT.boot)
    c.index <- cIndex(surv.model.fit.i, COHORT.optimised[test.set, ])
    calibration.score <-
      calibrationScoreWrapper(surv.model.fit.i, COHORT.optimised[test.set, ])
    
    data.frame(
      i,
      t(coef(surv.model.fit.i)),
      #t(var.imp.vector),
      c.index,
      calibration.score
    )
  }

# Save the fit object
write.csv(surv.model.params.boot, paste0(output.filename.base, '-surv-boot.csv'))

# Tidy up by removing the cluster
stopCluster(cl)

surv.model.fit.coeffs <-  bootStatsDf(surv.model.params.boot)
