bootstraps <- 200
split.rule <- 'logrank'
n.threads <- 8

# Cross-validation variables
ns.splits <- 0:20
ns.trees  <- c(500, 1000, 2000)
ns.imputations <- 1:3
cv.n.folds <- 3
n.data <- NA # This is of full dataset...further rows may be excluded in prep

calibration.filename <- paste0(output.filename.base, '-calibration.csv')

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )

untransformed.vars <- c('anonpatid', 'surv_time', 'imd_score', 'exclude')

source('../lib/shared.R')
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

# Process settings: don't touch anything!!
process.settings <-
  list(
    var       = c(untransformed.vars, continuous.vars),
    method    = rep(NA, length(untransformed.vars) + length(continuous.vars)),
    settings  = rep(NA, length(untransformed.vars) + length(continuous.vars))
  )

# If we've not already done a calibration, then do one
if(!file.exists(calibration.filename)) {
  # Create an empty data frame to aggregate stats per fold
  cv.performance <- data.frame()
  
  # Items to cross-validate over
  cv.vars <- expand.grid(ns.splits, ns.trees, ns.imputations)
  names(cv.vars) <- c('n.splits', 'n.trees', 'n.imputations')
  
  # prep the data (since we're not cross-validating on data prep this can be
  # done before the loop)
  
  # Prep the data
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
  
  # Finally, add missing flag columns, but leave the missing data intact because
  # rfsrc can do on-the-fly imputation
  COHORT.cv <- prepCoxMissing(COHORT.cv, missingReplace = NA)
  
  # Add on those column names we just created
  surv.predict <-
    c(surv.predict, names(COHORT.cv)[grepl('_missing', names(COHORT.cv))])
  
  # Run crossvalidations. No need to parallelise because rfsrc is parallelised
  for(i in 1:nrow(cv.vars)) {
    cat(
      'Calibration', i, '...\n'
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
          model.type = 'rfsrc',
          n.trees = cv.vars$n.trees[i],
          split.rule = split.rule,
          n.threads = n.threads,
          nsplit = cv.vars$n.splits[i],
          nimpute = cv.vars$n.imputations[i],
          na.action = 'na.impute'
        )
      time.learn <- handyTimer(time.start)
      
      time.start <- handyTimer()
      # Get C-indices for training and validation sets
      c.index.train <-
        cIndex(
          surv.model.fit, COHORT.cv[-cv.folds[[j]],],
          na.action = 'na.impute'
        )
      c.index.val <-
        cIndex(
          surv.model.fit, COHORT.cv[cv.folds[[j]],],
          na.action = 'na.impute'
        )
      time.predict <- handyTimer(time.start)
      
      # Append the stats we've obtained from this fold
      cv.fold.performance <-
        rbind(
          cv.fold.performance,
          data.frame(
            calibration = i,
            cv.fold = j,
            n.trees = cv.vars$n.trees[i],
            n.splits = cv.vars$n.splits[i],
            n.imputations = cv.vars$n.imputations[i],
            c.index.train,
            c.index.val,
            time.learn,
            time.predict
          )
        )
      
    } # End cross-validation loop (j)
    
    
    # rbind the performance by fold
    cv.performance <-
      rbind(
        cv.performance,
        cv.fold.performance
      )
    
    # Save output at the end of each loop
    write.csv(cv.performance, calibration.filename)
    
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

# Prep the data to fit and test with
COHORT.prep <-
  prepData(
    # Data for cross-validation excludes test set
    COHORT.use,
    cols.keep,
    process.settings,
    surv.time, surv.event,
    surv.event.yes,
    extra.fun = caliberExtraPrep
  )

# Finally, add missing flag columns, but leave the missing data intact because
# rfsrc can do on-the-fly imputation
COHORT.prep <- prepCoxMissing(COHORT.prep, missingReplace = NA)

#' ## Fit the final model
#' 
#' This may take some time, so we'll cache it if possible...

#+ fit_final_model

surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.prep[-test.set,],
    model.type = 'rfsrc',
    n.trees = cv.performance[best.calibration.row1, 'n.trees'],
    split.rule = split.rule,
    n.threads = n.threads,
    nsplit = cv.performance[best.calibration.row1, 'n.splits'],
    nimpute = cv.performance[best.calibration.row1, 'n.imputations'],
    na.action = 'na.impute'
  )

surv.model.fit.boot <-
  survivalBootstrap(
    surv.predict,
    COHORT.prep[-test.set,], # Training set
    COHORT.prep[test.set,],  # Test set
    model.type = 'rfsrc',
    n.trees = cv.performance[best.calibration.row1, 'n.trees'],
    split.rule = split.rule,
    n.threads = n.threads,
    nsplit = cv.performance[best.calibration.row1, 'n.splits'],
    nimpute = cv.performance[best.calibration.row1, 'n.imputations'],
    na.action = 'na.impute',
    bootstraps = bootstraps
  )

# Save the fit object
saveRDS(
  surv.model.fit.boot,
  paste0(output.filename.base, '-surv-model-bootstraps.rds')
)

# Get C-indices for training and test sets
surv.model.fit.coeffs <-  bootStats(surv.model.fit.boot)

# Save them to the all-models comparison table
varsToTable(
  data.frame(
    model = 'rfsrc',
    imputation = FALSE,
    discretised = FALSE,
    c.index = surv.model.fit.coeffs['c.test', 'val'],
    c.index.lower = surv.model.fit.coeffs['c.test', 'lower'],
    c.index.upper = surv.model.fit.coeffs['c.test', 'upper'],
    calibration.score = surv.model.fit.coeffs['calibration.score', 'val'],
    calibration.score.lower =
      surv.model.fit.coeffs['calibration.score', 'lower'],
    calibration.score.upper =
      surv.model.fit.coeffs['calibration.score', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)