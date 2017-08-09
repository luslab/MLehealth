bootstraps <- 3
split.rule <- 'logrank'
n.threads <- 20

# Cross-validation variables
ns.splits <- c(0, 5, 10, 15, 20, 30)
ms.try <- c(50, 100, 200, 300, 400)
n.trees.cv  <- 500
n.imputations <- 3
cv.n.folds <- 3
n.trees.final <- 2000
n.data <- NA # This is of full dataset...further rows may be excluded in prep

calibration.filename <- paste0(output.filename.base, '-calibration.csv')

# If we've not already done a calibration, then do one
if(!file.exists(calibration.filename)) {
  # Create an empty data frame to aggregate stats per fold
  cv.performance <- data.frame()
  
  # Items to cross-validate over
  cv.vars <- expand.grid(ns.splits, ms.try)
  names(cv.vars) <- c('n.splits', 'm.try')
  
  COHORT.cv <- COHORT.bigdata[-test.set, ]
  
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
          n.trees = n.trees.cv,
          split.rule = split.rule,
          n.threads = n.threads,
          nsplit = cv.vars$n.splits[i],
          nimpute = n.imputations,
          na.action = 'na.impute',
          mtry = cv.vars$m.try[i]
        )
      time.learn <- handyTimer(time.start)
      
      time.start <- handyTimer()
      # Get C-index on validation set
      c.index.val <-
        cIndex(
          surv.model.fit, COHORT.cv[cv.folds[[j]],],
          na.action = 'na.impute'
        )
      time.c.index <- handyTimer(time.start)
      
      time.start <- handyTimer()
      # Get C-index on validation set
      calibration.score <-
        calibrationScore(
          calibrationTable(
            surv.model.fit, COHORT.cv[cv.folds[[j]],], na.action = 'na.impute'
          )
        )
      time.calibration <- handyTimer(time.start)
      
      # Append the stats we've obtained from this fold
      cv.fold.performance <-
        rbind(
          cv.fold.performance,
          data.frame(
            calibration = i,
            cv.fold = j,
            n.splits = cv.vars$n.splits[i],
            m.try = cv.vars$m.try[i],
            c.index.val,
            time.learn,
            time.c.index,
            time.calibration
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

#' ## Fit the final model
#' 
#' This may take some time, so we'll cache it if possible...

#+ fit_final_model

surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.bigdata[-test.set,],
    model.type = 'rfsrc',
    n.trees = n.trees.final,
    split.rule = split.rule,
    n.threads = n.threads,
    nimpute = n.imputations,
    nsplit = cv.performance[best.calibration.row1, 'n.splits'],
    mtry = cv.performance[best.calibration.row1, 'm.try'],
    na.action = 'na.impute',
    importance = 'permute'
  )

# Save the fit object
saveRDS(
  surv.model.fit,
  paste0(output.filename.base, '-surv-model.rds')
)

surv.model.fit.boot <-
  survivalBootstrap(
    surv.predict,
    COHORT.bigdata[-test.set,], # Training set
    COHORT.bigdata[test.set,],  # Test set
    model.type = 'rfsrc',
    n.trees = n.trees.final,
    split.rule = split.rule,
    n.threads = n.threads,
    nimpute = n.imputations,
    nsplit = cv.performance[best.calibration.row1, 'n.splits'],
    mtry = cv.performance[best.calibration.row1, 'm.try'],
    na.action = 'na.impute',
    bootstraps = bootstraps
  )

# Save the fit object
saveRDS(
  surv.model.fit.boot,
  paste0(output.filename.base, '-surv-model-bootstraps.rds')
)

# Get C-indices for training and test sets
surv.model.fit.coeffs <- bootStats(surv.model.fit.boot, uncertainty = '95ci')

# Save them to the all-models comparison table
varsToTable(
  data.frame(
    model = 'rfbigdata',
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