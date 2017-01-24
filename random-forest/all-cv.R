#' # Cross-validating discretisation of input variables in a Cox model
#' 
#' 

data.filename <- '../../data/cohort-sanitised.csv'
calibration.filename <- '../../output/all-cvl-rf-try1.csv'
results.filename <- '../../output/all-cv-rf-results-try1.csv'

# What kind of model to fit to...currently 'cph' (Cox model), 'ranger' or
# 'rfsrc' (two implementations of random survival forests)
model.type <- 'ranger'

# n.trees is (obviously) only relevant for random forests
n.trees <- 500
# The following two variables are only relevant if the model.type is 'ranger'
split.rule <- 'logrank'
n.threads <- 8

# Cross-validation variables
input.n.bins <- 2:20
cv.n.folds <- 3
n.calibrations <- 100
n.data <- 30000

continuous.vars <-
  c(
    'age', 'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
  )

source('shared.R')

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- data.frame(fread(data.filename))

# If n.data was specified...
if(!is.na(n.data)){
  # Take a subset n.data in size
  COHORT.use <- sampledf(COHORT.full, n.data)
  rm(COHORT.full)
} else {
  # Use all the data
  COHORT.use <- COHORT.full
  rm(COHORT.full)
  # Without n.data, need a quick null preparation of the data to get its length
  COHORT.prep <-
    prepData(
      COHORT.use,
      cols.keep, discretise.settings, surv.time, surv.event,
      surv.event.yes, extra.fun = caliberExtraPrep, n.keep = n.data
    )
  n.data <- nrow(COHORT.prep)
}

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
          predict.vars,
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

# Fit to whole training set
surv.model.fit <-
  survivalFit(
    predict.vars,
    COHORT.cv[-cv.folds[[j]],],
    model.type = model.type,
    n.trees = n.trees,
    split.rule = split.rule,
    n.threads = n.threads
  )

# Get C-indices for training and test sets
c.index.train <- calcCoxCIndex(surv.model.fit, COHORT.optimised[-test.set, ])
c.index.test <- calcCoxCIndex(surv.model.fit, COHORT.optimised[test.set, ])

#' # Results
#' 
#' ## Performance
#' 
#' The C-index on the full training set is **`r round(c.index.train, 3)`**.
#' 
#' The C-index on the held-out test set is **`r round(c.index.test, 3)`**.
#' 
#' ## Model fit
#' 
#+ resulting_fit

print(surv.model.fit)


if(model.type == 'cph') {
  #' ## Cox coefficients
  #'
  #+ cox_coefficients_plot
  
  cph.coeffs <- cphCoeffs(fit, COHORT.optimised)
  
  ggplot(cph.coeffs, aes(x = var, y = beta, group = var)) +   
    geom_bar(aes(fill = var, group = var), width=0.9,position = position_dodge(width = 0.9), stat = "identity") +
    geom_text(aes(label = level), position = position_dodge(width = 0.9), angle = 90, hjust = 0) +
    theme(legend.position='none')
} else {
  # Something in here for random forests
}