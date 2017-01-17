#' # Cross-validating discretisation of input variables in a Cox model
#' 
#' 

filename <- '../../data/cohort-sanitised.csv'

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

# Define a survival model formula based on variables being used
surv.formula <-
  formula(
    paste0(
      'COHORT.surv~',
      paste(surv.predict, collapse='+')
    )
  )

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- data.frame(fread(filename))

# We need to do a quick null preparation of the data to get its length
COHORT.prep <-
  prepData(
    COHORT.full,
    cols.keep, discretise.settings, surv.time, surv.event,
    surv.event.yes, extra.fun = caliberExtraPrep, n.keep = n.data
  )

test.set <- sample(1:nrow(COHORT.prep), (1/3)*nrow(COHORT.prep))

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
      COHORT.full[-test.set, ],
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
    # Create survival object for training set
    COHORT.surv <- Surv(
      time  = COHORT.cv[-cv.folds[[j]], 'time_death'],
      event = COHORT.cv[-cv.folds[[j]], 'surv_event']
    )
    # Fit the Cox model to the training set
    fit.cph <- cph(surv.formula, COHORT.cv[-cv.folds[[j]],], surv = TRUE)
    time.learn.cph <- handyTimer(time.start)
    
    time.start <- handyTimer()
    # Get C-indices for training and validation sets
    c.index.cph.train <- calcCoxCIndex(fit.cph, COHORT.cv[-cv.folds[[j]],])
    c.index.cph.val <- calcCoxCIndex(fit.cph, COHORT.cv[cv.folds[[j]],])
    time.predict.cph <- handyTimer(time.start)
    
    # Append the stats we've obtained from this fold
    cv.performance <-
      rbind(
        cv.performance,
        data.frame(
          calibration = i,
          cv.fold = j,
          as.list(n.bins),
          c.index.cph.train,
          c.index.cph.val,
          time.learn.cph,
          time.predict.cph
        )
      )
    
    # Save output at each step
    write.csv(cv.performance, '../../output/cph-crossvalidation-try1.csv')
  }
  
}

# Finally, train the best obtained model on the full training set, and use it
# to predict the held-back test set