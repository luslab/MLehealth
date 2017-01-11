#' # Cross-validating discretisation of input variables
#' 
#' Random forests don't really have any tunable parameters beyond the choice of 
#' splitting criterion. However, our random forests do: because we're looking

filename <- '../../data/cohort-sample.csv'

# OK, first let's try constant numbers of bins for all variables, and also vary
# tod.rounding. Then, do some proper randomness based on the results..?
input.n.bins <- 2:10
tod.round.vals <- c(1, 0.5, 0.2, 0.1, 0.05, 0.01)
cv.n.folds <- 3
n.trees <- 50

n.data <- 1000

risk.time <- 5

source('shared.R')
require(ranger)

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
    surv.event.yes, extra.fun = caliberExtraPrep
  )

test.set <- sample(1:nrow(COHORT.prep), (1/3)*nrow(COHORT.prep))

var.combinations <-
  expand.grid(1:length(input.n.bins), 1:length(tod.round.vals))

for(i in 1:nrow(var.combinations)) {
  # prep the data given the variables provided
  COHORT.cv <-
    prepData(
      # Data for cross-validation excludes test set
      COHORT.full[-test.set, ],
      cols.keep,
      discretise.settings,
      surv.time, surv.event,
      surv.event.yes,
      default.quantiles =
        seq(
          # Quantiles are obviously between 0 and 1
          0, 1,
          # For n bins, you need n + 1 breaks
          length.out = input.n.bins[var.combinations[i, 1]] + 1
        ),
      extra.fun = caliberExtraPrep
    )
  
  # Bin the death times appropriately
  COHORT.cv$time_death_round <-
    round_any(
      COHORT.cv$time_death,
      tod.round.vals[input.n.bins[var.combinations[i, 2]]]
    )
  
  # Get folds for cross-validation
  cv.folds <- cvFolds(nrow(COHORT.cv), cv.n.folds)
  
  # Create an empty data frame to aggregate stats per fold
  cv.performance.rf <- data.frame()
  
  for(j in 1:cv.n.folds) {
    # Fit the random forest to the training set
    time.start <- handyTimer()
    COHORT.surv <- Surv(
      time  = COHORT.cv[cv.folds[[j]], 'time_death_round'],
      event = COHORT.cv[cv.folds[[j]], 'surv_event']
    )
    fit.rf <-
      ranger(
        surv.formula,
        COHORT.cv[cv.folds[[j]],], # Training set
        num.trees = n.trees,
        splitrule = 'logrank',
        num.threads = 8
      )
    time.learn.rf <- handyTimer(time.start)
    
    memory.used.rf <- as.numeric(object.size(fit.rf))
    
    # Predict the held-out validation set
    time.start <- handyTimer()
    predict.rf <- 
      predict(
        fit.rf,
        COHORT.cv[cv.folds[[j]],], # Validation set
        num.threads=8
      )
    risk.bin <- which.min(abs(predict.rf$unique.death.times - risk.time))
    # Get the chance of having died (ie 1 - survival) for all patients at that time (ie in that bin)
    rf.risk.proxy <- 1 - predict.rf$survival[, risk.bin]
    # Append C-index
    c.index.rf <-
      as.numeric(
        survConcordance(
          Surv(time_death, surv_event) ~ rf.risk.proxy,
          COHORT.cv[cv.folds[[j]],]
        )$concordance
      )
    time.predict.rf <- handyTimer(time.start)
    
    # Clear up and explicitly delete the forest just created, to be on the safe side
    rm(fit.rf)
    
    # Append the stats we've obtained from this fold
    cv.performance.rf <-
      rbind(
        cv.performance.rf,
        data.frame(
          discretise.bins = input.n.bins[var.combinations[i, 1]],
          tod.round =  input.n.bins[var.combinations[i, 2]],
          cv.fold = j,
          c.index = c.index.rf,
          time.learn = time.learn.rf,
          time.predict = time.predict.rf,
          memory.used = memory.used.rf
        )
      )
    
    # Save output at each step
    write.csv(cv.performance.rf, '../../output/rf-crossvalidation-try1.csv')
    
    # Fit the Cox model to the training set
    
    # Predict the Cox model based on the validation set
  }
  
}

# Finally, train the best obtained model on the full training set, and use it
# to predict the held-back test set