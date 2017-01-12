#' # Medical records survival model parameter selection
#'
#+ setup, message=FALSE

source('prep-data.R')
require(ranger)

# Number of data points to use of each of the training and test set
n.training <- 30000 # total 70096
n.test     <- 20000 # total 35048
risk.time <- 5

#tod.round.range <- c(0.01, 0.5)
tod.round <- 0.1
max.trees <- 5000
n.calibrations <- 30

# Subset training and test data frames to only use number of patients specified
COHORT.training <- COHORT.training.full[1:n.training,]
COHORT.test <- COHORT.test.full[1:n.test,]

COHORT.training$time_death_round <-
  round_any(COHORT.training$time_death, tod.round)
COHORT.test$time_death_round <- round_any(COHORT.test$time_death, tod.round)

COHORT.survival.round <- Surv(
  time  = COHORT.training[, 'time_death_round'],
  event = COHORT.training[, 'surv_event']
)

surv.formula <-
  formula(
    paste0(
      'COHORT.survival.round~',
      paste(surv.predict, collapse='+')
    )
  )

calibrations <- data.frame()

#' ## Number of survival trees needed

#+ calibration_loop

for(i in 1:n.calibrations) {
  n.trees <-
    sample(
      1:max.trees, 1,
      # weight towards smaller trees
      prob = 1/sqrt(1:max.trees)
    )

  time.start <- proc.time()['elapsed']
  fit.rf <-
    ranger(
      surv.formula,
      COHORT.training,
      num.trees = n.trees,
      #splitrule = 'C',
      num.threads = 8
    )
  time.end <- proc.time()['elapsed']

  time.learn <- time.end - time.start

  forest.memory <- as.numeric(object.size(fit.rf))

  time.start <- proc.time()['elapsed']
  predict.rf <- predict(fit.rf, COHORT.test, num.threads=8)
  risk.bin <- which.min(abs(predict.rf$unique.death.times - risk.time))
  # Get the chance of having died (ie 1 - survival) for all patients at that time (ie in that bin)
  rf.risk.proxy <- 1 - predict.rf$survival[, risk.bin]
  # Calculate C-index
  c.index <-
    survConcordance(
      Surv(time_death, surv_event) ~ rf.risk.proxy, COHORT.test
    )
  time.end <- proc.time()['elapsed']
  time.predict <- time.end - time.start

  calibrations <-
    rbind(
      calibrations,
      data.frame(n.trees, time.learn, time.predict, forest.memory, c.index = c.index$concordance)
    )

  # Save output at each step
  write.csv(calibrations, '../output/rf-ntrees.csv')
  # Explicitly delete the forest just created, to be on the safe side
  rm(fit.rf)
}


#' ### Processing time vs number of trees
#' See how performance scales with size of forest...these should both be linear.

#+ compare_times_to_ntrees

ggplot(calibrations, aes(x = n.trees)) +
  geom_point(aes(y = time.learn), colour = 'red') +
  geom_point(aes(y = time.predict), colour = 'blue')

#' ### Time to learn vs time to predict
#' (Expecting this to be basically linear!)

#+ compare_learning_prediction

ggplot(calibrations, aes(time.learn, time.predict)) +
  geom_point()

#' ### Number of trees vs performance
#' Let's see where this saturates...

#+ compare_ntrees_cindex

ggplot(calibrations, aes(n.trees, c.index)) +
  geom_point()

#' ### Time to learn vs performance
#' And what's the actual trade-off between time spent and performance?

#+ compare_time_cindex

ggplot(calibrations, aes(time.learn, c.index)) +
  geom_point()

#' ### Forest size vs memory usage
#' How does memory usage by the model relate to forest size?

#+ compare_ntrees_memory

ggplot(calibrations, aes(n.trees, forest.memory)) +
  geom_point()

#' ### Memory usage vs performance
#' And finally, what's the trade-off between memory used by the model and performance?

#+ compare_memory_cindex

ggplot(calibrations, aes(forest.memory, c.index)) +
  geom_point()