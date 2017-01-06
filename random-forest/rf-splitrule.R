#' # Medical records survival model parameter selection
#'
#+ setup, message=FALSE

source('prep-data.R')
require(ranger)
require(xtable)

# Number of data points to use of each of the training and test set
n.training <- 30000 # total 70096
n.test     <- 20000 # total 35048
risk.time <- 5

#tod.round.range <- c(0.01, 0.5)
tod.round <- 0.1
n.trees <- 500
splitrules <- c('C', 'logrank')
n.calibrations <- 10

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

#' ## Effect of splitting rule

#+ calibration_loop

calibrations <- data.frame()

splitrules <- rep(splitrules, each = n.calibrations)

for(splitrule in splitrules) {

  time.start <- handyTimer()
  fit.rf <-
    ranger(
      surv.formula,
      COHORT.training,
      num.trees = n.trees,
      splitrule = splitrule,
      num.threads = 8
    )
  time.learn <- handyTimer(time.start)

  forest.memory <- as.numeric(object.size(fit.rf))

  time.start <- handyTimer()
  predict.rf <- predict(fit.rf, COHORT.test, num.threads=8)
  risk.bin <- which.min(abs(predict.rf$unique.death.times - risk.time))
  # Get the chance of having died (ie 1 - survival) for all patients at that time (ie in that bin)
  rf.risk.proxy <- 1 - predict.rf$survival[, risk.bin]
  # Calculate C-index
  c.index <-
    survConcordance(
      Surv(time_death, surv_event) ~ rf.risk.proxy, COHORT.test
    )
  time.predict <- handyTimer(time.start)

  calibrations <-
    rbind(
      calibrations,
      data.frame(splitrule, time.learn, time.predict, forest.memory, c.index = as.numeric(c.index$concordance))
    )

  write.csv(calibrations, 'output/rf-splitrule.csv')

  # Clear up and explicitly delete the forest just created, to be on the safe side
  rm(fit.rf)
}

#+ prepare_plots

calibrations.summary <-
  # Needed because aggreagate creates columns which are themselves lists, this expands them
  as.data.frame(as.list(
    # Work out mean, sd and upper and lower limits for error bars
    aggregate(
      data = calibrations, . ~ splitrule,
      FUN = function(x) c(
        mean = mean(x),
        sd = sd(x),
        lower = mean(x) - sd(x),
        upper = mean(x) + sd(x)
        )
      )
    )
  )

#+ output_tables, results = 'asis'

print(
  xtable(
    calibrations
  ),
  type = 'html',
  include.rownames = FALSE
)

print(
  xtable(
    calibrations.summary
  ),
  type = 'html',
  include.rownames = FALSE
)

#' ### C-index by splitrule
#' Does splitting on C-index do better than a simple logrank split?

#+ splitrule_cindex

ggplot(calibrations.summary, aes(x = splitrule, y = c.index.mean)) +
  geom_bar(stat = "identity", aes(fill = splitrule)) +
  geom_errorbar(
    aes(ymin = c.index.lower, ymax = c.index.upper),
    width = 0.1
  )

#' ### Time to learn vs splitrule
#' The C-index split seems to take longer to do learning...but how much?

#+ splitrule_time_learn

ggplot(calibrations.summary, aes(x = splitrule, y = time.learn.mean)) +
  geom_bar(stat = "identity", aes(fill = splitrule)) +
    geom_errorbar(
      aes(ymin = time.learn.lower, ymax = time.learn.upper),
    width = 0.1
  )

#' ### Time to predict vs splitrule
#' This should surely be the same for both.

#+ splitrule_time_predict

ggplot(calibrations.summary, aes(x = splitrule, y = time.predict.mean)) +
  geom_bar(stat = "identity", aes(fill = splitrule)) +
  geom_errorbar(
    aes(ymin = time.predict.lower, ymax = time.predict.upper),
    width = 0.1
  )