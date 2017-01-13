#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy=FALSE)

#' # Medical records survival model comparison
#'
#' This document compares two different kinds of survival model to one-another:
#' a conventional Cox proportional hazards model, and a random survival forest.
#' 
#' In both cases, the baseline is the C-index calculated in
#' [Rapsomaniki _et al._ 2013](http://eurheartj.oxfordjournals.org/content/early/2013/12/17/eurheartj.eht533),
#' which used this dataset to build a prognostic model. The C-index on their
#' internal validation set was **0.811**.
#'
#' ## Data preparation
#'
#' We start by preparing the data for consistent comparisons. Since the data
#' contain missing values and we wish to preserve the information contained in
#' this missingness, the simplest way to do this is to discretise, with
#' the missing values as a separate category. By contrast, Rapsomaniki _et al._
#' used imputation to fill in missing values based on values in other patients.

#+ setup, message=FALSE

source('prep-data.R')

# Number of data points to use of each of the training and test set
n.training <- 30000 # total 70096
n.test     <- 20000 # total 35048

# Subset training and test data frames to only use number of patients specified
COHORT.training <- COHORT.training.full[1:n.training,]
COHORT.test <- COHORT.test.full[1:n.test,]

#' ## Discretising variables
#' This is a version of the script testing discretisation. The quantiles used
#' were `r print(discretise.quantiles)`.

#' ## Removing variables which cause an error
#'
#' For some reason, including the variable ``hx_mi`` throws the error
#' ``X matrix deemed to be singular; variable hx_mi=TRUE``. This is usually
#' caused by a variable which is perfectly correlated with the output you're
#' trying to model (see
#' [here](http://stackoverflow.com/questions/20977401/coxph-x-matrix-deemed-to-be-singular)),
#' which clearly isn't the case here (using the ``table`` function with
#' ``hx_mi`` and a variety of other variables yields no perfect correlations).
#' However, for now, I've excluded it from the remainder of the analysis.

#+ exclude_vars

exclude.vars <- c('hx_mi')
surv.predict <- surv.predict[!(surv.predict %in% exclude.vars)]

#' ## Cox proportional hazards model
#'
#' Start by fitting a conventional Cox model to the dataset...

#+ cox_model

# Make a survival object
COHORT.survival <- Surv(
  time  = COHORT.training[, c(surv.time)],
  event = COHORT.training[, 'surv_event']
)

fit.cph <- cph(
  # Build a formula for the survival model from the variables chosen
  formula(
    paste0(
      'COHORT.survival~',
      paste(surv.predict, collapse='+')
    )
  ),
  data = COHORT.training,
  surv = TRUE
)

#' Then, make some predictions on the test set and use them to work out the
#' C-index, to give a measure of performance...

#+ cox_c-index

predict.cph <- predict(fit.cph, COHORT.test)

c.index.cph <-
  survConcordance(
    Surv(time_death, surv_event) ~ predict.cph,
    COHORT.test
  )

print(c.index.cph)

#' Our ourput value of **`r round(c.index.cph$concordance, 3)`**
#' compares reasonably favourably to the 0.811 obtained by Rapsomaniki
#' _et al._, especially considering we've not done any complex per-variable
#' calibration as they did, nor even any imputation!
#'
#' Let's take a quick look at the fitted coefficients, just to make sure that
#' they make sense...

#+ cox_survival_curves_plot

# Take 100 random patients
random.patients <- sample(1:nrow(COHORT.test), 100)

# Get their survival curves
surv.curves.cph <- cphSurvivalCurves(COHORT.test[random.patients,], fit.cph)

ggplot(surv.curves.cph,
       aes(
         x = t,
         y = s,
         colour = time_death,
         group = id,
         alpha = surv_event
       )
) +
  geom_line()

#' ## Random survival forests
#' 
#' First, we round times of death for acceptable random forest performance...
#' 
#+ rf_setup, message=FALSE

n.trees <- 1000

tod.round <- 0.1

COHORT.training$time_death_round <-
  round_any(COHORT.training$time_death, tod.round)
COHORT.test$time_death_round <- round_any(COHORT.test$time_death, tod.round)

COHORT.survival.round <- Surv(
  time  = COHORT.training[, 'time_death_round'],
  event = COHORT.training[, 'surv_event']
)

#' 
#' ### Ranger
#'
#' Let's try using random survival forests from the R package ``ranger``.
#'
#+ ranger_setup, message=FALSE

require(ranger)

#+ ranger_learn, cache=cacheoption

time.start <- handyTimer()

fit.rf <-
  ranger(
    formula(
      paste0(
        'COHORT.survival.round~',
        paste(surv.predict, collapse='+')
      )
    ),
    COHORT.training,
    num.trees = n.trees,
    #splitrule = 'C',
    num.threads = 8
  )

time.taken <- handyTimer(time.start)

#' Training the forest is OK-fast but not instant (it took
#' `r round(time.taken/60, 1)` minutes). Let's see how we did...

#+ rf_predict, cache=cacheoption

time.start <- handyTimer()
c.index <- calcRFCIndex(fit.rf, COHORT.test, 5)
time.taken <- handyTimer(time.start)

#' *C = `r round(c.index, 3)`*
#' 
#' Calculating C-index took `r round(time.taken/60, 1)` minutes.

#+ rf_survival_curves_plot

predict.rf <- predict(fit.rf, COHORT.test[random.patients,], num.threads=8)

# Get their survival curves
surv.curves.rf <- rfSurvivalCurves(COHORT.test[random.patients,], predict.rf)

ggplot(surv.curves.rf,
       aes(
         x = t,
         y = s,
         colour = time_death,
         group = id,
         alpha = surv_event
       )
) +
  geom_line()

#' ### randomForestSRC
#'
#' How does a different package perform? Let's try ``randomForestSRC``...

#+ rfsrc_setup, message=FALSE

require(randomForestSRC)

#+ rfsrc_learn, cache=cacheoption

time.start <- handyTimer()

fit.rf <-
  rfsrc(
    formula(
      paste0(
        'COHORT.survival.round~',
        paste(surv.predict, collapse='+')
      )
    ),
    COHORT.training,
    ntree = n.trees
  )

time.taken <- handyTimer(time.start)

#' Training took `r round(time.taken/60, 1)` minutes). Let's see how we did...

#+ rfsrc_predict, cache=cacheoption

time.start <- handyTimer()
c.index <- calcRFCIndex(fit.rf, COHORT.test, 5)
time.taken <- handyTimer(time.start)

#' *C = `r round(c.index, 3)`*
#' 
#' Calculating C-index took `r round(time.taken/60, 1)` minutes.

#+ rfsrc_survival_curves_plot

predict.rf <- predict(fit.rf, COHORT.test[random.patients,], num.threads=8)

# Get their survival curves
surv.curves.rf <- rfSurvivalCurves(COHORT.test[random.patients,], predict.rf)

ggplot(surv.curves.rf,
       aes(
         x = t,
         y = s,
         colour = time_death,
         group = id,
         alpha = surv_event
       )
) +
  geom_line()