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

#+ cox_coefficients_plot

cph.coeffs <- cphCoeffs(fit.cph, COHORT.test)

ggplot(cph.coeffs, aes(x = var, y = beta, group = var)) +   
  geom_bar(aes(fill = var, group = var), width=0.9,position = position_dodge(width = 0.9), stat = "identity") +
  geom_text(aes(label = level), position = position_dodge(width = 0.9), angle = 90, hjust = 0) +
  theme(legend.position='none')

#' Finally, let's plot a few survival curves as a sanity check...

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


#' OK, so it hopefully looks from those like patients known to die early (the
#' darker, opaque lines) have worse predicted survival curves than patients
#' known to die late (lighter, opaque lines). Semi-transparent lines represent
#' censored patients and are coloured by their censorship time, so are included
#' for information only.
#'
#' It seems that the Cox model is doing sensible things and performing
#' reasonably well. Time to try random forests...

#' ## Random survival forests
#'
#' Let's try using random survival forests from the R package ``ranger``. We
#' start by including this package and defining the number of trees we're going
#' to be using...

#+ rf_setup, message=FALSE

require(ranger)
n.trees <- 1000

#' We're also going to define a new column called ``time_death_round`` which
#' does exactly what it says on the tin, rounding the time of death to the
#' nearest 0.1 (in this case). This is because a survival curve is defined at
#' every unique time of death present in the dataset: every time someone dies,
#' the Kaplan-Meier estimate drops by a tiny bit. Given that the times of death
#' are reported to very high accuracy, this means we have
#' `r length(unique(COHORT.training$time_death))` unique values in our training
#' set, which leads to hilariously prohibitive processing times and memory usage
#' (I've seen usage during training exceed 100 GB of RAM before implementing
#' this!) because every branch of every tree stores of order that many values.

#+ rf_tod_round

tod.round <- 0.1

COHORT.training$time_death_round <-
  round_any(COHORT.training$time_death, tod.round)
COHORT.test$time_death_round <- round_any(COHORT.test$time_death, tod.round)

#' We now have just `r length(unique(COHORT.training$time_death_round))` unique
#' values. Much better!
#'
#' Now, to learn that forest...

#+ rf_learn, cache=cacheoption

COHORT.survival.round <- Surv(
  time  = COHORT.training[, 'time_death_round'],
  event = COHORT.training[, 'surv_event']
)

time.start <- proc.time()['elapsed']

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

time.end <- proc.time()['elapsed']

#' Training the forest is OK-fast but not instant (it took
#' `r round((time.end - time.start)/60, 1)` minutes). Let's see how we did...

#+ rf_predict, cache=cacheoption

time.start <- proc.time()['elapsed']
predict.rf <- predict(fit.rf, COHORT.test, num.threads=8)
time.end <- proc.time()['elapsed']

#' Prediction took `r round((time.end - time.start)/60, 1)` minutes. Now to
#' calculate the concordance score...except it's not so simple as with a Cox
#' model...

#' ### Calculating concordance scores with random forests
#'
#' The
#' [survConcordance function](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/survConcordance.html)
#' from the survival package in R allows concordance scores to be calculated.
#' Its usage is commonly something like:
#'
#' ``survConcordance(Surv(time_death, surv_event) ~ predicted.risk, input.data.frame)``
#'
#' This creates a survival object from columns in ``input.data.frame`` and
#' compares that to scores, ``predicted.risk``. What's slightly confusing is
#' that the risk here is some scalar measure which, because it's a risk, goes up
#' as predicted lifespan goes down. Thus, if you input something like lifespan,
#' you'll actually get a result which is (1 - _c_) because pairs which should be
#' concordant are discordant, and vice-versa.
#'
#' It shouldn't (and doesn't seem to) matter what kind of number you input here,
#' because concordance simply checks whether one number is bigger than another.
#' Thus, anything where a bigger number is bad for the patient works.
#' The question then arises, in a nonparametric model, what to input. In
#' something like a Cox model, the multiplicative risk compared to the baseline
#' hazard has a known value. In a nonparametric model, you can imagine weird
#' survival curves which aren't monotonically worse or better, eg one which has
#' a larger chance of death at early times but less at late times, or similar.
#'
#' I've chosen to use time at which the survival fraction falls below a certain
#' value as a proxy for 'risk'. Given our data, that fraction needs to be quite
#' high (few survival curves make it below 0.5, for example).

#+ rf_c-index

surv.fraction <- 0.9

rf.risk.proxy <-
  # Must be ncol minus the position because the prediction value needs to be a
  # 'risk', ie bigger for earlier death...
  ncol(predict.rf$survival) -
    apply(
      abs(predict.rf$survival - surv.fraction),
      1, which.min
    )

c.index.rf <-
  survConcordance(
    Surv(time_death, surv_event) ~ rf.risk.proxy,
    COHORT.test
  )

print(c.index.rf)

#' The C-index is **`r round(c.index.rf$concordance, 3)`**. This seems to be
#' slightly worse than the conventional Cox model, or indeed the baseline we're
#' trying to beat...
#'
#' Let's try a few other values of ``surv.fraction``, just to be sure.

#+ concordance_table, results='asis'

surv.fractions <- c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5)
c.indices.rf.1 <- c()

for(surv.fraction in surv.fractions) {
  rf.risk.proxy <-
    ncol(predict.rf$survival) -
      apply(
        abs(predict.rf$survival - surv.fraction),
        1, which.min
      )

  c.indices.rf.1 <-
  c(
    c.indices.rf.1,
    survConcordance(
      Surv(time_death, surv_event) ~ rf.risk.proxy,
      COHORT.test
    )$concordance
  )
}

require(xtable)
print(
  xtable(
    data.frame(
      surv.fraction = surv.fractions,
      concordance = c.indices.rf.1
    ),
    digits = c(0,1,4)
  ),
  type = 'html',
  include.rownames = FALSE
)

#' This shows an issue with calculating the concordance this way: it
#' varies depending on the value you pick! This is probably because if you
#' choose too high a value, there's a lot of noise on who gets there first (see
#' the survival curves below...they're kinda jumpy) and, if you choose too high
#' a value, many curves never drop that far and you get a lot of tied risks
#' (which will count against the concordance where patients aren't tied in
#' survival time, which is often.)

#' Let's try something similar, but instead looking at 'risk' as defined by the
#' odds of surviving until a specific time point...

#+ concordance_table_risk_time, results='asis'

# Choose a time to get the risk at
risk.times <- 1:10 # years
c.indices.rf.2 <- c()

for(risk.time in risk.times) {
  # Find the closest bin to that time
  risk.bin <- which.min(abs(predict.rf$unique.death.times - risk.time))
  # Get the chance of having died (ie 1 - survival) for all patients at that time (ie in that bin)
  rf.risk.proxy <- 1 - predict.rf$survival[, risk.bin]
  # Calculate C-index
  c.index <-
    survConcordance(
      Surv(time_death, surv_event) ~ rf.risk.proxy, COHORT.test
    )
  c.indices.rf.2 <- c(c.indices.rf.2, c.index$concordance)
}

require(xtable)
print(
  xtable(
    data.frame(
      risk.time = risk.times,
      concordance = c.indices.rf.2
    ),
    digits = c(0,0,4)
  ),
  type = 'html',
  include.rownames = FALSE
)

#' This suffers from the same issue, but the variability is a little smaller
#' (a range of `r round(diff(range(c.indices.rf.2)), 3)` using this approach,
#' versus `r round(diff(range(c.indices.rf.1)), 3)` using the time to a given
#' survival fraction). Whatever happens though, the values seem to hover around
#' 0.79, which is respectable, but not as good as the Cox models.
#'
#' Finally, plot the same 100 random patients as before...

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