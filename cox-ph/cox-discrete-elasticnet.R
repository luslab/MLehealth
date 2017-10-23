#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- TRUE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Variable selection in data-driven health records with discretised
#' # Cox models
#' 
#' Having extracted around 600 variables which occur most frequently in patient
#' records, let's try to narrow these down using a methodology based on varSelRf
#' combined with survival modelling. We'll find the predictability of variables
#' as defined by the p-value of a logrank test on survival curves of different
#' categories within that variable, and then iteratively throw out unimportant
#' variables, cross-validating for optimum performance.
#' 
#' ## User variables
#' 
#+ user_variables

output.filename.base <- '../../output/cox-bigdata-varsellogrank-01'

cv.n.folds <- 3
vars.drop.frac <- 0.2 # Fraction of variables to drop at each iteration
bootstraps <- 100

n.data <- NA # This is after any variables being excluded in prep

n.threads <- 20

#' ## Data set-up
#' 
#+ data_setup

data.filename.big <- '../../data/cohort-datadriven-02.csv'

surv.predict.old <- c('age', 'smokstatus', 'imd_score', 'gender')
untransformed.vars <- c('time_death', 'endpoint_death', 'exclude')

source('../lib/shared.R')
require(xtable)

# Define these after shared.R or they will be overwritten!
exclude.vars <-
  c(
    # Entity type 4 is smoking status, which we already have
    "clinical.values.4_data1", "clinical.values.4_data5",
    "clinical.values.4_data6",
    # Entity 13 data2 is the patient's weight centile, and not a single one is
    # entered, but they come out as 0 so the algorithm, looking for NAs, thinks
    # it's a useful column
    "clinical.values.13_data2",
    # Entities 148 and 149 are to do with death certification. I'm not sure how 
    # it made it into the dataset, but since all the datapoints in this are
    # looking back in time, they're all NA. This causes rfsrc to fail.
    "clinical.values.148_data1", "clinical.values.148_data2",
    "clinical.values.148_data3", "clinical.values.148_data4",
    "clinical.values.148_data5",
    "clinical.values.149_data1", "clinical.values.149_data2",
    # These are all the same value except where NA, which causes issues with
    # discretisation
    "clinical.values.14_data2", "clinical.values.62_data1",
    "clinical.values.64_data1", "clinical.values.65_data1",
    "clinical.values.65_data1", "clinical.values.67_data1",
    "clinical.values.68_data2", "clinical.values.70_data1"
  )

COHORT <- fread(data.filename.big)

bigdata.prefixes <-
  c(
    'hes.icd.',
    'hes.opcs.',
    'tests.enttype.',
    'clinical.history.',
    'clinical.values.',
    'bnf.'
  )

bigdata.columns <-
  colnames(COHORT)[
    which(
      # Does is start with one of the data column names?
      startsWithAny(names(COHORT), bigdata.prefixes) &
        # And it's not one of the columns we want to exclude?
        !(colnames(COHORT) %in% exclude.vars)
    )
    ]

COHORT.bigdata <-
  COHORT[, c(
    untransformed.vars, surv.predict.old, bigdata.columns
  ),
  with = FALSE
  ]

# Get the missingness before we start removing missing values
missingness <- sort(sapply(COHORT.bigdata, percentMissing))
# Remove values for the 'untransformed.vars' above, which are the survival
# values plus exclude column
missingness <- missingness[!(names(missingness) %in% untransformed.vars)]

# Deal appropriately with missing data
# Most of the variables are number of days since the first record of that type
time.based.vars <-
  names(COHORT.bigdata)[
    startsWithAny(
      names(COHORT.bigdata),
      c('hes.icd.', 'hes.opcs.', 'clinical.history.')
    )
    ]
# We're dealing with this as a logical, so we want non-NA values to be TRUE,
# is there is something in the history
for (j in time.based.vars) {
  set(COHORT.bigdata, j = j, value = !is.na(COHORT.bigdata[[j]]))
}

# Again, taking this as a logical, set any non-NA value to TRUE.
prescriptions.vars <- names(COHORT.bigdata)[startsWith(names(COHORT.bigdata), 'bnf.')]
for (j in prescriptions.vars) {
  set(COHORT.bigdata, j = j, value = !is.na(COHORT.bigdata[[j]]))
}

# This leaves tests and clinical.values, which are test results and should be
# imputed.

# Manually fix clinical values items...
#
# "clinical.values.1_data1"  "clinical.values.1_data2"
# These are just blood pressure values...fine to impute
#
# "clinical.values.13_data1" "clinical.values.13_data3"
# These are weight and BMI...also fine to impute
#
# Entity 5 is alcohol consumption status, 1 = Yes, 2 = No, 3 = Ex, so should be
# a factor, and NA can be a factor level
COHORT.bigdata$clinical.values.5_data1 <-
  factorNAfix(factor(COHORT.bigdata$clinical.values.5_data1), NAval = 'missing')

# Both gender and smokstatus are factors...fix that
COHORT.bigdata$gender <- factor(COHORT.bigdata$gender)
COHORT.bigdata$smokstatus <-
  factorNAfix(factor(COHORT.bigdata$smokstatus), NAval = 'missing')

# Exclude invalid patients
COHORT.bigdata <- COHORT.bigdata[!COHORT.bigdata$exclude]
COHORT.bigdata$exclude <- NULL

# Remove negative survival times
COHORT.bigdata <- subset(COHORT.bigdata, time_death > 0)

# Define test set
test.set <- testSetIndices(COHORT.bigdata, random.seed = 78361)

# If n.data was specified, trim the data table down to size
if(!is.na(n.data)) {
  COHORT.bigdata <- sample.df(COHORT.bigdata, n.data)
}

# Create an appropraite survival column
COHORT.bigdata <- 
  prepSurvCol(
    data.frame(COHORT.bigdata), 'time_death', 'endpoint_death', 'Death'
  )

# Start by predicting survival with all the variables provided
surv.predict <- c(surv.predict.old, bigdata.columns)

# Set up a csv file to store calibration data, or retrieve previous data
calibration.filename <- paste0(output.filename.base, '-varselcalibration.csv')

# Create process settings

# Variables to leave alone, including those whose logrank p-value is NA because
# that means there is only one value in the column and so it can't be discretised
# properly anyway
vars.noprocess <- c('surv_time', 'surv_event')
process.settings <-
  list(
    var        = vars.noprocess,
    method     = rep(NA, length(vars.noprocess)),
    settings   = rep(list(NA), length(vars.noprocess))
  )
# Find continuous variables which will need discretising
continuous.vars <- names(COHORT.bigdata)[sapply(COHORT.bigdata, class) %in% c('integer', 'numeric')]
# Remove those variables already explicitly excluded, mainly for those whose
# logrank score was NA
continuous.vars <- continuous.vars[!(continuous.vars %in% process.settings$var)]
process.settings$var <- c(process.settings$var, continuous.vars)
process.settings$method <-
  c(process.settings$method,
    rep('binByQuantile', length(continuous.vars))
  )
process.settings$settings <-
  c(
    process.settings$settings,
    rep(
      list(
        seq(
          # Quantiles are obviously between 0 and 1
          0, 1,
          # All have the same number of bins
          length.out = 10
        )
      ),
      length(continuous.vars)
    )
  )


# Need a way to ID in advance those which are going to fail here, ie those where
# there are no quantiles. The 

COHORT.prep <-
  prepData(
    # Data for cross-validation excludes test set
    COHORT.bigdata,
    names(COHORT.bigdata),
    process.settings,
    'surv_time', 'surv_event',
    TRUE
  )

# Kludge...remove surv_time.1 and rename surv_event.1
COHORT.prep$surv_time.1 <- NULL
names(COHORT.prep)[names(COHORT.prep) == 'surv_event.1'] <- 'surv_event'

COHORT.bin <- convertFactorsToBinaryColumns(COHORT.prep)
# model.matrix renames logicals to varTRUE, so fix that for status
colnames(COHORT.bin)[colnames(COHORT.bin) == 'surv_eventTRUE'] <- 'surv_event'

test.set <- testSetIndices(COHORT.bin)

# Coxnet code, should you ever decide to go that route
# test <-
#   Coxnet(
#     data.matrix(COHORT.bin[-test.set, !(colnames(COHORT.bin) %in% c('time', 'status'))]),
#     data.matrix(COHORT.bin[-test.set, c('time', 'status')]),
#     penalty = 'Enet',
#     alpha = 0,
#     nlambda = 50, nfolds = 10, maxit = 1e+5
#   )

#' ## Elastic net regression
#' 
#' Run a loop over alphas running from LASSO to ridge regression, and see which
#' is best after tenfold cross-validation...

require(glmnet)
initParallel(8)

alphas <- seq(0, 1, length.out = 2)
mse <- c()

for(alpha in alphas) {
  cv.fit <-
    cv.glmnet(
      COHORT.bin[-test.set, !(colnames(COHORT.bin) %in% c('surv_time', 'surv_event'))],
      Surv(COHORT.bin[-test.set, 'surv_time'], COHORT.bin[-test.set, 'surv_event']),
      family = "cox",
      maxit = 1000,
      alpha = alpha,
      parallel = TRUE
    )
  best.lambda.i <- which(cv.fit$lambda == cv.fit$lambda.min) # should this be lambda.1se?
  mse <- c(mse, cv.fit$cvm[best.lambda.i])
}

alpha.best <- alphas[which.min(mse)]

# To avoid saving all the fits, let's just refit the best one
cv.fit <-
  cv.glmnet(
    COHORT.bin[-test.set, !(colnames(COHORT.bin) %in% c('surv_time', 'surv_event'))],
    Surv(COHORT.bin[-test.set, 'surv_time'], COHORT.bin[-test.set, 'surv_event']),
    family = "cox",
    maxit = 1000,
    alpha = alpha.best,
    parallel = TRUE
  )

#' The best alpha was `r alpha.best`, and the lambda with the lowest mean-square
#' error was `r cv.fit$lambda.min`. We'll be using the strictest lambda which is
#' within 1 se of the minimum, `r cv.fit$lambda.1se`.
#'
#' ## Performance
#' 
#' ### C-index
#' 
#' Calculate C-index manually. The glmnet interface requiring matrices is
#' sufficiently different to the usual one that I've not spent time integrating
#' it with the rest of the ``handymedical.R`` functions yet.
#' 
#+ c_index

glmnetCIndex <- function(model.fit, dm) {
  test.predictions <-
    getRisk(
      model.fit,
      dm[, !(colnames(dm) %in% c('surv_time', 'surv_event'))]
    )
  
  as.numeric(
    survConcordance(
      as.formula(paste0('Surv(surv_time, surv_event) ~ surv_event')),
      data.frame(
        surv_time = dm[, 'surv_time'],
        surv_event = dm[, 'surv_event'],
        risk = test.predictions
      )
    )$concordance
  )
}

c.index <- glmnetCIndex(cv.fit, COHORT.bin[test.set, ])

#' C-index is `r c.index`.
#'
#' ### Calibration
#' 
#' For now, calibration is manual too just to get it working. It's surprisingly
#' hard... A package called `c060` should contain a function `predictProb`, but
#' on loading it, is says function not found. So here is a very manual solution,
#' creating a dummy Cox model using the `survival` package, inspired by
#' [this](https://stat.ethz.ch/pipermail/r-help/2012-May/312029.html).
#' 
#+ calibration_plot

glmnetCalibrationTable <- function(model.fit, dm, test.set, risk.time = 5) {
  # Select the coefficients of the model which are greater than zero
  coef.non0 <- as.vector(abs(coef(model.fit, s = "lambda.1se")) > 0)
  # Values of coefficients
  selected.coef <- coef(model.fit, s = "lambda.1se")[coef.non0]
  # Names of coefficients
  selected.vars <-
    colnames(dm)[
      !(colnames(dm) %in% c('surv_time', 'surv_event'))
      ][coef.non0]
  
  # Make a dummy data frame from the model matrix for a dummy Cox model because
  # coxph needs a data frame not a matrix
  dummy.df <- data.frame(COHORT.bin[, c('surv_time', 'surv_event', selected.vars)])
  
  # Round the times in the dummy data frame to save memory and time
  dummy.df$surv_time <- round(dummy.df$surv_time, 1)
  
  dummy.cph <-
    coxph(
      as.formula(
        paste0(
          'Surv(surv_time, surv_event) ~ ',
          # use make.names because turning the binarised matrix into a data frame
          # converts its colnames
          paste0(make.names(selected.vars), collapse = '+')
        )),
      data = dummy.df[-test.set, ],
      init = selected.coef, iter=0
    )
  
  # Perform a fit to get survival curves for the test set
  sfit <- survfit(dummy.cph, newdata = dummy.df[test.set, ])
  
  risktime.col <- which.min(abs(sfit$time - risk.time))
  
  # Returns % survived, we want % dead
  risks <- 1 - as.vector(sfit$surv[risktime.col, ])
  
  calibration.table <-
    data.frame(
      surv_event = dummy.df$surv_event[test.set],
      surv_time = dummy.df$surv_time[test.set],
      risk = risks
    )
  
  # Was there an event? Start with NA, because default is unknown (ie censored)
  calibration.table$event <- NA
  # Event before risk.time
  calibration.table$event[
    calibration.table$surv_event & calibration.table$surv_time <= risk.time
    ] <- TRUE
  # Event after, whether censorship or not, means no event by risk.time
  calibration.table$event[calibration.table$surv_time > risk.time] <- FALSE
  # Otherwise, censored before risk.time, leave as NA
  
  # Drop unnecessary columns and return
  calibration.table[, c('risk', 'event')]
}

calibration.table <- glmnetCalibrationTable(cv.fit, COHORT.bin, test.set)

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' Calibration score is `r calibration.score`.

#' ### Coefficients
#' 
#' The elastic net regression generates coefficients for every factor/level
#' combination (for continuous values, this means that every decile gets its own
#' coefficient). The lambda penalty means that some amount of variable selection
#' is performed in this process, penalising too many large coefficients.
#' Depending on the value of alpha, quite a few of the coefficients can be
#' exactly zero. Let's have a look at what we got... 
#' 
#+ coefficients

# Make a data frame of the coefficients in decreasing order
cv.fit.coefficients.ordered <-
  data.frame(
    factorlevel = # Names are ordered in decreasing order of absolute value
      factorOrderedLevels(
        colnames(COHORT.bin)[
          order(abs(coef(cv.fit, s = "lambda.1se")), decreasing = TRUE)
        ]
      ),
    val =
      coef(cv.fit, s = "lambda.1se")[
        order(abs(coef(cv.fit, s = "lambda.1se")), decreasing = TRUE)
      ]
  )

# Get the variable names by removing TRUE, (x,y] or missing from the end
cv.fit.coefficients.ordered$var <-
  gsub('TRUE', '', cv.fit.coefficients.ordered$factorlevel)
cv.fit.coefficients.ordered$var <-
  # Can contain numbers, decimals, e+/- notation and commas separating bounds
  gsub('\\([0-9,.e\\+-]+\\]', '', cv.fit.coefficients.ordered$var)
cv.fit.coefficients.ordered$var <-
  gsub('missing', '', cv.fit.coefficients.ordered$var)
# Kludgey manual fix for 5_data1 which can take values 1, 2 or 3 and is
# therefore very hard to catch
cv.fit.coefficients.ordered$var <-
  gsub(
    'clinical.values.5_data1[0-9]', 'clinical.values.5_data1',
    cv.fit.coefficients.ordered$var
  )

# And then get human-readable descriptions
cv.fit.coefficients.ordered$desc <-
  lookUpDescriptions(cv.fit.coefficients.ordered$var)

#' #### Top 30 coefficients
#'
#+ coefficients_table, results='asis'
print(
  xtable(
    cv.fit.coefficients.ordered[1:30, c('desc', 'factorlevel', 'val')],
    digits = c(0, 0, 0, 3)
  ),
  type = 'html',
  include.rownames = FALSE
)

#' #### Graph of all coefficient values
#' 
#' Nonzero values are red, zero values are blue.

ggplot(cv.fit.coefficients.ordered, aes(x = factorlevel, y = val, colour = val == 0)) +
  geom_point() +
  theme(
    axis.title.x=element_blank(), axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  )

#' Overall, there are `r sum(cv.fit.coefficients.ordered$val != 0)` nonzero
#' coefficients out of `r nrow(cv.fit.coefficients.ordered)`. In the case of
#' multilevel factors or continuous values, multiple coefficients may result
#' from a single variable in the original data. Correcting for this, there are
#' `r length(unique(cv.fit.coefficients.ordered$var[cv.fit.coefficients.ordered$val != 0]))`
#' unique variables represented out of
#' `r length(unique(cv.fit.coefficients.ordered$var))` total variables.
