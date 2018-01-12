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

output.filename.base <- '../../output/cox-discrete-elasticnet-08'

bootstraps <- 100
bootstrap.filename <- paste0(output.filename.base, '-boot-all.csv')

n.data <- NA # This is after any variables being excluded in prep

n.threads <- 16

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
#' 
#+ elastic_net_full

require(glmnet)
initParallel(n.threads)

time.start <- handyTimer()

alphas <- seq(0, 1, length.out = 101)
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
  best.lambda.i <- which(cv.fit$lambda == cv.fit$lambda.min)
  mse <- c(mse, cv.fit$cvm[best.lambda.i])
}

time.cv <- handyTimer(time.start)

write.csv(
  data.frame(
    alphas, mse
  ),
  paste0(output.filename.base, '-alpha-calibration.csv')
)

#' `r length(alphas)` alpha values tested in `r time.cv` seconds!

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

# Save for future use
saveRDS(cv.fit, 'cv.fit.rds')

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
      as.formula(paste0('Surv(surv_time, surv_event) ~ risk')),
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
  # Work out risks at risk.time for the special case of a glmnet model
  
  # Derive baseline hazard from cv.glmnet model, heavily based on the
  # glmnet.survcurve and glmnet.basesurv functions in hdnom...
  
  # Get predictions from the training set, because it's the training set whose
  # baseline hazard we need
  
  # This is the relevant section of glmnet.survcurve from 02-hdnom-nomogram.R:
  # lp = as.numeric(predict(object, newx = data.matrix(x),
  #                         s = object$'lambda', type = 'link'))
  # lp means linear predictor from predict.glmnet, because type = 'link'
  lp <-
    as.numeric(
      predict(
        model.fit,
        newx =
          data.matrix(
            dm[-test.set, !(colnames(dm) %in% c('surv_time', 'surv_event'))]
          ),
        s = model.fit$lambda.1se,
        type = 'link'
      )
    )
  # At all unique times in the training set...
  t.unique <-
    # MUST sort these or the cumulative sum below will go crazy!
    sort(unique(dm[-test.set, 'surv_time'][dm[-test.set, 'surv_event'] == 1L]))
  
  alpha <- c()
  for (i in 1:length(t.unique)) {
    # ...loop over calculating the fraction of the population which dies at each
    # timepoint
    alpha[i] <-
      sum(
        # Training set 
        dm[-test.set, 'surv_time'][
          dm[-test.set, 'surv_event'] == 1
        ] == t.unique[i]
      ) /
      sum(
        exp(lp[dm[-test.set, 'surv_time'] >= t.unique[i]])
      )
  }
  
  # Get the cumulative hazard at risk.time by interpolating...
  baseline.cumhaz <-
    approx(
      t.unique, cumsum(alpha), yleft = 0, xout = risk.time, rule = 2
    )$y
  
  # Get predictions from the test set to modify the baseline hazard with
  lp.test <-
    as.numeric(
      predict(
        model.fit,
        newx =
          data.matrix(
            dm[test.set, !(colnames(dm) %in% c('surv_time', 'surv_event'))]
          ),
        s = model.fit$lambda.1se,
        type = 'link'
      )
    )
  
  # 1 minus to get % dead rather than alive
  risks <- 1 - exp(-exp(lp.test) * (baseline.cumhaz))
 
  calibration.table <-
    data.frame(
      surv_event = dm[test.set, 'surv_event'],
      surv_time = dm[test.set, 'surv_time'],
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

calibrationPlot(calibration.table, show.censored = TRUE, max.points = 10000)

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
#' 
#' ## Bootstrapping
#' 
#' Having got those results for a single run on all the data, now bootstrap to
#' find sample-induced variability in performance statistics. Again, because
#' glmnet requires a matrix rather than a data frame this would require a large
#' amount of code in ``handymedical.R``, so do this manually.
#'
#+ bootstrap_performance

time.start <- handyTimer()

# Instantiate a blank data frame
bootstrap.params <- data.frame()

for(i in 1:bootstraps) {
  # Take a bootstrap sample of the training set. We do this with COHORT.prep for
  # the variable importance calculations later.
  COHORT.prep.boot <- bootstrapSampleDf(COHORT.prep[-test.set, ])
  
  # Create a binary matrix for fitting
  COHORT.boot <- convertFactorsToBinaryColumns(COHORT.prep.boot)
  # model.matrix renames logicals to varTRUE, so fix that for status
  colnames(COHORT.boot)[colnames(COHORT.boot) == 'surv_eventTRUE'] <- 'surv_event'
  
  # Fit, but with alpha fixed on the optimal value
  cv.fit.boot <- #readRDS('cv.fit.rds')
    cv.glmnet(
      COHORT.boot[, !(colnames(COHORT.boot) %in% c('surv_time', 'surv_event'))],
      Surv(COHORT.boot[, 'surv_time'], COHORT.boot[, 'surv_event']),
      family = "cox",
      maxit = 1000,
      alpha = alpha.best
    )

  c.index.boot <- glmnetCIndex(cv.fit.boot, COHORT.bin[test.set,])
  
  calibration.table.boot <-
    glmnetCalibrationTable(
      cv.fit.boot, rbind(COHORT.boot, COHORT.bin[test.set, ]),
      test.set = (nrow(COHORT.boot) + 1):(nrow(COHORT.boot) + nrow(COHORT.bin[test.set, ]))
    )
  
  calibration.boot <- calibrationScore(calibration.table.boot)
  
  print(calibrationPlot(calibration.table.boot))
  
  var.imp.vector <- c()
  # Loop over variables to get variable importance
  for(
    var in
    colnames(
      COHORT.prep.boot[, !(colnames(COHORT.prep.boot) %in% c('surv_time', 'surv_event'))]
    )
  ) {
    # Create a dummy data frame and scramble the column var
    COHORT.vimp <- COHORT.prep.boot
    COHORT.vimp[, var] <- sample(COHORT.vimp[, var], replace = TRUE)
    # Make it into a model matrix for fitting
    COHORT.vimp <- convertFactorsToBinaryColumns(COHORT.vimp)
    # model.matrix renames logicals to varTRUE, so fix that for status
    colnames(COHORT.vimp)[colnames(COHORT.vimp) == 'surv_eventTRUE'] <- 'surv_event'
    # Calculate the new C-index
    c.index.vimp <- glmnetCIndex(cv.fit.boot, COHORT.vimp)

    # Append the difference between the C-index with scrambling and the original
    var.imp.vector <-
      c(
        var.imp.vector,
        c.index.boot - c.index.vimp
      )
  }

  names(var.imp.vector) <-
    paste0(
      'vimp.c.index.',
      colnames(
        COHORT.prep.boot[
          ,
          !(colnames(COHORT.prep.boot) %in% c('surv_time', 'surv_event'))
          ]
      )
    )
  
  bootstrap.params <-
    rbind(
      bootstrap.params,
      data.frame(
        t(var.imp.vector),
        c.index = c.index.boot,
        calibration.score = calibration.boot
      )
    )
  
  # Save the bootstrap parameters for later use
  write.csv(bootstrap.params, bootstrap.filename)
}

time.boot.final <- handyTimer(time.start)

#' `r bootstraps` bootstrap fits completed in `r time.boot.final` seconds!

# Get coefficients and variable importances from bootstrap fits
surv.model.fit.coeffs <- bootStatsDf(bootstrap.params)

print(surv.model.fit.coeffs)

# Save performance results
varsToTable(
  data.frame(
    model = 'cox-elnet',
    imputation = FALSE,
    discretised = TRUE,
    c.index = surv.model.fit.coeffs['c.index', 'val'],
    c.index.lower = surv.model.fit.coeffs['c.index', 'lower'],
    c.index.upper = surv.model.fit.coeffs['c.index', 'upper'],
    calibration.score = surv.model.fit.coeffs['calibration.score', 'val'],
    calibration.score.lower =
      surv.model.fit.coeffs['calibration.score', 'lower'],
    calibration.score.upper =
      surv.model.fit.coeffs['calibration.score', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)

#' The bootstrapped C-index is
#' **`r round(surv.model.fit.coeffs['c.index', 'val'], 3)`
#' (`r round(surv.model.fit.coeffs['c.index', 'lower'], 3)` - 
#' `r round(surv.model.fit.coeffs['c.index', 'upper'], 3)`)**
#' on the held-out test set.
#' 
#' The bootstrapped calibration score is
#' **`r round(surv.model.fit.coeffs['calibration.score', 'val'], 3)`
#' (`r round(surv.model.fit.coeffs['calibration.score', 'lower'], 3)` - 
#' `r round(surv.model.fit.coeffs['calibration.score', 'upper'], 3)`)**.
#' 
#' ### Variable importances
#' 
#' Top 20 most important variables from the most recent bootstrap. (This is
#' obviously indicative but just to plot a quick graph and get an idea.)
#' 
#+ bootstrap_var_imp

boot.var.imp.ordered <-
  data.frame(
    var = textAfter(names(var.imp.vector), 'vimp.c.index.'),
    val = var.imp.vector,
    stringsAsFactors = FALSE
  )

boot.var.imp.ordered$desc <- lookUpDescriptions(boot.var.imp.ordered$var)

ggplot(
    boot.var.imp.ordered[order(boot.var.imp.ordered$val[1:20], decreasing = TRUE), ],
    aes(x = var, y = val)
  ) +
  geom_bar(stat = 'identity')