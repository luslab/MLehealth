#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- TRUE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Summarising big data
#' 
#' Having extracted a huge number of variables, let's find out what we got...
#' 
#' ## Data set-up
#' 
#+ data_setup

output.filename.base <- '../../output/rf-bigdata-try7-ALL-cv-smallmtry'

data.filename.big <- '../../data/cohort-datadriven-02.csv'
model.type <- 'rfsrc'
n.data <- NA
n.vars <- NA
nsplit <- 20
mtry <- 500
n.trees <- 500

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


surv.predict.old <- c('age', 'smokstatus', 'imd_score', 'gender')
untransformed.vars <- c('time_death', 'endpoint_death', 'exclude')

source('../lib/shared.R')

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
    # Entities 148 and 149 are to do with death certification. I'm not sure how it made it
    # into the dataset, but since all the datapoints in this are looking back
    # in time, they're thankfully all NA. This causes rfsrc to fall over.
    "clinical.values.148_data1", "clinical.values.148_data2",
    "clinical.values.148_data3", "clinical.values.148_data4",
    "clinical.values.148_data5",
    "clinical.values.149_data1", "clinical.values.149_data2",
    # From the old exclude.vars
    "hx_mi"
  )

COHORT <- fread(data.filename.big)

percentMissing <- function(x) {
  sum(is.na(x))/length(x) * 100
}

missingness <- sapply(COHORT, percentMissing)

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

# If n.vars was not set as NA, then thin out the variables to that number
if(!is.na(n.vars)) {
  bigdata.columns <-
    bigdata.columns[order(missingness[bigdata.columns])[1:n.vars]]
}

COHORT.bigdata <-
  COHORT[, c(
    untransformed.vars, surv.predict.old, bigdata.columns
  ),
  with = FALSE
  ]

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

COHORT.bigdata <-
  prepSurvCol(data.frame(COHORT.bigdata), 'time_death', 'endpoint_death', 'Death')

test.set <- testSetIndices(COHORT.bigdata, random.seed = 78361)

surv.predict <- c(surv.predict.old, bigdata.columns)



# Copied from rfsrc-cv-mtry-nsplit-logical.R

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
            calibration.score,
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
  
  
  
}