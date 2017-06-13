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

output.filename.base <- '../../output/rf-bigdata-try6-top100'

data.filename.big <- '../../data/cohort-datadriven-02.csv'
model.type <- 'rfsrc'
n.data <- NA
n.vars <- NA

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
    # Entity 148 is to do with death certification. I'm not sure how it made it
    # into the dataset, but since all the datapoints in this are looking back
    # in time, they're thankfully all NA. This causes rfsrc to fall over.
    "clinical.values.148_data1", "clinical.values.148_data2",
    "clinical.values.148_data3", "clinical.values.148_data4",
    "clinical.values.148_data5",
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

#' # Fitting
#' 
#' Let's fit the model...
#+ rf_bigdata_fit, cache=cacheoption

time.start <- handyTimer()
surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.bigdata[-test.set,],
    model.type = 'rfsrc',
    n.trees = 2000,
    split.rule = split.rule,
    n.threads = 16,
    nimpute = 3,
    nsplit = 10,
    mtry = 500, # There are about 600 variables, so try most of them each time
    na.action = 'na.impute'
  )
time.fit <- handyTimer(time.start)

time.start <- handyTimer()
var.imp <- vimp(surv.model.fit, importance = 'permute.ensemble')
time.vimp <- handyTimer(time.start)

time.start <- handyTimer()
c.index <-
  cIndex(surv.model.fit, COHORT.bigdata[test.set,], na.action = 'na.impute')
time.cindex <- handyTimer(time.start)

#' ## Model fit
#' 
#+ resulting_fit

print(surv.model.fit)

#' ## Timing
#' 
#' The model took `r round(time.fit)` seconds to fit, `r round(time.vimp)`
#' seconds to compute variable importance and `r round(time.cindex)` seconds to
#' calculate the C-index.
#'
#' # Results
#' 
#' ## Performance
#' 
#' ### Discrimination
#' 
#' C-index is **`r round(c.index, 4)`** on the held-out test set.
#' 
#' ### Calibration
#' 
#' Does the model predict realistic probabilities of an event?
#' 
#+ calibration_plot

calibration.table <-
  calibrationTable(
    # Standard calibration options
    surv.model.fit, COHORT.bigdata[test.set,],
    # Always need to specify NA imputation for rfsrc
    na.action = 'na.impute'
  )

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' The area between the calibration curve and the diagonal is 
#' **`r round(calibration.score[['area']], 3)`** +/-
#' **`r round(calibration.score[['se']], 3)`**.
#' 
#' ## Variable importance
#' 

vimp.df <-
  data.frame(
    var = names(var.imp$importance),
    imp = normalise(var.imp$importance, max),
    stringsAsFactors = FALSE
  )


vimp.df <- lookUpDescriptions(vimp.df)

vimp.df$missingness <- missingness[vimp.df$var]

#' What are the most important variables?

print(
  head(
    vimp.df[order(vimp.df$imp, decreasing = TRUE),
            c('imp', 'description', 'missingness')],
    20
  )
)

#' ## Partial effects plots
#' 
#' Plot the partial effects of the top 10 most significant variables

# Add in a column for rounded time to death because that's the y-variable here
COHORT.bigdata$surv_time_round <- round_any(COHORT.bigdata$surv_time, 0.1)

risk.by.logical <- data.frame()

for(variable in vimp.df$var[order(vimp.df$imp, decreasing = TRUE)[1:10]]) {
  risk.by.variable <-
    partialEffectTable(
      surv.model.fit, COHORT.bigdata[-test.set, ], variable,
      na.action = 'na.impute'
    )
  
  if(is.logical(risk.by.variable[, variable])) {
    # If it's a logical value, then make FALSE the baseline and just give values
    # when true, and then append this data frame into a larger one storing all
    # the values from logical variables
    
    # First, get average responses on a per-variable basis
    
    
    risk.by.logical <-
      rbind(
        risk.by.logical,
        data.frame(
          # Store variable name
          variable,
          # And ratio of values associated with TRUEs and FALSEs
          tf.ratio =
            # Every 2nd one is a TRUE, because the values are sorted 
            risk.by.variable[
              seq.int(2, nrow(risk.by.variable), 2), 'risk.normalised'
              ] /
            # Conversely, all the odd rows are FALSE
            risk.by.variable[
              seq.int(1,nrow(risk.by.variable) - 1, 2), 'risk.normalised'
              ]
        )
      )
  }
  
  # Get the mean of the normalised risk for every value of the variable
  risk.aggregated <-
    aggregate(
      as.formula(paste0('risk.normalised ~ ', variable)),
      risk.by.variable, median
    )
  
  if(is.factor(risk.by.variable[, variable])) {
    # If it's a factor, draw a bar chart for the per-level partial effects
    print(
      ggplot(risk.by.variable, aes_string(x = variable, y = 'risk.normalised')) +
        geom_point(alpha = 0.01, position = position_jitter(w = 0.8, h = 0)) +
        geom_bar(data = risk.aggregated, colour = 'blue')
    )
  } else if(is.numeric(risk.by.variable[, variable])) {
    # If it's a number, draw a line plot across the range of values
    
    
    # work out the limits on the x-axis by taking the 1st and 99th percentiles
    x.axis.limits <-
      quantile(COHORT.bigdata[, variable], c(0.01, 0.99), na.rm = TRUE)
    print(
      ggplot(risk.by.variable, aes_string(x = variable, y = 'risk.normalised')) +
        geom_line(alpha=0.01, aes(group = id)) +
        geom_line(data = risk.aggregated, colour = 'blue') +
        coord_cartesian(xlim = c(x.axis.limits))
    )
  }
}

if(nrow(risk.by.logical) > 0) {
  # If there were some logical variables in there, let's plot them all on one
  # bar plot, much like the factor plot from before
  
  # Get the averaged risk for each variable
  risk.aggregated <- aggregate(tf.ratio ~ variable, risk.by.logical, median)
  
  print(
    ggplot(risk.by.logical, aes(x = variable, y = tf.ratio)) +
      geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
      coord_cartesian(ylim = c(0, 2))
  )
}