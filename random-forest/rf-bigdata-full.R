#' # Summarising big data
#' 
#' Having extracted a huge number of variables, let's find out what we got...
#' 
#+ data_setup

output.filename.base <- '../../output/rf-bigdata-try5'

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

if(!is.na(n.vars)) {
  bigdata.columns <- sort(missingness[bigdata.columns])[1:n.vars]
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

time.start <- handyTimer()
surv.model.fit <-
  survivalFit(
    surv.predict,
    COHORT.bigdata[-test.set,],
    model.type = 'rfsrc',
    n.trees = 16,
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
#' **`r round(calibration.score['area'], 3)`** +/-
#' **`r round(calibration.score['se'], 3)`**.
#' 
#' ## Variable importance
#' 

vimp.df <-
  data.frame(
    var = names(var.imp$importance),
    imp = normalise(var.imp$importance, max)
  )

vimp.df <- lookUpDescriptions(vimp.df)

#' What are the most important variables?

print(
  head(
    vimp.df[order(vimp.df$imp, decreasing = TRUE), c('description', 'imp')],
    20
  )
)

#' ## Partial effects plots
#' 
#' Plot the partial effects of the top 10 most significant variables

for(variable in vimp.df$var[order(vimp.df$imp, decreasing = TRUE)[1:10]]) {
  risk.by.variable <- generalEffectDf(surv.model.fit, COHORT.bigdata, variable)
  
  # Get the mean of the normalised risk for every value of the variable
  risk.aggregated <-
    aggregate(
      as.formula(paste0('risk.normalised ~ ', variable)),
      risk.by.variable, mean
    )
  
  # work out the limits on the x-axis by taking the 1st and 99th percentiles
  x.axis.limits <-
    quantile(COHORT.full[, variable], c(0.01, 0.99), na.rm = TRUE)
  
  print(
    ggplot(risk.by.variable, aes_string(x = variable, y = 'risk.normalised')) +
      geom_line(alpha=0.01, aes(group = id)) +
      geom_line(data = risk.aggregated, colour = 'blue') +
      coord_cartesian(xlim = c(x.axis.limits))
  )
}
