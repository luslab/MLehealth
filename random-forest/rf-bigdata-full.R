#' # Summarising big data
#' 
#' Having extracted a huge number of variables, let's find out what we got...

data.filename.big <- '../../data/cohort-datadriven-02.csv'
model.type <- 'rfsrc'
n.data <- NA

surv.predict.old <- c('age', 'smokstatus', 'imd_score')
untransformed.vars <- c('time_death', 'endpoint_death', 'imd_score', 'exclude')
exclude.vars <-
  c(
    # Entity type 4 is smoking status, which we already have
    "clinical.values.4_data1", "clinical.values.4_data5",
    "clinical.values.4_data6",
    # Entity 13 data2 is the patient's weight centile, and not a single one is
    # entered, but they come out as 0 so the algorithm, looking for NAs, thinks
    # it's a useful column
    "clinical.values.13_data2"
  )

source('../lib/shared.R')

COHORT <- fread(data.filename.big)

summary2 <- function(x) {
  if('data.frame' %in% class(x)) {
    lapply(x, summary2)
  } else {
    if(length(unique(x)) < 30) {
      if(length(unique(x)) < 10) {
        return(round(c(table(x))/length(x), 3)*100)
      } else {
        summ <- sort(table(x), decreasing = TRUE)
        return(
            round(
              c(
                summ[1:5],
                other = sum(summ[6:length(summ)]),
                missing = sum(is.na(x))
              # divide all by the length and turn into %
              )/length(x), 3)*100
        )
      }
    } else {
      return(
        c(
          min = min(x, na.rm = TRUE),
          max = max(x, na.rm = TRUE),
          median = median(x, na.rm = TRUE),
          missing = round(sum(is.na(x))/length(x), 3)*100
        )
      )
    }
  }
}

percentMissing <- function(x) {
  sum(is.na(x))/length(x) * 100
}

startsWithAny <- function(x, prefixes) {
  apply(
    # sapply startsWith over all prefixes, giving a table of logicals
    sapply(
      prefixes,
      function(prefix) {
        startsWith(x, prefix)
      }
    ),
    # ...then, apply a massive logical or over the rows of that table
    MARGIN = 1, FUN = any
  )
}


COHORT.summary <- summary2(COHORT)

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
  which(
    # Does is start with one of the data column names?
    startsWithAny(names(COHORT), bigdata.prefixes) &
      # And it's not one of the columns we want to exclude?
      !(names(COHORT) %in% exclude.vars)
  )

top.bigdata <- sort(missingness[bigdata.columns])[1:100]

COHORT.bigdata <-
  COHORT[, c(
    untransformed.vars, surv.predict.old, names(top.bigdata)
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
# Missing values can therefore be replaced with -1, because, if this is to be
# detected, it will be in the future
for (j in time.based.vars) {
  set(COHORT.bigdata, j = j, value = NA2val(COHORT.bigdata[[j]], -1))
}

# The drug data are number of packs prescribed recently, so missing data can be
# replaced with 0, because none were prescribed
prescriptions.vars <- names(COHORT.bigdata)[startsWith(names(COHORT.bigdata), 'bnf.')]
for (j in prescriptions.vars) {
  set(COHORT.bigdata, j = j, value = NA2val(COHORT.bigdata[[j]], 0))
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

