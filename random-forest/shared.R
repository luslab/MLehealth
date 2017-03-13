#' # prep-data.R
#' We start by preparing the data for reproducible comparisons...

# Set the random seed for reproducibility
random.seed <- 35498L

# Specify the data file containing the patient cohort
cohort.file <- '../../../data/cohort-sanitised.csv'

# The fraction of the data to use as the test set (1 - this will be used as the
# training set)
test.fraction <- 1/3

# If surv.predict wasn't already specified, use the defaults...
if(!exists('surv.predict')) {
  # Column names of variables to use for predictions
  surv.predict <- c(
    'age', 'gender', 'most_deprived', 'diagnosis', 'pci_6mo', 'cabg_6mo',
    'hx_mi', 'long_nitrate', 'smokstatus', 'hypertension', 'diabetes',
    'total_chol_6mo', 'hdl_6mo', 'heart_failure', 'pad', 'hx_af', 'hx_stroke',
    'hx_renal', 'hx_copd', 'hx_cancer', 'hx_liver', 'hx_depression',
    'hx_anxiety', 'pulse_6mo', 'crea_6mo', 'total_wbc_6mo','haemoglobin_6mo',
    'most_deprived'
  )
}

cols.keep <- c(surv.predict, 'exclude', 'imd_score')

exclude.vars <- c('hx_mi')
surv.predict <- surv.predict[!(surv.predict %in% exclude.vars)]

# Check to see if endpoint exists to avoid error
if(!exists('endpoint')) {
  # Default is all-cause mortality
  endpoint <- 'death'
}

# If we're looking at MI...
if(endpoint == 'mi'){
  surv.time      <- 'time_coronary'
  surv.event     <- 'endpoint_coronary'
  surv.event.yes <- c('Nonfatal MI', 'Coronary death')
# Default is all-cause mortality...
} else {
  # Name of column giving time for use in survival object
  surv.time    <- 'time_death'
  # Name of event column for survival object
  surv.event   <- 'endpoint_death' # Cannot be 'surv_event' or will break later!
  # Value of surv.event column if an event is recorded
  surv.event.yes <- 'Death'
}



# Quantile boundaries for discretisation
discretise.quantiles <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)

# Columns to discretise in specific ways, or not discretise at all. Those not
# listed here will be discretised by quantile with the default quantiles listed
# above.
discretise.settings <-
  list(
    var        = c('anonpatid', 'surv_time', 'imd_score', 'exclude'),
    method     = c(NA, NA, NA, NA),
    settings   = list(NA, NA, NA, NA)
  )

################################################################################
### END USER VARIABLES #########################################################
################################################################################

if(!is.na(random.seed)) {
  set.seed(random.seed)
}

source('../lib/handymedical.R', chdir = TRUE)

# Define a function of extra non-general prep to be done on this dataset
caliberExtraPrep <- function(df) {
  df <-
    df[
        # remove negative times to death
        df$surv_time > 0 &
        # remove patients who should be excluded
        !df$exclude
      ,
      ]
  # Remove the exclude column, which we don't need any more
  df$exclude <- NULL

  # Create most_deprived, as defined in the paper: the bottom 20%
  df$most_deprived <- df$imd_score > quantile(df$imd_score, 0.8, na.rm = TRUE)
  df$most_deprived <- factorNAfix(factor(df$most_deprived), NAval = 'missing')
  # Remove the imd_score, to avoid confusion later
  df$imd_score <- NULL
  
  df
}