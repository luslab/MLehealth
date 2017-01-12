#' # prep-data.R
#' We start by preparing the data for reproducible comparisons...

# Set the random seed for reproducibility
set.seed(35498)

# Specify the data file containing the patient cohort
cohort.file <- '../../data/cohort-sanitised.csv'

# The fraction of the data to use as the test set (1 - this will be used as the
# training set)
test.fraction <- 1/3

# Column names of variables to use for predictions
surv.predict <- c(
  'age', 'gender', 'most_deprived', 'diagnosis', 'pci_6mo', 'cabg_6mo', 'hx_mi',
  'long_nitrate', 'smokstatus', 'hypertension', 'diabetes', 'total_chol_6mo',
  'hdl_6mo', 'heart_failure', 'pad', 'hx_af', 'hx_stroke', 'hx_renal',
  'hx_copd', 'hx_cancer', 'hx_liver', 'hx_depression', 'hx_anxiety',
  'pulse_6mo', 'crea_6mo', 'total_wbc_6mo','haemoglobin_6mo'
  )

# Name of column giving time for use in survival object
surv.time    <- 'time_death'
# Name of event column for survival object
surv.event   <- 'endpoint_death' # Cannot be 'surv_event' or will break later!
# Value of surv.event column if an event is recorded
surv.event.yes <- 'Death'

# Columns which aren't to be used for prediction, but we want to keep track of
other.cols <- c('anonpatid', 'imd_score')

# Quantile boundaries for discretisation
discretise.quantiles <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1) #c(0, 0.5, 1)

# Columns to discretise in specific ways, or not discretise at all. Those not
# listed here will be discretised by quantile with the default quantiles listed
# above.
discretise.settings <-
  list(
    var        = c('anonpatid','time_death','imd_score'),
    method     = c(NA, NA, NA),
    settings   = list(NA, NA, NA),
    output.var = c(NA, NA, NA)
  )

################################################################################
### END USER VARIABLES #########################################################
################################################################################

source('../lib/handymedical.R', chdir = TRUE)

# make a vector of all columns we'll be using
all.cols <- unique(c(surv.time, surv.event, surv.predict,
                     other.cols))

# load the standard COHORT variable
COHORT <- fread(cohort.file)
# convert it from data table to data frame
COHORT.df <- data.frame(COHORT)

for(cl in all.cols){
  col <- COHORT.df[, which(colnames(COHORT.df)==cl)]
  #print(paste('1', cl,typeof(col), is.character(col), is.factor(col)))
  if(is.character(col)) {
    #print(paste('is char', cl))
    COHORT.df[, which(colnames(COHORT.df)==cl)] <- factor(col)
    col <- COHORT.df[, which(colnames(COHORT.df)==cl)]
  }
  #print(paste('2', cl, typeof(col), is.character(col), is.factor(col)))
  if(is.factor(col)){
    #print(paste('is factor', cl))
    if(anyNA(col)){
      COHORT.df[, which(colnames(COHORT.df)==cl)] <-
        factorNAfix(col, NAval = 'missing')
    }
  }
}
COHORT.to.use <-
  COHORT.df[
            # remove negative times to death
            COHORT.df$time_death > 0,
            # only include the columns we actually need...and some of all.cols
            # won't exist yet, so only include ones which are actually in the
            # data frame or it will throw an error
            all.cols[all.cols %in% names(COHORT.df)]
            ]

# discretise the data by going through column by column and acting appropriately
for(column in names(COHORT.to.use)) {
  # if we have a specific way to discretise this column, let's do it!
  if(column %in% discretise.settings$var) {
    j <- match(column, discretise.settings$var)
    # discretisation method being NA means don't discretise
    if(!is.na(discretise.settings$method[j])) {
      discretise.fun <- match.fun(discretise.settings$method[j])
      # create a column called output.var and populate it
      COHORT.to.use[, discretise.settings$output.var[j]] <-
        discretise.fun(COHORT.to.use[,column],
                       discretise.settings$settings[[j]]
        )
    }
    # if there's no specific method mentioned, but it's numerical, discretise in
    # the default way
  } else if(class(COHORT.to.use[,column]) %in% c('numeric', 'integer')) {
    COHORT.to.use[,column] <-
      binByQuantile(COHORT.to.use[,column], discretise.quantiles)
    # if it's logical, turn it into a two-level fator
  } else if(class(COHORT.to.use[,column]) == 'logical') {
    COHORT.to.use[,column] <- factor(COHORT.to.use[,column])
  }
}

# Create most_deprived, as defined in the paper
COHORT.to.use$most_deprived <-
  COHORT.to.use$imd_score > quantile(COHORT.to.use$imd_score, 0.8, na.rm = TRUE)
COHORT.to.use$most_deprived <-
  factorNAfix(
    factor(COHORT.to.use$most_deprived),
    NAval = 'missing'
  )
# Remove the imd_score, to avoid confusion later
COHORT.to.use$imd_score <- NULL

# Create a column denoting censorship or otherwise of events
COHORT.to.use$surv_event <- COHORT.to.use[, surv.event] == surv.event.yes

# Remove the surv.event column so we don't use it as a covariate later
COHORT.to.use <- COHORT.to.use[, names(COHORT.to.use) != surv.event]

test.rows <- sample(1:nrow(COHORT.to.use), test.fraction * nrow(COHORT.to.use))

COHORT.training.full <- COHORT.to.use[-test.rows,]
COHORT.test.full     <- COHORT.to.use[test.rows,]