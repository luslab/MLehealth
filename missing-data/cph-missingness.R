#' # A Cox model based purely on missingness
#' 
#' 

data.filename <- '../../data/cohort-sanitised.csv'
calibration.filename <- '../../output/cph-crossvalidation-try1.csv'
results.filename <- '../../output/cph-crossvalidation-results-try1.csv'

n.data <- 30000 # This is of full dataset...further rows may be excluded in prep

poss.missing.vars <-
  c(
    'total_chol_6mo', 'hdl_6mo', 'pulse_6mo', 'crea_6mo',
    'total_wbc_6mo', 'haemoglobin_6mo'
    )

source('../random-forest/shared.R')

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
COHORT.full <- data.frame(fread(data.filename))

# If n.data was specified...
if(!is.na(n.data)){
  # Take a subset n.data in size
  COHORT.use <- sample.df(COHORT.full, n.data)
  rm(COHORT.full)
} else {
  # Use all the data
  COHORT.use <- COHORT.full
  rm(COHORT.full)
}

# Prepare the data
COHORT.prep <-
  prepData(
    COHORT.use,
    cols.keep, discretise.settings, surv.time, surv.event,
    surv.event.yes, extra.fun = caliberExtraPrep, n.keep = n.data
  )
n.data <- nrow(COHORT.prep)

# Define indices of test set
test.set <- sample(1:n.data, (1/3)*n.data)

# Go through columns, dividing them into missing and nonmissing
COHORT.missing <-
  COHORT.prep[, c(poss.missing.vars, 'time_death', 'surv_event')]
for(poss.missing.col in poss.missing.vars) {
  COHORT.missing[, poss.missing.col] <-
    COHORT.prep[, poss.missing.col] == 'missing'
}

# Fit model to the training set
surv.model.fit <-
  survivalFit(
    poss.missing.vars,
    COHORT.missing[-test.set,],
    model.type = 'cph'
  )

# Get C-indices for training and test sets
c.index.train <-
  cIndex(surv.model.fit, COHORT.missing[-test.set, ], model.type = 'cph')
c.index.test <- 
  cIndex(surv.model.fit, COHORT.missing[test.set, ], model.type = 'cph')

#' # Results
#' 
#' ## Performance
#' 
#' The C-index on the full training set is **`r round(c.index.train, 3)`**.
#' 
#' The C-index on the held-out test set is **`r round(c.index.test, 3)`**.
#' 
#' ## Model fit
#' 
#+ resulting_fit

print(surv.model.fit)

#' ## Cox coefficients
#'
#+ cox_coefficients_plot

# TODO: doesn't work for some reason
cph.coeffs <- cphCoeffs(surv.model.fit, COHORT.missing)

ggplot(
  cph.coeffs,
  # TODO: needs variable putting in here too - check name when cph.coeffs works
  aes(x = paste0(level,), y = exp(beta))
) +   
  geom_bar(
    aes(fill = var),
    width=0.9,
    stat = "identity"
  ) +
  geom_text(aes(label = level),angle = 90, hjust = 0)


#' ## Number of missing values as a predictor
#' 

COHORT.missing$n.missing <- rowSums(COHORT.missing[, poss.missing.vars])

surv.model.n.missing <-
  survivalFit(
    'n.missing',
    COHORT.missing[-test.set,],
    model.type = 'cph'
  )

c.index.train <-
  cIndex(surv.model.n.missing, COHORT.missing[-test.set, ], model.type = 'cph')
c.index.test <- 
  cIndex(surv.model.n.missing, COHORT.missing[test.set, ], model.type = 'cph')