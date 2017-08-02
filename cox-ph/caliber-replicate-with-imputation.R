#+ knitr_setup, include = FALSE

# Whether to cache the intensive code sections. Set to FALSE to recalculate
# everything afresh.
cacheoption <- FALSE
# Disable lazy caching globally, because it fails for large objects, and all the
# objects we wish to cache are large...
opts_chunk$set(cache.lazy = FALSE)

#' # Replicating Rapsomaniki _et al._ 2014
#' 
#' ## User variables
#' 
#' First, define some variables...

#+ define_vars

imputed.data.filename <- '../../data/COHORT_complete.rds'
n.data <- NA # This is of full dataset...further rows may be excluded in prep
endpoint <- 'death.imputed'

old.coefficients.filename <- 'rapsomaniki-cox-values-from-paper.csv'

output.filename.base <- '../../output/caliber-replicate-imputed-survreg-3'


cox.var.imp.perm.filename <-
  '../../output/caliber-replicate-imputed-survreg-bootstrap-var-imp-perm-1.csv'
cox.var.imp.perm.missing.filename <-
  '../../output/caliber-replicate-with-missing-survreg-bootstrap-var-imp-perm-1.csv'


bootstraps <- 100
n.threads <- 10

#' ## Setup

#+ setup, message=FALSE

source('../lib/shared.R')
require(xtable)
require(ggrepel)

# Load the data and convert to data frame to make column-selecting code in
# prepData simpler
imputed.data <- readRDS(imputed.data.filename)

# Remove rows with death time of 0 to avoid fitting errors
for(i in 1:length(imputed.data)) {
  imputed.data[[i]] <- imputed.data[[i]][imputed.data[[i]][, surv.time] > 0, ]
}

# Define n.data based on the imputed data, which has already been preprocessed
n.data <- nrow(imputed.data[[1]])
# Define indices of test set
test.set <- testSetIndices(imputed.data[[1]], random.seed = 78361)

#' OK, we've now got **`r n.data`** patients, split into a training set of
#' `r n.data - length(test.set)` and a test set of `r length(test.set)`.
#'
#'
#' ## Transform variables
#' 
#' The model uses variables which have been standardised in various ways, so
#' let's go through and transform our input variables in the same way...

source('caliber-scale.R')

for(i in 1:length(imputed.data)) {
  imputed.data[[i]] <- caliberScale(imputed.data[[i]], surv.time, surv.event)
}

#' ## Survival fitting
#' 
#' Fit a Cox model to the preprocessed data. The paper uses a Cox model with an
#' exponential baseline hazard, as here. The standard errors were calculated
#' with 200 bootstrap samples, which we're also doing here.

#+ fit_cox_model, cache=cacheoption

surv.formula <-
  Surv(surv_time, surv_event) ~
    ### Sociodemographic characteristics #######################################
    ## Age in men, per year
    ## Age in women, per year
    ## Women vs. men
    # ie include interaction between age and gender!
    age*gender +
    ## Most deprived quintile, yes vs. no
    most_deprived +
    ### SCAD diagnosis and severity ############################################
    ## Other CHD vs. stable angina
    ## Unstable angina vs. stable angina
    ## NSTEMI vs. stable angina
    ## STEMI vs. stable angina
    diagnosis +
    #diagnosis_missing +
    ## PCI in last 6 months, yes vs. no
    pci_6mo +
    ## CABG in last 6 months, yes vs. no
    cabg_6mo +
    ## Previous/recurrent MI, yes vs. no
    hx_mi +
    ## Use of nitrates, yes vs. no
    long_nitrate +
    ### CVD risk factors #######################################################
    ## Ex-smoker / current smoker / missing data vs. never
    smokstatus +
    ## Hypertension, present vs. absent
    hypertension +
    ## Diabetes mellitus, present vs. absent
    diabetes_logical +
    ## Total cholesterol, per 1 mmol/L increase
    total_chol_6mo +
    ## HDL, per 0.5 mmol/L increase
    hdl_6mo +
    ### CVD co-morbidities #####################################################
    ## Heart failure, present vs. absent
    heart_failure +
    ## Peripheral arterial disease, present vs. absent
    pad +
    ## Atrial fibrillation, present vs. absent
    hx_af +
    ## Stroke, present vs. absent
    hx_stroke +
    ### Non-CVD comorbidities ##################################################
    ## Chronic kidney disease, present vs. absent
    hx_renal +
    ## Chronic obstructive pulmonary disease, present vs. absent
    hx_copd +
    ## Cancer, present vs. absent
    hx_cancer +
    ## Chronic liver disease, present vs. absent
    hx_liver +
    ### Psychosocial characteristics ###########################################
    ## Depression at diagnosis, present vs. absent
    hx_depression +
    ## Anxiety at diagnosis, present vs. absent
    hx_anxiety +
    ### Biomarkers #############################################################
    ## Heart rate, per 10 b.p.m increase
    pulse_6mo +
    ## Creatinine, per 30 μmol/L increase
    crea_6mo +
    ## White cell count, per 1.5 109/L increase
    total_wbc_6mo +
    ## Haemoglobin, per 1.5 g/dL increase
    haemoglobin_6mo

# Do a quick and dirty fit on a single imputed dataset, to draw calibration
# curve from
fit.exp <- survreg(
  formula = surv.formula,
  data = imputed.data[[1]][-test.set, ],
  dist = "exponential"
)

fit.exp.boot <- list()

# Perform bootstrap fitting for every multiply imputed dataset
for(i in 1:length(imputed.data)) {
  fit.exp.boot[[i]] <-
    boot(
      formula = surv.formula,
      data = imputed.data[[i]][-test.set, ],
      statistic = bootstrapFitSurvreg,
      R = bootstraps,
      parallel = 'multicore',
      ncpus = n.threads,
      test.data = imputed.data[[i]][test.set, ]
    )
}

# Save the fits, because it might've taken a while!
saveRDS(fit.exp.boot, paste0(output.filename.base, '-surv-boot-imp.rds'))

# Unpackage the uncertainties from the bootstrapped data
fit.exp.boot.ests <-  bootMIStats(fit.exp.boot)

# Save bootstrapped performance values
varsToTable(
  data.frame(
    model = 'cox',
    imputation = TRUE,
    discretised = FALSE,
    c.index = fit.exp.boot.ests['c.test', 'val'],
    c.index.lower = fit.exp.boot.ests['c.test', 'lower'],
    c.index.upper = fit.exp.boot.ests['c.test', 'upper'],
    calibration.score = fit.exp.boot.ests['calibration.score', 'val'],
    calibration.score.lower = fit.exp.boot.ests['calibration.score', 'lower'],
    calibration.score.upper = fit.exp.boot.ests['calibration.score', 'upper']
  ),
  performance.file,
  index.cols = c('model', 'imputation', 'discretised')
)

#' ## Performance
#' 
#' Having fitted the Cox model, how did we do? The c-indices were calculated as
#' part of the bootstrapping, so we just need to take a look at those...
#' 
#' C-indices are **`r round(fit.exp.boot.ests['c.train', 'val'], 3)`
#' (`r round(fit.exp.boot.ests['c.train', 'lower'], 3)` - 
#' `r round(fit.exp.boot.ests['c.train', 'upper'], 3)`)** on the training set and
#' **`r round(fit.exp.boot.ests['c.test', 'val'], 3)`
#' (`r round(fit.exp.boot.ests['c.test', 'lower'], 3)` - 
#' `r round(fit.exp.boot.ests['c.test', 'upper'], 3)`)** on the test set.
#' Not too bad!
#' 
#' 
#' ### Calibration
#' 
#' The bootstrapped calibration score is
#' **`r round(fit.exp.boot.ests['calibration.score', 'val'], 3)`
#' (`r round(fit.exp.boot.ests['calibration.score', 'lower'], 3)` - 
#' `r round(fit.exp.boot.ests['calibration.score', 'upper'], 3)`)**.
#' 
#' Let's draw a representative curve from the unbootstrapped fit... (It would be
#' better to draw all the curves from the bootstrap fit to get an idea of
#' variability, but I've not implemented this yet.)
#' 
#+ calibration_plot

calibration.table <-
  calibrationTable(fit.exp, imputed.data[[i]][test.set, ])

calibration.score <- calibrationScore(calibration.table)

calibrationPlot(calibration.table)

#' The area between the calibration curve and the diagonal is 
#' **`r round(calibration.score[['area']], 3)`** +/-
#' **`r round(calibration.score[['se']], 3)`**.
#'  
#' ## Coefficients
#' 
#' As well as getting comparable C-indices, it's also worth checking to see how
#' the risk coefficients calculated compare to those found in the original
#' paper. Let's compare...

# Load CSV of values from paper
old.coefficients <- read.csv(old.coefficients.filename)

# Get coefficients from this fit
new.coefficients <-
  bootMIStats(fit.exp.boot, uncertainty = '95ci', transform = negExp)
names(new.coefficients) <- c('our_value', 'our_lower', 'our_upper')
new.coefficients$quantity.level <- rownames(new.coefficients)

# Create a data frame comparing them
compare.coefficients <- merge(old.coefficients, new.coefficients)

# Kludge because age:genderWomen is the pure interaction term, not the risk for
# a woman per unit of advancing spline-transformed age
compare.coefficients[
  compare.coefficients$quantity.level == 'age:genderWomen', 'our_value'
  ] <-
  compare.coefficients[
    compare.coefficients$quantity.level == 'age:genderWomen', 'our_value'
    ] *
  compare.coefficients[
    compare.coefficients$quantity.level == 'age', 'our_value'
    ]

# Save CSV of results
write.csv(compare.coefficients, output.filename)

# Plot a graph by which to judge success
ggplot(compare.coefficients, aes(x = their_value, y = our_value)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = 1, colour = 'grey') +
  geom_vline(xintercept = 1, colour = 'grey') +
  geom_point() +
  geom_errorbar(aes(ymin = our_lower, ymax = our_upper)) +
  geom_errorbarh(aes(xmin = their_lower, xmax = their_upper)) +
  geom_text_repel(aes(label = long_name)) +
  theme_classic(base_size = 8)

#+ coefficients_table, results='asis'

print(
  xtable(
    data.frame(
      variable =
        paste(
          compare.coefficients$long_name, compare.coefficients$unit, sep=', '
        ),
      compare.coefficients[c('our_value', 'their_value')]
    ),
    digits = c(0,0,3,3)
  ),
  type = 'html',
  include.rownames = FALSE
)

#' ### Variable importance
#' 
#' Let's compare the variable importance from this method with accounting for
#' missing values explicitly. Slight kludge as it's only using one imputed
#' dataset and a fit based on another, but should give some idea.
#' 
#+ cox_variable_importance

cox.var.imp.perm <- 
  generalVarImp(
    fit.exp, imputed.data[[2]][test.set, ], model.type = 'survreg'
  )

write.csv(cox.var.imp.perm, cox.var.imp.perm.filename, row.names = FALSE)

cox.var.imp.perm.missing <- read.csv(cox.var.imp.perm.missing.filename)

cox.var.imp.comparison <-
  merge(
    cox.var.imp.perm,
    cox.var.imp.perm.missing,
    by = 'var',
    suffixes = c('', '.missing')
  )

ggplot(cox.var.imp.comparison, aes(x = var.imp.missing, y = var.imp)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

#' There's a good correlation!
